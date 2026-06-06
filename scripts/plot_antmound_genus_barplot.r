# Script: plot_antmound_genus_barplot.R
# Project: 16S rRNA Microbiome Analysis
#
# Purpose:
#   This script reads mothur/Krona genus-level taxonomy files from
#   ant mound microbiome samples, calculates genus-level relative
#   abundance, and generates a stacked barplot showing dominant genera
#   across native, R2A-cultured, and R5-cultured sample groups.
# Input:
#   - krona_files/
#       Directory containing mothur/Krona tabular taxonomy files
#       named in the format 0.03.<SampleID>.tabular
#   - antmound_metadata.tsv
#       Metadata table containing sample IDs and experimental groupings
#
# Output:
#   - antmound_genus_stacked_barplot.png
#   - antmound_genus_relative_abundance.tsv
# Author: Atiq Bacus

library(tidyverse)

# -------------------------
# Paths
# -------------------------
taxonomy_dir <- "krona_files"          # folder containing 0.03.*.tabular
metadata_file <- "antmound_metadata.tsv"

# -------------------------
# Read metadata
# -------------------------
meta <- read.delim(metadata_file, check.names = FALSE) %>%
  mutate(SampleID = as.character(SampleID))

# -------------------------
# List files
# -------------------------
files <- list.files(
  taxonomy_dir,
  pattern = "^0\\.03\\.[0-9]+\\.tabular$",
  full.names = TRUE
)

# -------------------------
# Function to read one file
# -------------------------
read_krona_file <- function(file) {
  sample_id <- str_match(basename(file), "^0\\.03\\.([0-9]+)\\.tabular$")[,2]
  
  df <- read.delim(
    file,
    header = FALSE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # first column = count
  # last column = genus
  count_col <- df[[1]]
  genus_col <- df[[ncol(df)]]
  
  tibble(
    SampleID = sample_id,
    Count = as.numeric(count_col),
    GenusRaw = genus_col
  ) %>%
    mutate(
      # remove confidence scores like "(100)"
      Genus = str_remove(GenusRaw, "\\([0-9]+\\)$"),
      Genus = if_else(is.na(Genus) | Genus == "" | str_detect(Genus, "unclassified|_unclassified"),
                      "Unclassified",
                      Genus)
    ) %>%
    group_by(SampleID, Genus) %>%
    summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop")
}

# -------------------------
# Read all files
# -------------------------
genus_counts <- map_dfr(files, read_krona_file)

# -------------------------
# Convert to relative abundance
# -------------------------
genus_rel <- genus_counts %>%
  group_by(SampleID) %>%
  mutate(RelAbundance = Count / sum(Count)) %>%
  ungroup()

# -------------------------
# Keep top genera overall
# -------------------------
top_n <- 12

top_genera <- genus_rel %>%
  group_by(Genus) %>%
  summarise(TotalAbundance = sum(RelAbundance), .groups = "drop") %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = top_n) %>%
  pull(Genus)

genus_plot <- genus_rel %>%
  mutate(GenusPlot = if_else(Genus %in% top_genera, Genus, "Other")) %>%
  group_by(SampleID, GenusPlot) %>%
  summarise(RelAbundance = sum(RelAbundance), .groups = "drop")

# -------------------------
# Merge metadata
# -------------------------
plot_df <- genus_plot %>%
  left_join(meta, by = "SampleID")

# order samples by metadata grouping
plot_df <- plot_df %>%
  mutate(
    Group = factor(Group, levels = c("no_culture", "R2A", "R5")),
    Aeration = factor(Aeration, levels = c("none", "static", "aerated")),
    SampleID = factor(SampleID, levels = meta$SampleID)
  )

# -------------------------
# Plot
# -------------------------
p <- ggplot(plot_df, aes(x = SampleID, y = RelAbundance, fill = GenusPlot)) +
  geom_col() +
  facet_grid(. ~ Group, scales = "free_x", space = "free_x") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Dominant Genera Across Ant Mound Samples",
    x = "Sample",
    y = "Relative Abundance",
    fill = "Genus"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave("antmound_genus_stacked_barplot.png", p, width = 14, height = 7, dpi = 300)

# -------------------------
# Save processed table
# -------------------------
write.table(
  plot_df,
  file = "antmound_genus_relative_abundance.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


