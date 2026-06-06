library(ggplot2)
library(dplyr)
library(readr)
library(ggpubr)

# -------------------------------
# Read input files
# -------------------------------
shannon_file <- readline("Enter path to shannon tabular file: ")
metadata_file <- readline("Enter path to metadata TSV file: ")

shannon_df <- read_tsv(shannon_file, show_col_types = FALSE)
meta_df <- read_tsv(metadata_file, show_col_types = FALSE)

# -------------------------------
# Check columns
# -------------------------------
cat("\nShannon file columns:\n")
print(colnames(shannon_df))

cat("\nMetadata file columns:\n")
print(colnames(meta_df))

# -------------------------------
# Clean Shannon data
# -------------------------------
required_shannon_cols <- c("group", "shannon", "shannon_lci", "shannon_hci")

missing_shannon <- setdiff(required_shannon_cols, colnames(shannon_df))
if (length(missing_shannon) > 0) {
  stop(
    paste(
      "Missing required Shannon columns:",
      paste(missing_shannon, collapse = ", ")
    )
  )
}

shannon_df <- shannon_df %>%
  mutate(
    group = as.character(group),
    shannon = as.numeric(shannon),
    shannon_lci = as.numeric(shannon_lci),
    shannon_hci = as.numeric(shannon_hci)
  )

# -------------------------------
# Clean metadata
# -------------------------------
required_meta_cols <- c("SampleID", "Group", "Medium", "Aeration")

missing_meta <- setdiff(required_meta_cols, colnames(meta_df))
if (length(missing_meta) > 0) {
  stop(
    paste(
      "Missing required metadata columns:",
      paste(missing_meta, collapse = ", ")
    )
  )
}

meta_df <- meta_df %>%
  mutate(
    group = as.character(SampleID),
    Condition = case_when(
      Group == "no_culture" ~ "Native",
      Group == "R2A" ~ "R2A",
      Group == "R5" ~ "R5",
      TRUE ~ NA_character_
    )
  ) %>%
  select(group, Condition, Aeration)

# -------------------------------
# Merge
# -------------------------------
plot_df <- shannon_df %>%
  left_join(meta_df, by = "group")

cat("\nMerged data preview:\n")
print(plot_df)

if (any(is.na(plot_df$Condition))) {
  warning("Some samples did not match metadata and have NA Condition values.")
}

plot_df <- plot_df %>%
  filter(!is.na(Condition))

plot_df$Condition <- factor(plot_df$Condition, levels = c("Native", "R2A", "R5"))

plot_df <- plot_df %>%
  mutate(group_num = suppressWarnings(as.numeric(group))) %>%
  arrange(Condition, group_num)

plot_df$group <- factor(plot_df$group, levels = plot_df$group)

# -------------------------------
# Bar plot with error bars
# -------------------------------
bar_plot <- ggplot(plot_df, aes(x = group, y = shannon, fill = Condition)) +
  geom_col() +
  geom_errorbar(aes(ymin = shannon_lci, ymax = shannon_hci), width = 0.2) +
  scale_fill_manual(values = c(
    "Native" = "darkgreen",
    "R2A" = "steelblue",
    "R5" = "firebrick"
  )) +
  labs(
    title = "Alpha Diversity (Shannon Index)",
    x = "Sample",
    y = "Shannon Diversity"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  )

print(bar_plot)
ggsave("antmound_shannon_barplot.png", bar_plot, width = 11, height = 6, dpi = 300)

# -------------------------------
# Boxplot with jitter + significance
# -------------------------------
box_plot <- ggplot(plot_df, aes(x = Condition, y = shannon, fill = Condition)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.12, size = 2.5, alpha = 0.85) +
  stat_compare_means(
    comparisons = list(
      c("Native", "R2A"),
      c("Native", "R5"),
      c("R2A", "R5")
    ),
    method = "wilcox.test",
    label = "p.signif"
  ) +
  scale_fill_manual(values = c(
    "Native" = "firebrick",
    "R2A" = "darkgreen",
    "R5" = "steelblue"
  )) +
  labs(
    title = "Shannon Diversity Across Conditions",
    x = "Condition",
    y = "Shannon Diversity"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none"
  )

print(box_plot)
ggsave("antmound_shannon_boxplot_stats.png", box_plot, width = 8, height = 6, dpi = 300)

# -------------------------------
# Pairwise stats output
# -------------------------------
cat("\nPairwise Wilcoxon tests:\n")
pairwise_results <- pairwise.wilcox.test(
  x = plot_df$shannon,
  g = plot_df$Condition,
  p.adjust.method = "BH"
)

print(pairwise_results)

# -------------------------------
# Summary stats
# -------------------------------
cat("\nSummary statistics by condition:\n")
summary_stats <- plot_df %>%
  group_by(Condition) %>%
  summarise(
    n = n(),
    mean_shannon = mean(shannon, na.rm = TRUE),
    median_shannon = median(shannon, na.rm = TRUE),
    sd_shannon = sd(shannon, na.rm = TRUE)
  )

print(summary_stats)
