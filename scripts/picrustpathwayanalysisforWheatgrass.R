# =========================
# PICRUSt2 pathway analysis
# Wheatgrass dataset
# Chronological sample order
# =========================

library(tidyverse)
library(pheatmap)
library(RColorBrewer)

# -------------------------
# File paths
# -------------------------
pathway_file <- "path_abun_unstrat.tsv.gz"
metadata_file <- "wheatgrass_metadata.tsv"

# -------------------------
# Load pathway table
# -------------------------
pathways <- read.delim(pathway_file, check.names = FALSE)
colnames(pathways)[1] <- "Pathway"

# -------------------------
# Transpose so samples are rows
# -------------------------
pathway_mat <- pathways %>%
  column_to_rownames("Pathway") %>%
  as.matrix()

pathway_t <- t(pathway_mat) %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  mutate(SampleID = as.character(SampleID))

# -------------------------
# Load metadata
# -------------------------
meta <- read.delim(metadata_file, check.names = FALSE) %>%
  mutate(
    SampleID = as.character(SampleID),
    Condition = as.factor(Condition),
    Time = as.factor(Time)
  ) %>%
  arrange(as.numeric(SampleID))

# -------------------------
# Merge metadata with pathway table
# -------------------------
dat <- left_join(meta, pathway_t, by = "SampleID")

# Check for unmatched samples
if (any(is.na(dat[, 4]))) {
  warning("Some metadata rows may not have matched sample names in the pathway table.")
}

# -------------------------
# Extract pathway-only matrix
# -------------------------
meta_cols <- c("SampleID", "Condition", "Time")

pathway_only <- dat %>%
  select(-any_of(meta_cols)) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.data.frame()

rownames(pathway_only) <- dat$SampleID

# -------------------------
# Normalize and log-transform
# -------------------------
pathway_rel <- sweep(as.matrix(pathway_only), 1, rowSums(pathway_only), "/")
pathway_rel[is.na(pathway_rel)] <- 0

pathway_log <- log10(pathway_rel + 1e-6)

# -------------------------
# PCA
# -------------------------
pca_res <- prcomp(pathway_log, center = TRUE, scale. = TRUE)

pca_df <- as.data.frame(pca_res$x[, 1:2]) %>%
  rownames_to_column("SampleID") %>%
  left_join(meta, by = "SampleID")

p <- ggplot(pca_df, aes(PC1, PC2, color = Condition, shape = Time)) +
  geom_point(size = 4) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    legend.position = "right"
  ) +
  labs(
    title = "PCA of Wheatgrass PICRUSt2 Pathway Profiles",
    x = paste0("PC1 (", round(100 * summary(pca_res)$importance[2, 1], 1), "%)"),
    y = paste0("PC2 (", round(100 * summary(pca_res)$importance[2, 2], 1), "%)")
  )

print(p)

ggsave(
  "picrust2_wheatgrass_pca.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)

# -------------------------
# Heatmap of most variable pathways
# -------------------------
pathway_vars <- apply(pathway_log, 2, var)

top_n <- 30
top_pathways <- names(sort(pathway_vars, decreasing = TRUE))[1:top_n]

heatmap_mat <- pathway_log[, top_pathways]

# Z-score by pathway
heatmap_scaled <- scale(heatmap_mat)
heatmap_scaled <- t(heatmap_scaled)

# -------------------------
# Force chronological sample order
# -------------------------
ordered_samples <- meta %>%
  arrange(as.numeric(SampleID)) %>%
  pull(SampleID)

# Keep only samples found in the heatmap matrix
ordered_samples <- ordered_samples[ordered_samples %in% colnames(heatmap_scaled)]

heatmap_scaled <- heatmap_scaled[, ordered_samples]

annotation_df <- meta %>%
  filter(SampleID %in% ordered_samples) %>%
  arrange(as.numeric(SampleID)) %>%
  column_to_rownames("SampleID") %>%
  select(Condition, Time)

annotation_df <- annotation_df[ordered_samples, ]

# -------------------------
# Chronological heatmap
# -------------------------
pheatmap(
  heatmap_scaled,
  annotation_col = annotation_df,
  fontsize_row = 8,
  fontsize_col = 10,
  main = "Top Variable Wheatgrass PICRUSt2 Pathways Ordered by SampleID",
  clustering_method = "complete",
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  gaps_col = c(2, 5),
  show_colnames = TRUE,
  show_rownames = TRUE,
  filename = "picrust2_wheatgrass_heatmap_chronological.png",
  width = 10,
  height = 12
)

# -------------------------
# Save processed tables
# -------------------------
write.table(
  pathway_rel,
  file = "picrust2_wheatgrass_pathway_relative_abundance.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

write.table(
  pca_df,
  file = "picrust2_wheatgrass_pca_coordinates.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

