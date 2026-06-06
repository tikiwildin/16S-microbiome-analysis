library(ape)
library(ggplot2)
library(dplyr)

file_path <- readline("Enter path to braycurtis.0.03.lt: ")

lines <- readLines(file_path)
lines <- trimws(lines)

# Number of samples
n <- as.numeric(lines[1])

# Sample rows
sample_lines <- lines[2:(n + 1)]

# Extract labels = first token from each row
labels <- sapply(strsplit(sample_lines, "\\s+"), `[`, 1)

# Initialize matrix
mat <- matrix(0, n, n)
rownames(mat) <- labels
colnames(mat) <- labels

# Fill lower triangle from each row
for (i in seq_len(n)) {
  parts <- strsplit(sample_lines[i], "\\s+")[[1]]
  values <- as.numeric(parts[-1])
  
  if (length(values) > 0) {
    mat[i, 1:length(values)] <- values
    mat[1:length(values), i] <- values
  }
}

diag(mat) <- 0
dist_obj <- as.dist(mat)

# PCoA
pcoa <- ape::pcoa(dist_obj)

pcoa_df <- data.frame(
  Sample = labels,
  PC1 = pcoa$vectors[, 1],
  PC2 = pcoa$vectors[, 2],
  stringsAsFactors = FALSE
)

# ---- Sample metadata for all 29 samples ----
sample_info <- data.frame(
  Sample = c(
    "1","2","3","4","5","6","7","8","9","10","11","13","14","15","16","17","18",
    "31","32","33","34","35","36","37","38","39","40","41","42"
  ),
  SampleName = c(
    "R5-1St-052819","R5-2St-052819","R5-3St-052819",
    "R5-1A-052819","R5-2A-052819","R5-3A-052819",
    "R2A-1(2)St-052819","R2A-2(2)St-052819","R2A-3(2)St-052819",
    "R2A-1(2)A-052819","R2A-2(2)A-052819",
    "R5-1(2)St-052819","R5-2(2)St-052819","R5-3(2)St-052819",
    "R5-1(2)A-052819","R5-2(2)A-052819","R5-3(2)A-052819",
    "AM2-1","AM2-2","AM2-3","AM3-1","AM3-2","AM3-3",
    "R2A-1St-052419","R2A-2St-052419","R2A-3St-052419",
    "R2A-1A-052419","R2A-2A-052419","R2A-3A-052419"
  ),
  Condition = c(
    rep("R5", 6),
    rep("R2A", 5),
    rep("R5", 6),
    rep("Native", 6),
    rep("R2A", 6)
  ),
  Oxygen = c(
    rep("Static", 3), rep("Aerated", 3),
    rep("Static", 3), rep("Aerated", 2),
    rep("Static", 3), rep("Aerated", 3),
    rep("None", 6),
    rep("Static", 3), rep("Aerated", 3)
  ),
  stringsAsFactors = FALSE
)

# Merge metadata onto ordination coordinates
pcoa_df <- left_join(pcoa_df, sample_info, by = "Sample")

# Keep condition order nice in legend
pcoa_df$Condition <- factor(pcoa_df$Condition, levels = c("Native", "R2A", "R5"))
pcoa_df$Oxygen <- factor(pcoa_df$Oxygen, levels = c("None", "Static", "Aerated"))

# Plot
ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Condition, shape = Oxygen, label = Sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.8, size = 3) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Beta Diversity (Bray-Curtis PCoA)",
    x = paste0("PC1 (", round(100 * pcoa$values$Relative_eig[1], 1), "%)"),
    y = paste0("PC2 (", round(100 * pcoa$values$Relative_eig[2], 1), "%)")
  )