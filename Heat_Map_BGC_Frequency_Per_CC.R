# Load required libraries
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(readr)
library(dplyr)
library(viridis)

# Load data
df <- read_csv("AntiSMASH_results for all samples.csv")

# Ensure 'Sample' is unique rownames
df <- as.data.frame(df)
df$Sample <- make.unique(as.character(df$Sample))
rownames(df) <- df$Sample
df$Sample <- NULL

# Standardize CC labels
df$CC <- as.character(df$CC)
df$CC[is.na(df$CC)] <- "NA"
df$CC <- as.factor(df$CC)

# Separate metadata and BGC count columns
bgc_cols <- setdiff(colnames(df), "CC")

# Aggregate total BGC counts per CC
bgc_summary <- df %>%
  group_by(CC) %>%
  summarise(across(all_of(bgc_cols), sum, na.rm = TRUE)) %>%
  ungroup()

# Convert to matrix: rows = BGC types, columns = CCs
bgc_matrix <- as.data.frame(bgc_summary)
rownames(bgc_matrix) <- bgc_matrix$CC
bgc_matrix$CC <- NULL
#bgc_matrix <- t(as.matrix(bgc_matrix))  # Transpose to get BGCs on y-axis

# Define color scale
max_val <- max(bgc_matrix)
mid_val <- round(max_val / 2)

#col_fun <- colorRamp2(c(0, mid_val, max_val),c("#f7fbff", "#6baed6", "#08306b"))
col_fun <- colorRamp2(c(0, mid_val, max_val), c("#deebf7", "#4292c6", "#000000"))


# Draw the heatmap
ht <- Heatmap(
  bgc_matrix,
  name = "BGC Count",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "Biosynthetic Gene Clusters (BGCs)",
  row_title_side = "right",
  row_title_gp = gpar(fontsize = 25, fontface = "bold"),
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 18),
  row_title = "Clonal Complexes (CCs)",
  column_title_side = "bottom",
  column_title_gp = gpar(fontsize = 25, fontface = "bold"),
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 20),
  border = TRUE,
  rect_gp = gpar(col = "grey80", lwd = 0.5),
  
  # ✅ Add this to control the legend (color key)
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 20, fontface = "bold"),
    labels_gp = gpar(fontsize = 16)
  )
)



# Save to PDF
pdf("antiSMASH_BGCtype_vs_CC_totalcounts_title_spaced_fixed.pdf", width = 14, height = 10)
draw(ht, heatmap_legend_side = "left", annotation_legend_side = "bottom")
dev.off()

# Save to PNG
png("antiSMASH_BGCtype_vs_CC_totalcounts_title_spaced_fixed.png", width = 4000, height = 4000, res = 200)
draw(ht, heatmap_legend_side = "left", annotation_legend_side = "bottom")
dev.off()

