if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")


# Load required libraries
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(readr)
library(viridis)
library(grid)

# Load data
df <- read_csv("*.csv")

# Ensure rownames are unique and from Sample column
df <- as.data.frame(df)
df$Sample <- make.unique(as.character(df$Sample))
rownames(df) <- df$Sample
df$Sample <- NULL

# Handle CC: convert NA to "NA" string label and treat as factor
df$CC <- as.character(df$CC)
df$CC[is.na(df$CC)] <- "NA"
df$CC <- as.factor(df$CC)

# Extract CC labels
cc_labels <- df$CC
df$CC <- NULL  # remove CC column for matrix conversion

# Convert rest of data to numeric matrix (BGC counts per sample)
mat <- as.matrix(df)

# Color scale
max_val <- ceiling(max(mat, na.rm = TRUE))
mid_val <- round(max_val / 2)
#col_fun <- colorRamp2(c(0, mid_val, max_val), c("#f7fbff", "#6baed6", "#08306b"))
col_fun <- colorRamp2(c(0, mid_val, max_val), c("#deebf7", "#4292c6", "#000000"))

# CC annotation colors
cc_levels <- levels(as.factor(cc_labels))
cc_palette <- brewer.pal(n = max(3, length(cc_levels)), "Set3")
names(cc_palette) <- cc_levels

# Row annotation for Clonal Complexes with legend font settings
row_anno <- rowAnnotation(
  CC = cc_labels,
  col = list(CC = cc_palette),
  show_annotation_name = FALSE,
  
  # ✅ Add legend font control here
  annotation_legend_param = list(
    title_gp = gpar(fontsize = 16, fontface = "bold"),
    labels_gp = gpar(fontsize = 13)
  )
)

# Build the heatmap
ht <- Heatmap(
  mat,
  name = "Primary metabolism GCs Count per sample",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 16),
  column_names_gp = gpar(fontsize = 16),
  left_annotation = row_anno,
  column_title = "Primary metabolism GCs",
  row_title_gp = gpar(fontsize = 23, fontface = "bold"),
  row_title = "Samples",
  column_title_gp = gpar(fontsize = 23, fontface = "bold"),
  row_title_side = "right",
  column_title_side = "bottom",
  column_names_rot = 28,
  border = TRUE,
  rect_gp = gpar(col = "grey80", lwd = 0.5),
  
  # ✅ Heatmap legend font control
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 18, fontface = "bold"),
    labels_gp = gpar(fontsize = 14)
  )
)

# Save to PDF
pdf("gutSMASH_heatmap_BGC_Frequency_per_samplel.pdf", width = 14, height = 14)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "bottom")
dev.off()

# Save to PNG
png("gutSMASH_heatmap_BGC_Frequency_per_samplel.png", width = 5000, height = 4000, res = 200)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "bottom")
dev.off()


