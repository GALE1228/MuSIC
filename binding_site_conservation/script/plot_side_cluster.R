#!/usr/bin/env Rscript
suppressMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(stringr)
  library(dplyr)
  library(reshape2)
  library(scales)
  library(patchwork)
  library(pheatmap)  # Add pheatmap for heatmap hierarchical clustering
})

# Disable default -h (help) short option to avoid conflicts
parser <- OptionParser(add_help_option=FALSE,
                       option_list = list(
                         make_option(c("--help"), action="store_true", default=FALSE,
                                     help="Display help information then exit"),
                         make_option(c("-H", "--heatmap"), type="character", default=NULL,
                                     help="Heatmap input: first column is 'RBP_Transcript.bed', remaining columns are 0/1 values for each species", metavar="file"),
                         make_option(c("-o", "--output"), type="character", default="heatmap_cluster_plot.pdf",
                                     help="Output heatmap PDF path [default %default]", metavar="file")
                       ))
opt <- parse_args(parser)

if (opt$help) {
  print_help(parser); q(status=0)
}
if (is.null(opt$heatmap)) {
  stop("Must specify -H/--heatmap input file")
}

#### 1. Process heatmap data ####
dt_h <- fread(opt$heatmap, header=TRUE, sep="\t", data.table=FALSE)
colnames(dt_h)[1] <- "File"
dt_h$RBP <- str_split_fixed(basename(dt_h$File), "_", 2)[,1]
species_cols <- setdiff(colnames(dt_h), c("File","RBP"))
dt_h[species_cols] <- lapply(dt_h[species_cols], as.integer)

# Aggregate 0/1 values by RBP for each species
agg_h <- aggregate(dt_h[, species_cols, drop=FALSE],
                   by = list(RBP = dt_h$RBP), FUN = sum)
mat_h <- as.matrix(agg_h[, species_cols])
rownames(mat_h) <- agg_h$RBP

desired_order <- c("HUMAN","PONAB","MACFA","MOUSE","RAT",
                   "CHICK","XENLA","DANRE","DROME","ARATH","YEAST")
col_order <- intersect(desired_order, colnames(mat_h))
if (length(col_order)==0) {
  stop("Species columns in heatmap input file don't match preset order")
}

# Sort RBPs: by sum of all species from low to high
row_sums_h <- rowSums(mat_h[, col_order, drop=FALSE])
row_order <- names(sort(row_sums_h, decreasing = FALSE))

# Convert to long format and calculate percentage: Count / transcript_count * 100
long_h <- melt(mat_h, varnames=c("RBP","Species"), value.name="Count")
long_h$RBP <- factor(long_h$RBP, levels=row_order)
long_h$Species <- factor(long_h$Species, levels=col_order)
long_h$Percent <- long_h$Count / row_sums_h[as.character(long_h$RBP)] * 100

#### 2. Hierarchical clustering ####
# Use pheatmap package to draw heatmap and perform hierarchical clustering in the process
pheatmap(mat_h, 
         cluster_rows = TRUE,   # Hierarchical clustering of RBP rows
         color = colorRampPalette(c("#60a9cf", "#e8b6ae"))(100),  # Set colors
         main = "RBP Heatmap with Hierarchical Clustering",
         show_rownames = TRUE,  # Show row names (RBP)
         show_colnames = TRUE,  # Show column names (species)
         fontsize = 8,
         cellwidth = 20,        # Cell width
         cellheight = 10,       # Cell height
         border_color = "black", # Border color
         filename = opt$output,  # Output file path
         cutree_rows = 3,        # Cut rows into 3 classes
)

cat(sprintf("Heatmap saved to %s\n", opt$output))
