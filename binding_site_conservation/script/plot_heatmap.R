#!/usr/bin/env Rscript
suppressMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(stringr)
  library(reshape2)
  library(scales)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input statistics result file, first column is 'RBP_Transcript.bed', remaining columns are 0/1 values for each species", metavar="file"),
  make_option(c("-o", "--output"), type="character", default="heatmap.pdf",
              help="Output heatmap PDF path [default %default]", metavar="file")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input)) {
  stop("Must specify input file path -i")
}

# 1. Read input, keep column names
dt <- fread(opt$input, header=TRUE, sep="\t", data.table=FALSE)

# 2. Extract RBP name (remove underscore and everything after)
colnames(dt)[1] <- "File"
dt$RBP <- str_split_fixed(basename(dt$File), "_", 2)[,1]

# 3. Extract species columns (except File and RBP), convert to integers
species_cols <- setdiff(colnames(dt), c("File","RBP"))
dt[species_cols] <- lapply(dt[species_cols], as.integer)

# 4. Aggregate by RBP: sum each species column, generate RBP × Species matrix
agg_dt <- aggregate(dt[, species_cols, drop=FALSE],
                    by = list(RBP = dt$RBP),
                    FUN = sum)
mat <- as.matrix(agg_dt[, species_cols])
rownames(mat) <- agg_dt$RBP

# 5. Fix species order (HUMAN first, rest by established order)
desired_order <- c("HUMAN", "PONAB", "MACFA", "MOUSE", "RAT",
                   "CHICK", "XENLA", "DANRE", "DROME", "ARATH", "YEAST")
col_order <- intersect(desired_order, colnames(mat))
if (length(col_order) == 0) {
  stop("Species columns in input file don't match predefined order, please check column names")
}

# 6. Sort RBPs by sum of all species from low to high (low first, high last)
row_sums <- rowSums(mat[, col_order, drop=FALSE])
row_order <- names(sort(row_sums, decreasing = FALSE))

# 7. Convert to long format: columns "RBP", "Species", "Count", and set factors in order
long_df <- melt(mat, varnames = c("RBP", "Species"), value.name = "Count")
long_df$RBP <- factor(long_df$RBP, levels = row_order)
long_df$Species <- factor(long_df$Species, levels = col_order)

# 8. Dynamically adjust plot size and limit to 50 inches
n_rbp <- length(row_order)
n_species <- length(col_order)
# More RBP rows, use narrower height: about 0.15 inch per row, minimum 4 inches, maximum 50
height_in <- min(max(4, n_rbp * 0.15 + 1), 50)
# Fewer species columns, about 0.6 inch per column for better display, minimum 4 inches, maximum 50
width_in  <- min(max(4, n_species * 0.6 + 1), 50)

# 9. Draw heatmap, put species names on top, reduce font size
p <- ggplot(long_df, aes(x = Species, y = RBP, fill = Count)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue",
                      name = "Overlap Count", labels = comma) +
  scale_x_discrete(position = "top") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 6),
    axis.text.x.bottom = element_blank(),
    axis.text.y = element_text(size = 5),
    axis.title   = element_blank(),
    panel.grid   = element_blank()
  )

# 10. Save to PDF
ggsave(opt$output, plot = p, width = width_in, height = height_in, limitsize = TRUE)

cat(sprintf("Heatmap sorted by total from low to high with adjusted species names and size saved to %s (width=%.1f\", height=%.1f\")\n",
            opt$output, width_in, height_in))
