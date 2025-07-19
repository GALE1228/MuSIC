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
})

# Disable default -h (help) short option to avoid conflicts
parser <- OptionParser(add_help_option=FALSE,
                       option_list = list(
                         make_option(c("--help"), action="store_true", default=FALSE,
                                     help="Show help information and exit"),
                         make_option(c("-H", "--heatmap"), type="character", default=NULL,
                                     help="Heatmap input: first column is 'RBP_Transcript.bed', remaining columns are 0/1 values for each species", metavar="file"),
                         make_option(c("-B", "--boxplot"), type="character", default=NULL,
                                     help="Boxplot input: format <path>/<RBP>_<transcript>.bed <val1> <val2> <val3>", metavar="file"),
                         make_option(c("-c", "--col"), type="integer", default=1,
                                     help="Boxplot column selection (1-3) [default %default]", metavar="int"),
                         make_option(c("-o", "--output"), type="character", default="combined_plot.pdf",
                                     help="Output combined plot PDF path [default %default]", metavar="file")
                       ))
opt <- parse_args(parser)

if (opt$help) {
  print_help(parser); q(status=0)
}
if (is.null(opt$heatmap) || is.null(opt$boxplot)) {
  stop("Must specify both -H/--heatmap and -B/--boxplot input files")
}
if (!(opt$col %in% 1:3)) {
  stop("Parameter -c must be 1, 2, or 3")
}

#### 1. Process heatmap data ####
dt_h <- fread(opt$heatmap, header=TRUE, sep="\t", data.table=FALSE)
colnames(dt_h)[1] <- "File"
dt_h$RBP <- str_split_fixed(basename(dt_h$File), "_", 2)[,1]
species_cols <- setdiff(colnames(dt_h), c("File","RBP"))
dt_h[species_cols] <- lapply(dt_h[species_cols], as.integer)

# Count total transcripts for each RBP as denominator
transcript_counts <- as.integer(table(dt_h$RBP))
names(transcript_counts) <- names(table(dt_h$RBP))

# Aggregate 0/1 values by RBP for each species
agg_h <- aggregate(dt_h[, species_cols, drop=FALSE],
                   by = list(RBP = dt_h$RBP), FUN = sum)
mat_h <- as.matrix(agg_h[, species_cols])
rownames(mat_h) <- agg_h$RBP

desired_order <- c("HUMAN","PONAB","MACFA","MOUSE","RAT",
                   "CHICK","XENLA","DANRE","DROME","ARATH","YEAST")
col_order <- intersect(desired_order, colnames(mat_h))
if (length(col_order)==0) {
  stop("Species columns in heatmap input file do not match preset order")
}

# Sort RBPs: by sum of all species from low to high
row_sums_h <- rowSums(mat_h[, col_order, drop=FALSE])
row_order <- names(sort(row_sums_h, decreasing = FALSE))

# Convert to long format and calculate percentage: Count / transcript_count * 100
long_h <- melt(mat_h, varnames=c("RBP","Species"), value.name="Count")
long_h$RBP <- factor(long_h$RBP, levels=row_order)
long_h$Species <- factor(long_h$Species, levels=col_order)
long_h$Percent <- long_h$Count / transcript_counts[as.character(long_h$RBP)] * 100

#### 2. Process boxplot data ####
df_b <- fread(opt$boxplot, header=FALSE, sep="\t", data.table=FALSE,
              col.names=c("path","V1","V2","V3"))
df_b$basename <- basename(df_b$path)
df_b$basename <- str_replace(df_b$basename, "\\.bed$", "")
parts_b <- str_split_fixed(df_b$basename, "_", 2)
df_b$RBP <- parts_b[,1]

df_b$norm1 <- df_b$V1/11*100
df_b$norm2 <- df_b$V2/10*100
df_b$norm3 <- df_b$V3/10*100
col_name <- paste0("norm", opt$col)
plot_b <- data.frame(RBP = df_b$RBP, value = df_b[[col_name]], stringsAsFactors=FALSE)

extra_rb <- setdiff(unique(plot_b$RBP), row_order)
combined_order <- c(row_order, extra_rb)
plot_b$RBP <- factor(plot_b$RBP, levels=combined_order)

#### 3. Build ggplot objects ####
# Calculate mean for each RBP group
group_means <- plot_b %>%
  group_by(RBP) %>%
  summarise(mean_value = mean(value))

# Merge means to original data frame
plot_b <- left_join(plot_b, group_means, by = "RBP")

# Boxplot: RBP on x-axis, fill color based on mean
p_boxplot <- ggplot(plot_b, aes(x = RBP, y = value, fill = mean_value)) +
  geom_boxplot(outlier.size = 0.5, median.size = 0) +
  # Set fill color as gradient
  scale_fill_gradient(low = "#60a9cf", high = "#e8b6ae", 
                      name = "Mean Value", 
                      guide = "colorbar") +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 270, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  ylab(sprintf("Percent of species by method %d (%%)", opt$col))

# Heatmap: RBP on x-axis, fill color using Percent
p_heatmap <- ggplot(long_h, aes(x = RBP, y = Species, fill = Percent)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "#60a9cf", high = "#e8b6ae",
                      name = "Percentage of conservation (%)", labels = function(x) sprintf("%.1f%%", x)) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 270, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )

#### 4. Calculate dimensions and arrange vertically (keep RBP aligned) ####
n_rbp <- length(combined_order)
n_species <- length(col_order)

# Boxplot height: n_rbp*0.03 + 1, minimum 2, maximum 8
box_h   <- min(max(2, n_rbp * 0.03 + 1), 5)
# Heatmap height: n_species*0.3 + 1, minimum 4, maximum 20
heat_h  <- min(max(4, n_species * 0.3 + 1), 20)
total_h <- box_h + heat_h

# Width based on RBP count, each RBP*0.2 + 2, minimum 6, maximum 50
total_w <- min(max(6, n_rbp * 0.2 + 2), 50)

combined_plot <- p_boxplot / p_heatmap +
  plot_layout(heights = c(box_h, heat_h)) &
  theme(plot.margin = margin(5,5,5,5))

ggsave(opt$output, combined_plot, width = total_w, height = total_h, limitsize=TRUE)

# Generate TSV file path
tsv_path <- sub("\\.pdf$", ".tsv", opt$output)

# Output RBP average conservation values to TSV file
write.table(group_means, file = tsv_path, sep = "\t", na = "nan", 
            quote = FALSE, row.names = FALSE)

cat(sprintf("Combined plot with RBP as x-axis saved to %s (width=%.1f\", height=%.1f\")\n",
            opt$output, total_w, total_h))
cat(sprintf("RBP average conservation values saved to %s\n", tsv_path))