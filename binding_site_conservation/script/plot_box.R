#!/usr/bin/env Rscript
suppressMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(stringr)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input file, format: <path>/<RBP>_<transcript>.bed <val1> <val2> <val3>", metavar="file"),
  make_option(c("-c", "--col"), type="integer", default=1,
              help="Select column for plotting (1-3) [default %default]", metavar="int"),
  make_option(c("-o", "--output"), type="character", default="boxplot.pdf",
              help="Output PDF file path [default %default]", metavar="file")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input)) {
  stop("Must specify input file path -i")
}
if (!(opt$col %in% 1:3)) {
  stop("Parameter -c must be 1, 2, or 3")
}

# Read data
df <- fread(opt$input, header=FALSE, sep="\t", data.table=FALSE,
            col.names = c("path", "V1", "V2", "V3"))

# Extract RBP and transcript
df$basename <- basename(df$path)
df$basename <- str_replace(df$basename, "\\.bed$", "")
parts <- str_split_fixed(df$basename, "_", 2)
df$RBP <- parts[,1]
df$transcript <- parts[,2]

# Normalize
df$norm1 <- df$V1 / 11 * 100
df$norm2 <- df$V2 / 10 * 100
df$norm3 <- df$V3 / 10 * 100

col_name <- paste0("norm", opt$col)
plot_df <- data.frame(RBP = df$RBP, value = df[[col_name]], stringsAsFactors = FALSE)

# Sort by mean for each RBP
mean_by_rbp <- aggregate(value ~ RBP, data = plot_df, FUN = mean)
ordered_rbps <- mean_by_rbp[order(mean_by_rbp$value), "RBP"]

# Set factor levels
plot_df$RBP <- factor(plot_df$RBP, levels = ordered_rbps)

# Adjust width by number of RBPs: about 0.2 inch per category, minimum 6 inches
n_rbp <- length(ordered_rbps)
plot_width <- max(6, n_rbp * 0.2)
plot_height <- 4

# Draw boxplot, use mean line instead of median
p <- ggplot(plot_df, aes(x = RBP, y = value)) +
  geom_boxplot(fill = "steelblue", outlier.size = 0.8, median.size = 0) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.6, color = "darkred", fatten = 0) +
  theme_bw() +
  xlab("RBP") +
  ylab(sprintf("Percent of species by method %d (%%)", opt$col)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save to PDF (limit size to max 50 inches)
ggsave(opt$output, plot = p, width = plot_width, height = plot_height, limitsize = TRUE)

cat(sprintf("Boxplot with mean line saved to %s\n", opt$output))
