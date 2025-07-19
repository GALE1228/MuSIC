#!/usr/bin/env Rscript
suppressMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(stringr)
  library(tidyverse)    # for file handling (grep, list.files)
})

# Parse parameters
option_list <- list(
  make_option(c("-i", "--input"),   type="character", help="Input multiple species files, format: <path>/species_1.txt,<path>/species_2.txt, supports regular expressions", metavar="file"),
  make_option(c("-c", "--col"),     type="integer",   default=1, help="Select column (1-3) [default %default]", metavar="int"),
  make_option(c("-p", "--pattern"), type="character", default=".*", help="RBP name regex, only plot matching ones [default all]", metavar="regex"),
  make_option(c("-o", "--output"),  type="character", default="violin_plot.pdf", help="Output PDF path", metavar="file")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$input)) stop("Must specify -i")
if (!(opt$col %in% 1:3)) stop("Parameter -c must be 1, 2, or 3")

# 1. Use regular expressions to read multiple files
input_files <- list.files(path = dirname(opt$input), pattern = basename(opt$input), full.names = TRUE)
if (length(input_files) == 0) {
  stop("No files match the regular expression")
}

# 2. Read each file and merge
species_list <- lapply(input_files, fread, header=FALSE, sep="\t", data.table=FALSE)

# 3. Set column names
col_names <- c("bed", "col1", "col2", "col3")
df_list <- lapply(species_list, function(file_data) {
  colnames(file_data) <- col_names
  return(file_data)
})

# Merge all species files
df <- do.call(rbind, df_list)

# 4. Normalize processing
normalize_col <- function(df, col_index, num_species) {
  return(df[[paste0("col", col_index)]] / num_species * 100)
}

df$value <- normalize_col(df, opt$col, length(input_files))

# 5. Regex filtering
plot_df <- df[grepl(opt$pattern, df$bed), ]
if (nrow(plot_df) == 0) stop("No RBP matches regex")

# 6. Classify by RBP and draw violin plot for each RBP's each file
# Build data by listing all values for each file and RBP
plot_df$RBP <- factor(plot_df$bed)

# 7. Draw violin plot
p <- ggplot(plot_df, aes(x = RBP, y = value, fill = RBP)) +
  geom_violin(trim = TRUE, show.legend = FALSE) +  # Remove legend
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "darkred") +
  theme_bw() +
  xlab(sprintf("Percent (method %d)", opt$col)) +
  ylab("RBP") +
  theme(axis.text.x = element_text(size = 8),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.border = element_blank(),
        legend.position = "none") +  # Hide legend
  facet_wrap(~bed, scales = "free_y")  # Each RBP one row, each row shows violin plots for multiple files

# Output settings
n_rbp <- length(unique(plot_df$RBP))
height <- min(max(4, n_rbp * 0.3), 50)

# Save results
pdf(opt$output, width = 10, height = height)
print(p)
dev.off()

cat(sprintf("Violin plot saved to %s (width=10\", height=%.1f\")\n", opt$output, height))
