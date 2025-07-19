#!/usr/bin/env Rscript
suppressMessages({
  library(data.table)
  library(ggplot2)
  library(optparse)
})

# Define command line arguments
parser <- OptionParser(add_help_option = FALSE,
                       option_list = list(
                         make_option(c("-h", "--help"), action = "store_true", default = FALSE,
                                     help = "Display help information and exit"),
                         make_option(c("-e", "--edgelist"), type = "character", default = "/home/bioinfo/02_project/05_RBP/05_homo/similarity/merge/example_allRBP_allRNA__edgelist.csv",
                                     help = "Edgelist file path [default %default]", metavar = "file"),
                         make_option(c("-t", "--tsv"), type = "character", default = "/home/bioinfo/02_project/05_RBP/05_homo/figs/count_side_v2.tsv",
                                     help = "TSV file path [default %default]", metavar = "file"),
                         make_option(c("-o", "--output"), type = "character", default = "/home/bioinfo/02_project/05_RBP/05_homo/figs/correlation_plot.pdf",
                                     help = "Output PDF path for the plot [default %default]", metavar = "file")
                       ))

opt <- parse_args(parser)

# If user requests help, display help information and exit
if (opt$help) {
  print_help(parser)
  quit(status = 0)
}

# Read data
edgelist <- fread(opt$edgelist)
tsv_data <- fread(opt$tsv)

# Calculate importance of each RBP in the network (sum of weights)
rbp_importance_from <- edgelist[, .(importance = sum(adjusted_weight)), by = from]
setnames(rbp_importance_from, "from", "RBP")
rbp_importance_to <- edgelist[, .(importance = sum(adjusted_weight)), by = to]
setnames(rbp_importance_to, "to", "RBP")
rbp_importance <- rbind(rbp_importance_from, rbp_importance_to)[, .(importance = sum(importance)), by = RBP]

# Merge data
merged_data <- merge(rbp_importance, tsv_data, by = "RBP", all = TRUE)

# Calculate correlation
cor_result <- cor.test(merged_data$importance, merged_data$mean_value)
r_value <- round(cor_result$estimate, 2)
p_value <- format.pval(cor_result$p.value, digits = 2)

# Calculate number of points N
N <- nrow(merged_data)

# Prepare correlation annotation text, include N, and label as Pearson's
cor_label <- paste0("Pearson's r = ", r_value, "\np = ", p_value, "\nN = ", N)

# Draw correlation analysis plot, add regression line
p <- ggplot(merged_data, aes(x = importance, y = mean_value)) +
  geom_point() +
  # Add linear regression line, use lm method, don't show confidence interval
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  # Explicitly set x, y axis English names, reduce font size
  labs(x = "Weighted Degree Centrality of RBP", y = "Cross-species Conservation of RBP Binding Sites") +
  theme_minimal(base_size = 8) +
  # Set x-axis ticks
  scale_x_continuous(breaks = c(0, 10, 20)) +
  # Set y-axis ticks
  scale_y_continuous(breaks = c(0, 50, 100)) +
  # Add correlation statistics annotation to bottom right
  annotate("text", x = max(merged_data$importance, na.rm = TRUE), 
           y = min(merged_data$mean_value, na.rm = TRUE), 
           label = cor_label, hjust = 1, vjust = 0) +
  # Remove grid lines, keep xy axes
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )

# Save plot
ggsave(opt$output, p, width = 4, height = 3)

cat(sprintf("Correlation plot has been saved to %s\n", opt$output))