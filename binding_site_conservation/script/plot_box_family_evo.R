#!/usr/bin/env Rscript

# Load libraries
suppressMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(stringr)
  library(ggtree)
  library(ggtreeExtra)  # Note: ggtreeExtra (case sensitive)
  library(cluster)      # for hclust
  library(ggdendro)     # for dendro_data
  library(dplyr)        # for group_by and filter
})

# Parse parameters
option_list <- list(
  make_option(c("-i", "--input"),   type="character", help="Input file, format: <path>/<RBP>_<transcript>.bed <v1> <v2> <v3>", metavar="file"),
  make_option(c("-c", "--col"),     type="integer",   default=1, help="Select column (1-3) [default %default]", metavar="int"),
  make_option(c("-p", "--pattern"), type="character", default=".*", help="RBP name regex, only plot matching ones [default all]", metavar="regex"),
  make_option(c("-o", "--output"),  type="character", default="boxplot_cluster_ggtree.pdf", help="Output PDF path", metavar="file")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$input)) stop("Must specify -i")
if (!(opt$col %in% 1:3)) stop("Parameter -c must be 1, 2, or 3")

# 1. Read data and extract RBP
df <- fread(opt$input, header=FALSE, sep="\t",
            col.names=c("path","V1","V2","V3"), data.table=FALSE)

# Extract RBP name from filename
df$basename <- basename(df$path)
df$basename <- str_replace(df$basename, "\\.bed$", "")
parts <- str_split_fixed(df$basename, "_", 2)
df$RBP <- parts[,1]

# 2. Normalize
norm_funcs <- list(
  function(x) x / 11 * 100,
  function(x) x / 10 * 100,
  function(x) x / 10 * 100
)
df$value <- norm_funcs[[opt$col]](df[[paste0("V", opt$col)]])
plot_df <- df[, c("RBP", "value")]

# 3. Regex filtering
plot_df <- plot_df[grepl(opt$pattern, plot_df$RBP), ]
if (nrow(plot_df) == 0) stop("No RBP matches regex")

# 4. Aggregate means for clustering
mean_df <- aggregate(value ~ RBP, data = plot_df, FUN = mean)

# If only one RBP, draw simple boxplot
if (nrow(mean_df) < 2) {
  p <- ggplot(plot_df, aes(x = value, y = RBP)) +
    geom_boxplot(width = 0.5, fill = "steelblue", outlier.size = 0.5) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "darkred") +
    theme_bw() +
    xlab(sprintf("Percent (method %d)", opt$col)) +
    ylab("") +
    theme(axis.text.y = element_text(size = 8),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.border = element_blank())
  
  ggsave(opt$output, p, width = 6, height = 1.5 + 0.2 * nrow(mean_df), limitsize = TRUE)
  cat("Only 1 RBP, no clustering tree, simple boxplot output\n")
  quit(status = 0)
}

# 5. Build clustering tree
mat <- as.matrix(mean_df$value)
rownames(mat) <- mean_df$RBP
dist_mat <- dist(mat)
hclust_tree <- hclust(dist_mat, method = "average")

# Convert to dendrogram type for ggtree
dend <- as.dendrogram(hclust_tree)

# 6. Construct data frame for ggtree
tree_data <- ggdendro::dendro_data(dend, type = "rectangle")
label_data <- tree_data$label

# 7. Add boxplot as "fruit" to leaf nodes (vertical)
# First filter out RBPs with only one unique value or zero variance
filtered_plot_df <- plot_df %>%
  group_by(RBP) %>%
  filter(n_distinct(value) > 1 & var(value, na.rm = TRUE) > 0) %>%  # ✅ Removed !any(value == 0)
  ungroup() %>%
  na.omit()

# ✅ Debug output: check filtered data
cat("Number of RBPs after filtering:", length(unique(filtered_plot_df$RBP)), "\n")
print(head(filtered_plot_df))

# ✅ New: if no data after filtering, exit and prompt
if (nrow(filtered_plot_df) == 0) {
  cat("❌ No valid data after filtering, please check if data is all 0 or has only one unique value\n")
  quit(save = "no", status = 1)
}

# Add boxplot
p <- ggtree(tree_data, layout = "rectangular") +
  geom_tiplab(aes(label = label), offset = 0.05, align = TRUE, linesize = 0.3) +
  geom_fruit(
    data = filtered_plot_df,
    geom = geom_boxplot,
    mapping = aes(x = RBP, y = value),
    width = 0.6,
    fill = "steelblue",
    outlier.size = 0.5
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_tree() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Add mean markers (vertical)
p <- p + geom_fruit(
  data = mean_df,
  geom = geom_crossbar,
  mapping = aes(ymin = value, ymax = value, x = RBP),
  color = "darkred",
  width = 0.4
)

# Output settings
n_rbp <- nrow(mean_df)
height <- min(max(4, n_rbp * 0.5), 50)

# Save results
pdf(opt$output, width = 10, height = height)
print(p)
dev.off()

cat(sprintf("Clustered horizontal boxplot saved to %s (width=10\", height=%.1f\")\n", opt$output, height))