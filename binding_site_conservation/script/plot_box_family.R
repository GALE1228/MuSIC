#!/usr/bin/env Rscript
suppressMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(stringr)
  library(ggdendro)
  library(ggtree)
  library(cowplot)
})

option_list <- list(
  make_option(c("-i", "--input"),   type="character", help="Input file, format: <path>/<RBP>_<transcript>.bed <v1> <v2> <v3>", metavar="file"),
  make_option(c("-c", "--col"),     type="integer",   default=1, help="Select column (1-3) [default %default]", metavar="int"),
  make_option(c("-p", "--pattern"), type="character", default=".*", help="RBP name regex, only plot matching ones [default all]", metavar="regex"),
  make_option(c("-o", "--output"),  type="character", default="boxplot_cluster.pdf", help="Output PDF path", metavar="file")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input)) stop("Must specify -i")
if (!(opt$col %in% 1:3))  stop("Parameter -c must be 1, 2, or 3")

# 1. Read data and extract RBP
df <- fread(opt$input, header=FALSE, sep="\t",
            col.names=c("path","V1","V2","V3"), data.table=FALSE)
df$basename <- basename(df$path)
df$basename <- str_replace(df$basename, "\\.bed$", "")
parts      <- str_split_fixed(df$basename, "_", 2)
df$RBP     <- parts[,1]

# 2. Normalize
df$norm1 <- df$V1/11*100
df$norm2 <- df$V2/10*100
df$norm3 <- df$V3/10*100
col_name  <- paste0("norm", opt$col)
plot_df   <- data.frame(RBP=df$RBP, value=df[[col_name]], stringsAsFactors=FALSE)

# 3. Regex filtering
plot_df <- plot_df[grepl(opt$pattern, plot_df$RBP), ]
if (nrow(plot_df)==0) stop("No RBP matches regex")

# 4. Aggregate means, prepare for clustering
mean_df <- aggregate(value~RBP, data=plot_df, FUN=mean)
rbs    <- mean_df$RBP

# If only one RBP, output horizontal boxplot directly
if (length(rbs)<2) {
  p <- ggplot(plot_df, aes(x=value, y=RBP)) +
    geom_boxplot(width=0.5, fill="steelblue", outlier.size=0.5) +
    stat_summary(fun=mean, geom="crossbar", width=0.5, color="darkred", fatten=0) +
    theme_bw(base_size=12) +
    xlab(sprintf("Percent (method %d)", opt$col)) +
    ylab("") +
    theme(axis.text.y=element_text(size=8),
          panel.grid.major.y=element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.border = element_blank())  # Remove plot border
  ggsave(opt$output, p, width=6, height=1.5 + 0.2*length(rbs), limitsize=TRUE)
  cat("Only 1 RBP, no clustering tree, simple boxplot output\n")
  quit(status=0)
}

# 5. Hierarchical clustering
mat    <- as.matrix(mean_df[,"value",drop=FALSE])
rownames(mat) <- mean_df$RBP
dxy    <- dist(mat)
hc     <- hclust(dxy, method="average")
order_rb <- hc$labels[hc$order]

# 6. Factor order
plot_df$RBP <- factor(plot_df$RBP, levels=order_rb)

# 7. Build boxplot, no outer border
p_box <- ggplot(plot_df, aes(x=value, y=RBP)) +
  geom_boxplot(width=0.6, fill="steelblue", outlier.size=0.5) +
  stat_summary(fun=mean, geom="crossbar", width=0.6, color="darkred", fatten=0) +
  theme_bw(base_size=12) +
  xlab(sprintf("Percent (method %d)", opt$col)) +
  ylab("") +
  theme(axis.text.y=element_text(size=6),  # Show RBP names
        axis.ticks.y=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.border = element_blank())  # Remove plot border

# 8. Use ggtree to draw dendrogram
ddata <- dendro_data(hc, type="rectangle")
p_den <- ggtree(hc, layout = "rectangular") +
  scale_x_continuous(expand = c(0, 0)) +
  theme_tree() +
  theme(legend.position="none")

# 9. Combine and save
n_rbp <- length(order_rb)
h <- min(max(4, n_rbp*0.2 + 1), 50)
combined <- plot_grid(p_den, p_box, nrow=1, rel_widths=c(1,4))
ggsave(opt$output, combined, width=8, height=h, limitsize=TRUE)

cat(sprintf("Clustered horizontal boxplot saved to %s (width=8\", height=%.1f\")\n", opt$output, h))
