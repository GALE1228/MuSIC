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
  make_option(c("-i", "--input"),   type="character", help="输入文件，格式：<path>/<RBP>_<transcript>.bed <v1> <v2> <v3>", metavar="file"),
  make_option(c("-c", "--col"),     type="integer",   default=1, help="选择列 (1-3) [default %default]", metavar="int"),
  make_option(c("-p", "--pattern"), type="character", default=".*", help="RBP 名正则，只画匹配的 [default all]", metavar="regex"),
  make_option(c("-o", "--output"),  type="character", default="boxplot_cluster.pdf", help="输出 PDF 路径", metavar="file")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input)) stop("必须指定 -i")
if (!(opt$col %in% 1:3))  stop("参数 -c 必须是 1、2 或 3")

# 1. 读取数据并提取 RBP
df <- fread(opt$input, header=FALSE, sep="\t",
            col.names=c("path","V1","V2","V3"), data.table=FALSE)
df$basename <- basename(df$path)
df$basename <- str_replace(df$basename, "\\.bed$", "")
parts      <- str_split_fixed(df$basename, "_", 2)
df$RBP     <- parts[,1]

# 2. 归一化
df$norm1 <- df$V1/11*100
df$norm2 <- df$V2/10*100
df$norm3 <- df$V3/10*100
col_name  <- paste0("norm", opt$col)
plot_df   <- data.frame(RBP=df$RBP, value=df[[col_name]], stringsAsFactors=FALSE)

# 3. 正则过滤
plot_df <- plot_df[grepl(opt$pattern, plot_df$RBP), ]
if (nrow(plot_df)==0) stop("没有 RBP 符合正则")

# 4. 汇总均值，准备聚类
mean_df <- aggregate(value~RBP, data=plot_df, FUN=mean)
rbs    <- mean_df$RBP

# 如果只有一个 RBP，直接输出水平箱线图
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
          panel.border = element_blank())  # 移除图边框
  ggsave(opt$output, p, width=6, height=1.5 + 0.2*length(rbs), limitsize=TRUE)
  cat("仅 1 个 RBP，无聚类树，已输出简单箱线图\n")
  quit(status=0)
}

# 5. 层次聚类
mat    <- as.matrix(mean_df[,"value",drop=FALSE])
rownames(mat) <- mean_df$RBP
dxy    <- dist(mat)
hc     <- hclust(dxy, method="average")
order_rb <- hc$labels[hc$order]

# 6. 因子顺序
plot_df$RBP <- factor(plot_df$RBP, levels=order_rb)

# 7. 构建箱线图，不加外框
p_box <- ggplot(plot_df, aes(x=value, y=RBP)) +
  geom_boxplot(width=0.6, fill="steelblue", outlier.size=0.5) +
  stat_summary(fun=mean, geom="crossbar", width=0.6, color="darkred", fatten=0) +
  theme_bw(base_size=12) +
  xlab(sprintf("Percent (method %d)", opt$col)) +
  ylab("") +
  theme(axis.text.y=element_text(size=6),  # 显示RBP名
        axis.ticks.y=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.border = element_blank())  # 移除图边框

# 8. 使用 ggtree 绘制树状图
ddata <- dendro_data(hc, type="rectangle")
p_den <- ggtree(hc, layout = "rectangular") +
  scale_x_continuous(expand = c(0, 0)) +
  theme_tree() +
  theme(legend.position="none")

# 9. 合并并保存
n_rbp <- length(order_rb)
h <- min(max(4, n_rbp*0.2 + 1), 50)
combined <- plot_grid(p_den, p_box, nrow=1, rel_widths=c(1,4))
ggsave(opt$output, combined, width=8, height=h, limitsize=TRUE)

cat(sprintf("已保存聚类水平箱线图至 %s (宽度=8\", 高度=%.1f\")\n", opt$output, h))
