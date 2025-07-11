#!/usr/bin/env Rscript

# 加载库
suppressMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(stringr)
  library(ggtree)
  library(ggtreeExtra)  # 注意这里是 ggtreeExtra（大小写敏感）
  library(cluster)      # for hclust
  library(ggdendro)     # 用于 dendro_data
})

# 解析参数
option_list <- list(
  make_option(c("-i", "--input"),   type="character", help="输入文件，格式：<path>/<RBP>_<transcript>.bed <v1> <v2> <v3>", metavar="file"),
  make_option(c("-c", "--col"),     type="integer",   default=1, help="选择列 (1-3) [default %default]", metavar="int"),
  make_option(c("-p", "--pattern"), type="character", default=".*", help="RBP 名正则，只画匹配的 [default all]", metavar="regex"),
  make_option(c("-o", "--output"),  type="character", default="boxplot_cluster_ggtree.pdf", help="输出 PDF 路径", metavar="file")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$input)) stop("必须指定 -i")
if (!(opt$col %in% 1:3)) stop("参数 -c 必须是 1、2 或 3")

# 1. 读取数据并提取 RBP
df <- fread(opt$input, header=FALSE, sep="\t",
            col.names=c("path","V1","V2","V3"), data.table=FALSE)

# 提取文件名中的 RBP 名称
df$basename <- basename(df$path)
df$basename <- str_replace(df$basename, "\\.bed$", "")
parts <- str_split_fixed(df$basename, "_", 2)
df$RBP <- parts[,1]

# 2. 归一化
norm_funcs <- list(
  function(x) x / 11 * 100,
  function(x) x / 10 * 100,
  function(x) x / 10 * 100
)
df$value <- norm_funcs[[opt$col]](df[[paste0("V", opt$col)]])
plot_df <- df[, c("RBP", "value")]

# 3. 正则过滤
plot_df <- plot_df[grepl(opt$pattern, plot_df$RBP), ]
if (nrow(plot_df) == 0) stop("没有 RBP 符合正则")

# 4. 汇总均值用于聚类
mean_df <- aggregate(value ~ RBP, data = plot_df, FUN = mean)

# 如果只有一个 RBP，绘制简单箱线图
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
  cat("仅 1 个 RBP，无聚类树，已输出简单箱线图\n")
  quit(status = 0)
}

# 5. 构建聚类树
mat <- as.matrix(mean_df$value)
rownames(mat) <- mean_df$RBP
dist_mat <- dist(mat)
hclust_tree <- hclust(dist_mat, method = "average")

# 转换为 dendrogram 类型供 ggtree 使用
dend <- as.dendrogram(hclust_tree)

# 6. 构造数据框用于 ggtree
tree_data <- ggdendro::dendro_data(dend, type = "rectangle")
label_data <- tree_data$label

# 7. 箱线图作为“水果”添加到叶节点
p <- ggtree(tree_data, layout = "rectangular") +
  geom_tiplab(aes(label = label), offset = 0.05, align = TRUE, linesize = 0.3) +
  geom_fruit(
    data = plot_df,
    geom = geom_boxplot,
    mapping = aes(x = value, y = RBP),  # ✅ 添加 y = RBP 映射
    width = 0.6,
    fill = "steelblue",
    outlier.size = 0.5
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_tree() +
  theme(legend.position = "none")

# 添加平均值标记
p <- p + geom_fruit(
  data = mean_df,
  geom = geom_crossbar,
  mapping = aes(xmin = value, xmax = value, y = RBP),  # ✅ y 已定义
  color = "darkred"
)

# 输出设置
n_rbp <- nrow(mean_df)
height <- min(max(4, n_rbp * 0.3), 50)

# 保存结果
pdf(opt$output, width = 10, height = height)
print(p)
dev.off()

cat(sprintf("已保存聚类水平箱线图至 %s (宽度=10\", 高度=%.1f\")\n", opt$output, height))