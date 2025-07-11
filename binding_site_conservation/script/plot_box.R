#!/usr/bin/env Rscript
suppressMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(stringr)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="输入文件，格式：<path>/<RBP>_<transcript>.bed <val1> <val2> <val3>", metavar="file"),
  make_option(c("-c", "--col"), type="integer", default=1,
              help="选择绘图的列 (1-3) [default %default]", metavar="int"),
  make_option(c("-o", "--output"), type="character", default="boxplot.pdf",
              help="输出PDF文件路径 [default %default]", metavar="file")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input)) {
  stop("必须指定输入文件路径 -i")
}
if (!(opt$col %in% 1:3)) {
  stop("参数 -c 必须是 1、2 或 3")
}

# 读取数据
df <- fread(opt$input, header=FALSE, sep="\t", data.table=FALSE,
            col.names = c("path", "V1", "V2", "V3"))

# 提取 RBP 与 transcript
df$basename <- basename(df$path)
df$basename <- str_replace(df$basename, "\\.bed$", "")
parts <- str_split_fixed(df$basename, "_", 2)
df$RBP <- parts[,1]
df$transcript <- parts[,2]

# 归一化
df$norm1 <- df$V1 / 11 * 100
df$norm2 <- df$V2 / 10 * 100
df$norm3 <- df$V3 / 10 * 100

col_name <- paste0("norm", opt$col)
plot_df <- data.frame(RBP = df$RBP, value = df[[col_name]], stringsAsFactors = FALSE)

# 按每个 RBP 的均值排序
mean_by_rbp <- aggregate(value ~ RBP, data = plot_df, FUN = mean)
ordered_rbps <- mean_by_rbp[order(mean_by_rbp$value), "RBP"]

# 设置因子水平
plot_df$RBP <- factor(plot_df$RBP, levels = ordered_rbps)

# 根据RBP数量调整宽度：每个类别约0.2英寸，最小6英寸
n_rbp <- length(ordered_rbps)
plot_width <- max(6, n_rbp * 0.2)
plot_height <- 4

# 绘制 boxplot，用均值线替代中位数线
p <- ggplot(plot_df, aes(x = RBP, y = value)) +
  geom_boxplot(fill = "steelblue", outlier.size = 0.8, median.size = 0) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.6, color = "darkred", fatten = 0) +
  theme_bw() +
  xlab("RBP") +
  ylab(sprintf("Percent of species by method %d (%%)", opt$col)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 保存到 PDF（限制尺寸不超过50英寸）
ggsave(opt$output, plot = p, width = plot_width, height = plot_height, limitsize = TRUE)

cat(sprintf("已将替换为平均值线的 boxplot 保存至 %s\n", opt$output))
