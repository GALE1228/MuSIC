#!/usr/bin/env Rscript
suppressMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(stringr)
  library(tidyverse)    # for file handling (grep, list.files)
})

# 解析参数
option_list <- list(
  make_option(c("-i", "--input"),   type="character", help="输入多个物种文件，格式：<path>/species_1.txt,<path>/species_2.txt，支持正则表达式", metavar="file"),
  make_option(c("-c", "--col"),     type="integer",   default=1, help="选择列 (1-3) [default %default]", metavar="int"),
  make_option(c("-p", "--pattern"), type="character", default=".*", help="RBP 名正则，只画匹配的 [default all]", metavar="regex"),
  make_option(c("-o", "--output"),  type="character", default="violin_plot.pdf", help="输出 PDF 路径", metavar="file")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$input)) stop("必须指定 -i")
if (!(opt$col %in% 1:3)) stop("参数 -c 必须是 1、2 或 3")

# 1. 使用正则表达式读取多个文件
input_files <- list.files(path = dirname(opt$input), pattern = basename(opt$input), full.names = TRUE)
if (length(input_files) == 0) {
  stop("没有符合正则表达式的文件")
}

# 2. 读取每个文件并合并
species_list <- lapply(input_files, fread, header=FALSE, sep="\t", data.table=FALSE)

# 3. 设置列名
col_names <- c("bed", "col1", "col2", "col3")
df_list <- lapply(species_list, function(file_data) {
  colnames(file_data) <- col_names
  return(file_data)
})

# 合并所有物种文件
df <- do.call(rbind, df_list)

# 4. 归一化处理
normalize_col <- function(df, col_index, num_species) {
  return(df[[paste0("col", col_index)]] / num_species * 100)
}

df$value <- normalize_col(df, opt$col, length(input_files))

# 5. 正则过滤
plot_df <- df[grepl(opt$pattern, df$bed), ]
if (nrow(plot_df) == 0) stop("没有 RBP 符合正则")

# 6. 按照RBP分类并绘制每个RBP的每个文件的violin图
# 通过列出每个文件和RBP的所有值，来构建数据
plot_df$RBP <- factor(plot_df$bed)

# 7. 绘制小提琴图
p <- ggplot(plot_df, aes(x = RBP, y = value, fill = RBP)) +
  geom_violin(trim = TRUE, show.legend = FALSE) +  # 去除图例
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "darkred") +
  theme_bw() +
  xlab(sprintf("Percent (method %d)", opt$col)) +
  ylab("RBP") +
  theme(axis.text.x = element_text(size = 8),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.border = element_blank(),
        legend.position = "none") +  # 隐藏图例
  facet_wrap(~bed, scales = "free_y")  # 每个RBP一行，每行显示多个文件的小提琴图

# 输出设置
n_rbp <- length(unique(plot_df$RBP))
height <- min(max(4, n_rbp * 0.3), 50)

# 保存结果
pdf(opt$output, width = 10, height = height)
print(p)
dev.off()

cat(sprintf("已保存小提琴图至 %s (宽度=10\", 高度=%.1f\")\n", opt$output, height))
