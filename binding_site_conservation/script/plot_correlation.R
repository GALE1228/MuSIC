#!/usr/bin/env Rscript
suppressMessages({
  library(data.table)
  library(ggplot2)
  library(optparse)
})

# 定义命令行参数
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

# 如果用户请求帮助，显示帮助信息并退出
if (opt$help) {
  print_help(parser)
  quit(status = 0)
}

# 读取数据
edgelist <- fread(opt$edgelist)
tsv_data <- fread(opt$tsv)

# 计算每个 RBP 在网络中的重要性（累加权重）
rbp_importance_from <- edgelist[, .(importance = sum(adjusted_weight)), by = from]
setnames(rbp_importance_from, "from", "RBP")
rbp_importance_to <- edgelist[, .(importance = sum(adjusted_weight)), by = to]
setnames(rbp_importance_to, "to", "RBP")
rbp_importance <- rbind(rbp_importance_from, rbp_importance_to)[, .(importance = sum(importance)), by = RBP]

# 合并数据
merged_data <- merge(rbp_importance, tsv_data, by = "RBP", all = TRUE)

# 计算相关性
cor_result <- cor.test(merged_data$importance, merged_data$mean_value)
r_value <- round(cor_result$estimate, 2)
p_value <- format.pval(cor_result$p.value, digits = 2)

# 计算点数 N
N <- nrow(merged_data)

# 准备相关性标注文本，包含 N，并标注 Pearson's
cor_label <- paste0("Pearson's r = ", r_value, "\np = ", p_value, "\nN = ", N)

# 绘制相关分析图，添加回归线
p <- ggplot(merged_data, aes(x = importance, y = mean_value)) +
  geom_point() +
  # 添加线性回归线，使用 lm 方法，不显示置信区间
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  # 明确设置 x、y 轴英文名称，缩小字体
  labs(x = "Weighted Degree Centrality of RBP", y = "Cross-species Conservation of RBP Binding Sites") +
  theme_minimal(base_size = 8) +
  # 设置横坐标刻度
  scale_x_continuous(breaks = c(0, 10, 20)) +
  # 设置纵坐标刻度
  scale_y_continuous(breaks = c(0, 50, 100)) +
  # 添加相关性统计指标标注到右下角
  annotate("text", x = max(merged_data$importance, na.rm = TRUE), 
           y = min(merged_data$mean_value, na.rm = TRUE), 
           label = cor_label, hjust = 1, vjust = 0) +
  # 去掉纵横辅助线，保留 xy 轴线
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )

# 保存图形
ggsave(opt$output, p, width = 4, height = 3)

cat(sprintf("Correlation plot has been saved to %s\n", opt$output))