#!/usr/bin/env Rscript

# 加载必要的包
library(dplyr)
library(ggplot2)
library(optparse)

# 设置命令行参数
option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "输入文件路径，多个文件以逗号分隔"),
  make_option(c("-o", "--output"), type = "character", default = "forest_plot_output.pdf", help = "输出的 PDF 文件名 [默认: %default]"),
  make_option(c("-c", "--column"), type = "integer", default = 1, help = "选择用于 y 轴的列（1, 2, 或 3）[默认: %default]")
)

# 解析命令行参数
opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

# 读取输入文件路径
input_files <- strsplit(opts$input, ",")[[1]]

# 初始化一个空的数据框
all_data <- data.frame()

# 定义 species_name 的映射
species_map <- c(
  "species_1" = "PONAB",
  "species_2" = "MACFA",
  "species_3" = "MOUSE",
  "species_4" = "RAT",
  "species_5" = "CHICK",
  "species_6" = "XENLA",
  "species_7" = "DANRE",
  "species_8" = "DROME",
  "species_9" = "ARATH",
  "species_10" = "YEAST"
)

# 处理每个文件
for (file_path in input_files) {
  # 提取文件名和文件路径
  file_name <- basename(file_path)
  
  # 提取分组信息
  group_info <- sub(".*_(species_\\d+).*", "\\1", file_name)  # 获取文件名中的分组信息
  species_number <- sub("species_(\\d+)", "\\1", group_info)  # 提取 species_x 中的数字部分

  # 读取文件内容
  data <- read.table(file_path, header = FALSE)
  colnames(data) <- c("path", "col1", "col2", "col3")  # 给列命名
  
  # 将 col1, col2, col3 根据 species_number 进行处理，并转换为百分比
  data$col1 <- (data$col1 / (as.numeric(species_number) + 1)) * 100
  data$col2 <- (data$col2 / as.numeric(species_number)) * 100
  data$col3 <- (data$col3 / as.numeric(species_number)) * 100

  # 添加 RBP 名称和分组信息到数据框
  data$RBP <- sub("^.*/([A-Za-z0-9_]+)_.*$", "\\1", data$path)
  data$group <- group_info
  
  # 添加 species_name 列
  data$species_name <- species_map[group_info]
  
  data <- data %>%
    group_by(RBP, species_name) %>%
    mutate(mean_y = mean(.data[[paste("col", opts$column, sep = "")]]),
           sorted_col = sort(.data[[paste("col", opts$column, sep = "")]]),
           head_25_mean = mean(head(sorted_col, floor(length(sorted_col) * 0.25))),
           tail_25_mean = mean(tail(sorted_col, floor(length(sorted_col) * 0.25)))) %>%
    ungroup()  # 取消分组

  # 将数据绑定到 all_data
  all_data <- rbind(all_data, data)
}

# 确保选择的列存在
y_col <- paste("col", opts$column, sep = "")
if (!y_col %in% colnames(all_data)) {
  stop(paste("Column", opts$column, "does not exist in the data."))
}

# 确保x轴按照物种顺序
species_order <- c("PONAB", "MACFA", "MOUSE", "RAT", "CHICK", "XENLA", "DANRE", "DROME", "ARATH", "YEAST")
all_data$species_name <- factor(all_data$species_name, levels = species_order)

# 绘制森林图
p <- ggplot(all_data, aes(x = species_name, y = mean_y, color = mean_y)) +
  geom_point(position = position_dodge(width = 0.8), size = 3) +  # 绘制均值点
  geom_errorbar(aes(ymin = head_25_mean, ymax = tail_25_mean), width = 0.3, position = position_dodge(width = 0.8)) +  # 误差条为头部25%均值到尾部25%均值
  facet_grid(RBP ~ ., scales = "free_y") +
  scale_color_gradient(low = "#60a9cf", high = "#e8b6ae", limits = c(0, 100)) +  # 填充色根据 y 均值梯度变化
  scale_y_continuous(
    limits = c(0, 100),
    breaks = c(0, 50, 100),               # 显示的刻度标签
    minor_breaks = c(25, 75)              # 显示但不标数字的次刻度线
  ) +
  theme_minimal() +
  labs(
    title = "IGF2BP family of RBPs",  # 图表标题英文
    x = "Species (Accumulated)",  # x轴标题英文
    y = "Extent of Conservation of Binding Sites",  # y轴标题更新
    color = "Average Conservancy"  # 图例标题英文
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  # 保留四边边框
    panel.grid.major = element_blank(),  # 取消主要参考线
    panel.grid.minor = element_blank(),   # 取消次要参考线
    panel.spacing = unit(1.2, "lines")
  )

# 输出为 PDF 文件
ggsave(opts$output, plot = p, width = 8, height = 6)
