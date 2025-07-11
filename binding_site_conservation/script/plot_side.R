#!/usr/bin/env Rscript
suppressMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(stringr)
  library(dplyr)
  library(reshape2)
  library(scales)
  library(patchwork)
})

# 禁用默认 -h（帮助）短选项，避免冲突
parser <- OptionParser(add_help_option=FALSE,
                       option_list = list(
                         make_option(c("--help"), action="store_true", default=FALSE,
                                     help="显示帮助信息然后退出"),
                         make_option(c("-H", "--heatmap"), type="character", default=NULL,
                                     help="热图输入：第一列为“RBP_Transcript.bed”，其余列是各物种0/1值", metavar="file"),
                         make_option(c("-B", "--boxplot"), type="character", default=NULL,
                                     help="箱线图输入：格式 <path>/<RBP>_<transcript>.bed <val1> <val2> <val3>", metavar="file"),
                         make_option(c("-c", "--col"), type="integer", default=1,
                                     help="箱线图选择的列 (1-3) [default %default]", metavar="int"),
                         make_option(c("-o", "--output"), type="character", default="combined_plot.pdf",
                                     help="输出合并图PDF路径 [default %default]", metavar="file")
                       ))
opt <- parse_args(parser)

if (opt$help) {
  print_help(parser); q(status=0)
}
if (is.null(opt$heatmap) || is.null(opt$boxplot)) {
  stop("必须同时指定 -H/--heatmap 和 -B/--boxplot 输入文件")
}
if (!(opt$col %in% 1:3)) {
  stop("参数 -c 必须是 1、2 或 3")
}

#### 1. 处理热图数据 ####
dt_h <- fread(opt$heatmap, header=TRUE, sep="\t", data.table=FALSE)
colnames(dt_h)[1] <- "File"
dt_h$RBP <- str_split_fixed(basename(dt_h$File), "_", 2)[,1]
species_cols <- setdiff(colnames(dt_h), c("File","RBP"))
dt_h[species_cols] <- lapply(dt_h[species_cols], as.integer)

# 统计每个 RBP 的转录本总数，作为分母
transcript_counts <- as.integer(table(dt_h$RBP))
names(transcript_counts) <- names(table(dt_h$RBP))

# 按 RBP 汇总各物种 0/1 值求和
agg_h <- aggregate(dt_h[, species_cols, drop=FALSE],
                   by = list(RBP = dt_h$RBP), FUN = sum)
mat_h <- as.matrix(agg_h[, species_cols])
rownames(mat_h) <- agg_h$RBP

desired_order <- c("HUMAN","PONAB","MACFA","MOUSE","RAT",
                   "CHICK","XENLA","DANRE","DROME","ARATH","YEAST")
col_order <- intersect(desired_order, colnames(mat_h))
if (length(col_order)==0) {
  stop("热图输入文件中的物种列与预设顺序不匹配")
}

# 排序 RBP：按所有物种之和从低到高
row_sums_h <- rowSums(mat_h[, col_order, drop=FALSE])
row_order <- names(sort(row_sums_h, decreasing = FALSE))

# 转换为长格式，并计算百分比：Count / transcript_count * 100
long_h <- melt(mat_h, varnames=c("RBP","Species"), value.name="Count")
long_h$RBP <- factor(long_h$RBP, levels=row_order)
long_h$Species <- factor(long_h$Species, levels=col_order)
long_h$Percent <- long_h$Count / transcript_counts[as.character(long_h$RBP)] * 100

#### 2. 处理箱线图数据 ####
df_b <- fread(opt$boxplot, header=FALSE, sep="\t", data.table=FALSE,
              col.names=c("path","V1","V2","V3"))
df_b$basename <- basename(df_b$path)
df_b$basename <- str_replace(df_b$basename, "\\.bed$", "")
parts_b <- str_split_fixed(df_b$basename, "_", 2)
df_b$RBP <- parts_b[,1]

df_b$norm1 <- df_b$V1/11*100
df_b$norm2 <- df_b$V2/10*100
df_b$norm3 <- df_b$V3/10*100
col_name <- paste0("norm", opt$col)
plot_b <- data.frame(RBP = df_b$RBP, value = df_b[[col_name]], stringsAsFactors=FALSE)

extra_rb <- setdiff(unique(plot_b$RBP), row_order)
combined_order <- c(row_order, extra_rb)
plot_b$RBP <- factor(plot_b$RBP, levels=combined_order)

#### 3. 构建 ggplot 对象 ####
# 计算每个 RBP 分组的均值
group_means <- plot_b %>%
  group_by(RBP) %>%
  summarise(mean_value = mean(value))

# 将均值合并到原数据框
plot_b <- left_join(plot_b, group_means, by = "RBP")

# 箱线图：RBP 在 x 轴，根据均值填充颜色
p_boxplot <- ggplot(plot_b, aes(x = RBP, y = value, fill = mean_value)) +
  geom_boxplot(outlier.size = 0.5, median.size = 0) +
  # 设置填充颜色为梯度色
  scale_fill_gradient(low = "#60a9cf", high = "#e8b6ae", 
                      name = "Mean Value", 
                      guide = "colorbar") +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 270, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  ylab(sprintf("Percent of species by method %d (%%)", opt$col))

# 热图：RBP 在 x 轴，填色用 Percent
p_heatmap <- ggplot(long_h, aes(x = RBP, y = Species, fill = Percent)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "#60a9cf", high = "#e8b6ae",
                      name = "Percentage of conservation (%)", labels = function(x) sprintf("%.1f%%", x)) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 270, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )

#### 4. 计算尺寸并上下排列（保持 RBP 对齐） ####
n_rbp <- length(combined_order)
n_species <- length(col_order)

# 箱线图高度： n_rbp*0.03 + 1，最小2，上限8
box_h   <- min(max(2, n_rbp * 0.03 + 1), 5)
# 热图高度： n_species*0.3 + 1，最小4，上限20
heat_h  <- min(max(4, n_species * 0.3 + 1), 20)
total_h <- box_h + heat_h

# 宽度按 RBP 数量，每 RBP*0.2 + 2，最小6，上限50
total_w <- min(max(6, n_rbp * 0.2 + 2), 50)

combined_plot <- p_boxplot / p_heatmap +
  plot_layout(heights = c(box_h, heat_h)) &
  theme(plot.margin = margin(5,5,5,5))

ggsave(opt$output, combined_plot, width = total_w, height = total_h, limitsize=TRUE)

# 生成 TSV 文件路径
tsv_path <- sub("\\.pdf$", ".tsv", opt$output)

# 输出 RBP 的平均保守性值到 TSV 文件
write.table(group_means, file = tsv_path, sep = "\t", na = "nan", 
            quote = FALSE, row.names = FALSE)

cat(sprintf("已将上下排列且 RBP 作为横坐标的组合图保存至 %s (宽度=%.1f\", 高度=%.1f\")\n",
            opt$output, total_w, total_h))
cat(sprintf("已将 RBP 的平均保守性值保存至 %s\n", tsv_path))