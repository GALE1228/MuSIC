#!/usr/bin/env Rscript
suppressMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(stringr)
  library(reshape2)
  library(scales)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="输入统计结果文件，第一列为“RBP_Transcript.bed”，其余列是各物种0/1值", metavar="file"),
  make_option(c("-o", "--output"), type="character", default="heatmap.pdf",
              help="输出热图PDF路径 [default %default]", metavar="file")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input)) {
  stop("必须指定输入文件路径 -i")
}

# 1. 读取输入，保留列名
dt <- fread(opt$input, header=TRUE, sep="\t", data.table=FALSE)

# 2. 提取 RBP 名称（去掉下划线及之后部分）
colnames(dt)[1] <- "File"
dt$RBP <- str_split_fixed(basename(dt$File), "_", 2)[,1]

# 3. 取出物种列（除 File 和 RBP 外），转换为整数
species_cols <- setdiff(colnames(dt), c("File","RBP"))
dt[species_cols] <- lapply(dt[species_cols], as.integer)

# 4. 按 RBP 汇总：各物种列求和，生成 RBP × Species 的矩阵
agg_dt <- aggregate(dt[, species_cols, drop=FALSE],
                    by = list(RBP = dt$RBP),
                    FUN = sum)
mat <- as.matrix(agg_dt[, species_cols])
rownames(mat) <- agg_dt$RBP

# 5. 固定物种顺序（HUMAN 放最前，其余按既定次序）
desired_order <- c("HUMAN", "PONAB", "MACFA", "MOUSE", "RAT",
                   "CHICK", "XENLA", "DANRE", "DROME", "ARATH", "YEAST")
col_order <- intersect(desired_order, colnames(mat))
if (length(col_order) == 0) {
  stop("输入文件中的物种列与预定义顺序不匹配，请检查列名")
}

# 6. 对 RBP 按所有物种之和从低到高排序（先低后高）
row_sums <- rowSums(mat[, col_order, drop=FALSE])
row_order <- names(sort(row_sums, decreasing = FALSE))

# 7. 转换为长格式：列 “RBP”, “Species”, “Count”，并按顺序设置因子
long_df <- melt(mat, varnames = c("RBP", "Species"), value.name = "Count")
long_df$RBP <- factor(long_df$RBP, levels = row_order)
long_df$Species <- factor(long_df$Species, levels = col_order)

# 8. 动态调整画幅尺寸并限制在50英寸以内
n_rbp <- length(row_order)
n_species <- length(col_order)
# RBP 行多，使用较窄高度：每行约0.15英寸，最小4英寸，上限50
height_in <- min(max(4, n_rbp * 0.15 + 1), 50)
# 物种列较少，每列约0.6英寸以便更好展示，最小4英寸，上限50
width_in  <- min(max(4, n_species * 0.6 + 1), 50)

# 9. 绘制热图，将物种名放到顶部，减小字体
p <- ggplot(long_df, aes(x = Species, y = RBP, fill = Count)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue",
                      name = "Overlap Count", labels = comma) +
  scale_x_discrete(position = "top") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 6),
    axis.text.x.bottom = element_blank(),
    axis.text.y = element_text(size = 5),
    axis.title   = element_blank(),
    panel.grid   = element_blank()
  )

# 10. 保存到 PDF
ggsave(opt$output, plot = p, width = width_in, height = height_in, limitsize = TRUE)

cat(sprintf("已将按总和从低到高排序并调整物种名和尺寸的热图保存至 %s (width=%.1f\", height=%.1f\")\n",
            opt$output, width_in, height_in))
