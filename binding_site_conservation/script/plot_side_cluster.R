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
  library(pheatmap)  # 添加 pheatmap 用于热图层次聚类
})

# 禁用默认 -h（帮助）短选项，避免冲突
parser <- OptionParser(add_help_option=FALSE,
                       option_list = list(
                         make_option(c("--help"), action="store_true", default=FALSE,
                                     help="显示帮助信息然后退出"),
                         make_option(c("-H", "--heatmap"), type="character", default=NULL,
                                     help="热图输入：第一列为“RBP_Transcript.bed”，其余列是各物种0/1值", metavar="file"),
                         make_option(c("-o", "--output"), type="character", default="heatmap_cluster_plot.pdf",
                                     help="输出热图PDF路径 [default %default]", metavar="file")
                       ))
opt <- parse_args(parser)

if (opt$help) {
  print_help(parser); q(status=0)
}
if (is.null(opt$heatmap)) {
  stop("必须指定 -H/--heatmap 输入文件")
}

#### 1. 处理热图数据 ####
dt_h <- fread(opt$heatmap, header=TRUE, sep="\t", data.table=FALSE)
colnames(dt_h)[1] <- "File"
dt_h$RBP <- str_split_fixed(basename(dt_h$File), "_", 2)[,1]
species_cols <- setdiff(colnames(dt_h), c("File","RBP"))
dt_h[species_cols] <- lapply(dt_h[species_cols], as.integer)

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
long_h$Percent <- long_h$Count / row_sums_h[as.character(long_h$RBP)] * 100

#### 2. 层次聚类 ####
# 使用 pheatmap 包来进行热图绘制，并在此过程中执行层次聚类
pheatmap(mat_h, 
         cluster_rows = TRUE,   # RBP行的层次聚类
         color = colorRampPalette(c("#60a9cf", "#e8b6ae"))(100),  # 设置颜色
         main = "RBP Heatmap with Hierarchical Clustering",
         show_rownames = TRUE,  # 显示行名（RBP）
         show_colnames = TRUE,  # 显示列名（物种）
         fontsize = 8,
         cellwidth = 20,        # 单元格宽度
         cellheight = 10,       # 单元格高度
         border_color = "black", # 边框颜色
         filename = opt$output,  # 输出文件路径
         cutree_rows = 3,        # 按行切割为3类
)

cat(sprintf("已将热图保存至 %s\n", opt$output))
