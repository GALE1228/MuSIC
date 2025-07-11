#!/usr/bin/env Rscript

suppressMessages(library(optparse))
library(stringr)

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="输入TSV文件路径，无表头，包含id,f,kc,f,lenQ,lenR等列", metavar="file"),
  make_option(c("-o", "--output"), type="character", default="output_maxscore.tsv",
              help="输出文件路径 [default= %default]", metavar="file"),
  make_option(c("-w", "--weight"), type="double", default=0.8,
              help="id*kc权重，长度权重为1-w，范围[0,1] [default= %default]")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("必须提供输入文件路径 --input", call.=FALSE)
}

# 1. 读取原始数据（无表头）
df <- read.table(opt$input, header=FALSE, sep="\t", stringsAsFactors=FALSE, fill=TRUE)

# 2. 根据实际数据调整以下列号（示例根据你的示例调整）
# 假设：
# V1 = Query transcript ID
# V3 = lenQ (查询序列长度)
# V8 = lenR (参考序列长度)
# 最后两列分别是 id:f:0.xxxxxx 和 kc:f:0.xxxxxx

id_col <- ncol(df) - 1
kc_col <- ncol(df)

# 提取 id 和 kc 数值
df$id <- as.numeric(str_extract(df[[id_col]], "(?<=id:f:)\\d*\\.?\\d+"))
df$kc <- as.numeric(str_extract(df[[kc_col]], "(?<=kc:f:)\\d*\\.?\\d+"))

# 提取长度列
df$lenQ <- as.numeric(df$V3)
df$lenR <- as.numeric(df$V8)

# 3. 计算score函数
calculate_score_weighted <- function(id, kc, lenQ, lenR, w1) {
  maxLenDiff <- pmax(pmin(2*lenR, 2000), 1)  # 保证至少为1，避免除0
  norm_len_diff <- pmin(abs(lenQ - lenR), maxLenDiff) / maxLenDiff
  length_factor <- 1 - norm_len_diff
  score <- w1 * (id * kc) + (1 - w1) * length_factor
  return(score)
}


df$score <- mapply(calculate_score_weighted, df$id, df$kc, df$lenQ, df$lenR,
                   MoreArgs = list(w1 = opt$weight))

# 4. 取每对转录本组合(Query和Reference ID)的最大score记录
# 假设 Query ID 在V1，Reference ID 在V6
library(dplyr)

df_max <- df %>%
  group_by(RefID = V6) %>%
  slice_max(score, n = 1, with_ties = FALSE) %>%
  ungroup()


# 5. 保存结果，包含所有原始列 + score
write.table(df_max, file=opt$output, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

cat(sprintf("处理完成，最大score结果保存至 %s，记录数 %d\n", opt$output, nrow(df_max)))
