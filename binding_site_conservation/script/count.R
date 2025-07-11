#!/usr/bin/env Rscript
suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(stringr))

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="输入TSV文件路径，无表头，至少5列：ID,start,end,length,score", metavar="file")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input)) {
  stop("必须指定输入文件路径 -i")
}

# 读取数据
df <- fread(opt$input, header=FALSE, sep="\t", data.table=FALSE,
            col.names = c("ID","start","end","length","score"))

# 1. 第一列唯一元素个数
unique_ids <- unique(df$ID)
cat(sprintf("%s\t%d\t", opt$input, length(unique_ids)))

# 2. 分组：ENST开头与非ENST开头
enst_rows <- grepl("^ENST", df$ID)
df_enst <- df[enst_rows, , drop=FALSE]
df_non_enst <- df[!enst_rows, , drop=FALSE]

# 定义区间是否重叠函数，区间[start1,end1]和[start2,end2]
is_overlap <- function(start1, end1, start2, end2) {
  return(!(end1 < start2 || end2 < start1))
}

# 统计非ENST行是否与任意ENST行重叠
overlap_non_enst <- logical(nrow(df_non_enst))
for (i in seq_len(nrow(df_non_enst))) {
  non_st <- df_non_enst$start[i]
  non_en <- df_non_enst$end[i]
  # 判断是否与任意ENST行重叠
  overlap_any <- any(mapply(is_overlap,
                            non_st, non_en,
                            df_enst$start, df_enst$end))
  overlap_non_enst[i] <- overlap_any
}
# 计算重叠的非ENST行对应ID个数
ids_overlap_non_enst <- unique(df_non_enst$ID[overlap_non_enst])
cat(sprintf("%d\t", length(ids_overlap_non_enst)))

# 3. 对于每个ENST行，统计与之重叠的非ENST行ID唯一数，取最大值
max_overlap_count <- 0
for (j in seq_len(nrow(df_enst))) {
  enst_st <- df_enst$start[j]
  enst_en <- df_enst$end[j]

  overlap_flag <- mapply(is_overlap,
                         df_non_enst$start, df_non_enst$end,
                         MoreArgs = list(start2 = enst_st, end2 = enst_en))
  ids_overlap <- unique(df_non_enst$ID[overlap_flag])
  if (length(ids_overlap) > max_overlap_count) {
    max_overlap_count <- length(ids_overlap)
  }
}
cat(sprintf("%d\n", max_overlap_count))
