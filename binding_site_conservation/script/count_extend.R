#!/usr/bin/env Rscript
suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(stringr))

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="输入 TSV 文件路径，无表头，至少 5 列：ID,start,end,length,score", metavar="file"),
  make_option(c("-e", "--extend"), type="integer", default=0,
              help="范围拓宽值，区域左右各向外拓宽 N 后再判断重叠 [default: %default]", metavar="N")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input)) {
  stop("必须指定输入文件路径 -i")
}

# 读取数据
df <- fread(opt$input, header=FALSE, sep="\t", data.table=FALSE,
            col.names = c("ID","start","end","length","score"))

# 先拓宽区间
df$start_ext <- pmax(0, df$start - opt$extend)
df$end_ext   <- df$end   + opt$extend

# 1. 第一列唯一元素个数
unique_ids <- unique(df$ID)
cat(sprintf("%s\t%d\t", opt$input, length(unique_ids)))

# 2. 分组：ENST 开头与非 ENST 开头
enst_rows    <- grepl("^ENST", df$ID)
df_enst      <- df[enst_rows, , drop=FALSE]
df_non_enst  <- df[!enst_rows, , drop=FALSE]

# 定义区间是否重叠函数
is_overlap <- function(start1, end1, start2, end2) {
  !(end1 < start2 || end2 < start1)
}

# 统计非ENST行是否与任意ENST行重叠
overlap_non_enst <- vapply(seq_len(nrow(df_non_enst)), function(i) {
  any(mapply(is_overlap,
             df_non_enst$start_ext[i], df_non_enst$end_ext[i],
             df_enst$start_ext,    df_enst$end_ext))
}, logical(1))

# 计算重叠的非ENST行对应 ID 个数
ids_overlap_non_enst <- unique(df_non_enst$ID[overlap_non_enst])
cat(sprintf("%d\t", length(ids_overlap_non_enst)))

# 3. 对于每个 ENST 行，统计与之重叠的非ENST行 ID 唯一数，取最大值
max_overlap_count <- 0
for (j in seq_len(nrow(df_enst))) {
  st_j <- df_enst$start_ext[j]
  en_j <- df_enst$end_ext[j]
  overlap_flag <- mapply(is_overlap,
                         df_non_enst$start_ext, df_non_enst$end_ext,
                         MoreArgs = list(start2 = st_j, end2 = en_j))
  ids_ov <- unique(df_non_enst$ID[overlap_flag])
  if (length(ids_ov) > max_overlap_count) {
    max_overlap_count <- length(ids_ov)
  }
}
cat(sprintf("%d\n", max_overlap_count))
