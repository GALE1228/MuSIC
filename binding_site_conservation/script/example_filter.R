#!/usr/bin/env Rscript
suppressMessages({
  library(optparse)
  library(data.table)
  library(stringr)
  library(Biostrings)
})

option_list <- list(
  make_option(c("-u", "--updated"), type="character", default=NULL,
              help="更新后（基于多序列比对坐标）的BED/TSV文件，至少5列：ID, start_upd, end_upd, extra1, extra2", metavar="file"),
  make_option(c("-r", "--raw"), type="character", default=NULL,
              help="原始BED文件（convert.R运行前），格式“ID|start_raw|end_raw|extra1|extra2”，无表头", metavar="file"),
  make_option(c("-f", "--fasta"), type="character", default=NULL,
              help="多序列比对FASTA文件路径，用于坐标还原", metavar="file"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="输出文件路径，仅包含与最佳ENST条目重叠的原始BED条目（按ID、坐标排序）", metavar="file")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$updated) || is.null(opt$raw) || is.null(opt$fasta)) {
  stop("必须同时指定 -u 更新后文件, -r 原始BED文件, -f FASTA")
}

# 1. 读取更新后BED/TSV，假定无表头，\t分隔
df_upd <- fread(opt$updated, header=FALSE, sep="\t", data.table=FALSE,
                col.names = c("ID", "start_upd", "end_upd", "extra1", "extra2"))
if (!all(c("ID","start_upd","end_upd") %in% colnames(df_upd))) {
  stop("更新后文件需包含至少三列：ID, start_upd, end_upd")
}

# 2. 读取原始BED，假定“|”分隔，无表头
raw_df <- fread(opt$raw, header=FALSE, sep="|", data.table=FALSE,
                col.names = c("ID", "start_raw", "end_raw", "extra1", "extra2"))
if (nrow(df_upd) != nrow(raw_df)) {
  stop("更新后文件与原始BED行数需一致，确保行一一对应")
}

# 3. 载入MSA，用于坐标还原
msa <- readDNAMultipleAlignment(opt$fasta, format="fasta")
seq_names <- rownames(msa)
if (length(seq_names) == 0) {
  stop("无法从FASTA中读取任何序列")
}

# 辅助：MSA反向映射到原始坐标
get_orig_pos <- function(msa_seq, msa_pos) {
  chars <- as.character(msa_seq)
  sub   <- substring(chars, 1, msa_pos)
  return(nchar(gsub("-", "", sub)))
}
find_msa_id <- function(id) {
  idx <- grep(paste0("^", id), seq_names)
  if (length(idx) == 0) return(NA_character_)
  return(seq_names[idx[1]])
}

# 4. 分离ENST与非ENST
is_enst <- grepl("^ENST", df_upd$ID)
df_enst <- df_upd[is_enst, , drop=FALSE]
df_non  <- df_upd[!is_enst, , drop=FALSE]

# 5. 区间重叠判断
is_overlap <- function(s1,e1,s2,e2) !(e1 < s2 || e2 < s1)

# 6. 计算每个ENST与不同非ENST的重叠数
enst_indices  <- which(is_enst)
n_enst        <- length(enst_indices)
overlap_count <- integer(n_enst)
for (i in seq_len(n_enst)) {
  idx_i <- enst_indices[i]
  st_i  <- df_upd$start_upd[idx_i]
  en_i  <- df_upd$end_upd[idx_i]
  mask  <- mapply(is_overlap,
                  s1 = df_non$start_upd, e1 = df_non$end_upd,
                  s2 = st_i,          e2 = en_i)
  overlap_ids <- unique(df_non$ID[mask])
  overlap_count[i] <- length(overlap_ids)
}

# 7. 找到最大重叠数，选取所有达到该值的ENST行
max_overlap <- if (n_enst > 0) max(overlap_count) else 0
if (max_overlap == 0) {
  stop("没有ENST行与任何非ENST行重叠")
}
best_enst_rows <- enst_indices[overlap_count == max_overlap]

# 8. 对每个最佳ENST行，收集该行及其所有重叠的非ENST行
all_indices <- integer(0)
for (idx_i in best_enst_rows) {
  st_i <- df_upd$start_upd[idx_i]
  en_i <- df_upd$end_upd[idx_i]
  mask <- mapply(is_overlap,
                 s1 = df_non$start_upd, e1 = df_non$end_upd,
                 s2 = st_i,          e2 = en_i)
  non_indices <- which(!is_enst)[mask]
  all_indices <- c(all_indices, idx_i, non_indices)
}
all_indices <- unique(all_indices)

# 9. 逐行还原原始坐标并收集输出
out_list <- list()
for (idx in all_indices) {
  row_upd <- df_upd[idx, ]
  id      <- row_upd$ID
  msa_id  <- find_msa_id(id)
  if (is.na(msa_id)) {
    warning(sprintf("未找到MSA序列 %s，跳过", id))
    next
  }
  msa_idx <- which(seq_names == msa_id)
  if (length(msa_idx) == 0) {
    warning(sprintf("MSA中未找到序列 %s，跳过", msa_id))
    next
  }
  msa_seq <- msa@unmasked[msa_idx][[1]]
  orig_st <- get_orig_pos(msa_seq, row_upd$start_upd)
  orig_en <- get_orig_pos(msa_seq, row_upd$end_upd)
  if (is.na(orig_st) || is.na(orig_en)) {
    warning(sprintf("ID %s 的坐标映射失败，跳过", id))
    next
  }
  raw_row <- raw_df[idx, ]
  out_list[[length(out_list)+1]] <- data.frame(
    ID     = id,
    start  = orig_st,
    end    = orig_en,
    extra1 = raw_row$extra1,
    extra2 = raw_row$extra2,
    stringsAsFactors = FALSE
  )
}

if (length(out_list) == 0) {
  stop("未生成任何输出条目")
}
out_df <- rbindlist(out_list)

# 10. 按ID字母和start升序排序
out_df <- out_df[order(out_df$ID, out_df$start), ]

# 11. 写出结果
if (!is.null(opt$output)) {
  fwrite(out_df, file = opt$output, sep="|", quote=FALSE, col.names=FALSE, row.names=FALSE)
  cat(sprintf("已输出 %d 条排序后的原始BED条目到 %s\n", nrow(out_df), opt$output))
} else {
  fwrite(out_df, file="", sep="|", quote=FALSE, col.names=FALSE, row.names=FALSE)
}
