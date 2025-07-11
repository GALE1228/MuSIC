#!/usr/bin/env Rscript
suppressMessages({
  library(optparse)
  library(Biostrings)
  library(data.table)
  library(stringr)
})

option_list <- list(
  make_option(c("-b", "--bed"), type="character", default=NULL,
              help="输入bed文件，格式ID|start|end|...，无表头", metavar="file"),
  make_option(c("-f", "--fasta"), type="character", default=NULL,
              help="多序列比对FASTA文件路径", metavar="file"),
  make_option(c("-o", "--outbed"), type="character", default="updated.bed",
              help="输出更新后bed文件路径", metavar="file")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$bed) || is.null(opt$fasta)) {
  stop("必须同时指定 -b 输入bed文件 和 -f 多序列比对FASTA文件")
}

# 读取bed文件
bed_raw <- fread(opt$bed, sep="|", header=FALSE, col.names=c("ID", "start", "end", "extra1", "extra2"))
cat(sprintf("读取 %d 条BED位点\n", nrow(bed_raw)))

# 读取多序列比对FASTA为DNAMultipleAlignment对象
msa <- readDNAMultipleAlignment(opt$fasta, format="fasta")
seq_names <- rownames(msa)
cat(sprintf("读取 %d 条多序列比对序列\n", length(seq_names)))

# gap-aware坐标映射函数
get_msa_pos <- function(seq, orig_pos) {
  count <- 0
  for (i in seq_along(seq)) {
    if (substr(as.character(seq), i, i) != "-") count <- count + 1
    if (count == orig_pos) return(i)
  }
  return(NA)
}

# 简单前缀匹配找到比对序列名
find_msa_seq <- function(id) {
  matched <- grep(paste0("^", id), seq_names, value=TRUE)
  if (length(matched) == 0) return(NA)
  return(matched[1])
}

bed_updated <- copy(bed_raw)
for (i in seq_len(nrow(bed_raw))) {
  id <- bed_raw$ID[i]
  msa_id <- find_msa_seq(id)
  if (is.na(msa_id)) {
    warning(sprintf("未找到比对序列 %s，保留原坐标", id))
    next
  }
  idx <- which(rownames(msa) == msa_id)
  if (length(idx) == 0) {
    warning(sprintf("未找到比对序列索引 %s，保留原坐标", msa_id))
    next
  }
  seq_set <- msa@unmasked[idx]
  seq <- seq_set[[1]]
  start_orig <- bed_raw$start[i]
  end_orig <- bed_raw$end[i]
  start_msa <- get_msa_pos(seq, start_orig)
  end_msa <- get_msa_pos(seq, end_orig)
  if (is.na(start_msa) || is.na(end_msa)) {
    warning(sprintf("ID %s 的位点转换失败，保留原坐标", id))
    next
  }
  bed_updated$start[i] <- start_msa
  bed_updated$end[i] <- end_msa
}

fwrite(bed_updated, file=opt$outbed, sep="\t", col.names=FALSE)
cat(sprintf("更新后的bed文件已保存至：%s\n", opt$outbed))
