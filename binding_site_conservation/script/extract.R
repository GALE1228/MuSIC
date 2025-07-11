#!/usr/bin/env Rscript
suppressMessages(library(optparse))
suppressMessages(library(Biostrings))
suppressMessages(library(data.table))
suppressMessages(library(stringr))

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="过滤后TSV文件路径，无表头，包含Query ID和Ref ID列", metavar="file"),
  make_option(c("-f", "--fasta"), type="character", default=NULL,
              help="包含所有query和ref序列的FASTA文件路径", metavar="file"),
  make_option(c("-o", "--outdir"), type="character", default="output_fasta_by_ref",
              help="输出目录，默认 'output_fasta_by_ref'", metavar="dir")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$input) || is.null(opt$fasta)) {
  cat("必须指定输入TSV和FASTA文件路径\n")
  quit(status=1)
}

if (!dir.exists(opt$outdir)) {
  dir.create(opt$outdir, recursive = TRUE)
}

df <- fread(opt$input, header=FALSE, sep="\t", data.table=FALSE)

query_col <- 1
ref_col <- 6

ref_to_queries <- split(df[[query_col]], df[[ref_col]])

seqs <- readDNAStringSet(opt$fasta)
seq_names <- names(seqs)

# 建立一个映射：简化ID(以 | 或 空格分割的第1段) -> 完整seq名称
simplify_id <- function(full_id) {
  str_split(full_id, "\\||\\s")[[1]][1]
}

id_map <- setNames(seq_names, sapply(seq_names, simplify_id))

# 根据简化ID找完整seq名
find_seq_name <- function(id) {
  matched_names <- id_map[names(id_map) == id]
  if (length(matched_names) == 0) {
    # 尝试前缀匹配：匹配以id开头的seq_name
    prefix_matches <- grep(paste0("^", id), seq_names, value=TRUE)
    if(length(prefix_matches) > 0) {
      return(prefix_matches[1])
    }
    return(NA)
  }
  return(matched_names[1])
}

for (ref_id in names(ref_to_queries)) {
  out_fasta <- file.path(opt$outdir, paste0(ref_id, ".fasta"))
  out_list <- file.path(opt$outdir, paste0(ref_id, ".list"))
  
  records <- list()

  # 查找ref对应的完整seq名
  ref_seq_name <- find_seq_name(ref_id)
  if (!is.na(ref_seq_name)) {
    records[[ref_seq_name]] <- seqs[[ref_seq_name]]
  } else {
    warning(sprintf("参考序列 %s 未找到，跳过该文件", ref_id))
    next
  }

  query_ids <- unique(ref_to_queries[[ref_id]])
  for (q_id in query_ids) {
    q_seq_name <- find_seq_name(q_id)
    if (!is.na(q_seq_name)) {
      records[[q_seq_name]] <- seqs[[q_seq_name]]
    } else {
      warning(sprintf("查询序列 %s 未找到", q_id))
    }
  }


  if (length(records) > 0) {
    xstring_set <- DNAStringSet(records)
    writeXStringSet(xstring_set, filepath = out_fasta)
    cat(sprintf("写入 %s，共 %d 条序列\n", out_fasta, length(records)))

    # 生成对应的.list文件，写入简化ID
    simple_ids <- sapply(names(records), simplify_id)
    writeLines(simple_ids, con = out_list)
    cat(sprintf("写入 %s，包含 %d 个ID\n", out_list, length(simple_ids)))
  } else {
    warning(sprintf("%s 无有效序列写入", ref_id))
  }

}
