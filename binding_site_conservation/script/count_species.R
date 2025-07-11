#!/usr/bin/env Rscript
suppressMessages({
  library(optparse)
  library(data.table)
  library(stringr)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="待统计的BED/TSV文件路径，第一列为ID", metavar="file"),
  make_option(c("-m", "--map"), type="character", default="id_species.txt",
              help="物种ID前缀映射文件，格式：<Species> <ID_prefixs>，用空白分隔，无表头，多个前缀用“/”分隔", metavar="file"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="输出文件路径（若不指定，则打印到标准输出）", metavar="file"),
  make_option(c("--header"), action="store_true", default=FALSE,
              help="是否输出表头 [default: 不输出]")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input)) {
  stop("必须指定输入文件路径：-i")
}

# 1. 读取 id_species 映射文件，用 read.table 读取任意空白分隔
map_df <- read.table(opt$map, header=FALSE, sep="", stringsAsFactors=FALSE,
                     col.names = c("Species", "Prefixes"))

# 2. 读取待统计文件，提取第一列 ID
bed_df <- fread(opt$input, header=FALSE, sep="\t", data.table=FALSE)
ids <- bed_df[[1]]

# 3. 对每个物种，根据可能的多个前缀判断是否存在匹配
results <- integer(nrow(map_df))
for (j in seq_len(nrow(map_df))) {
  # 拆分多个前缀
  raw_prefixes <- str_trim(map_df$Prefixes[j])
  prefix_list <- unlist(strsplit(raw_prefixes, "/", fixed=TRUE))
  # 构建正则，例如 "^(NM|NR|XM|XR)"
  pattern <- paste0("^(", paste(prefix_list, collapse="|"), ")")
  present <- any(str_detect(ids, pattern))
  results[j] <- ifelse(present, 1L, 0L)
}

# 4. 准备输出：第一列为输入文件名，后续列依次为各物种对应的0/1
file_base <- basename(opt$input)
out_df <- data.frame(
  File = file_base,
  t(results),
  stringsAsFactors = FALSE
)
colnames(out_df) <- c("File", map_df$Species)

# 5. 写出或打印，依据 --header 参数决定是否输出列名
if (!is.null(opt$output)) {
  fwrite(out_df, file = opt$output, sep = "\t", quote = FALSE,
         row.names = FALSE, col.names = opt$header)
} else {
  fwrite(out_df, file = "", sep = "\t", quote = FALSE,
         row.names = FALSE, col.names = opt$header)
}
