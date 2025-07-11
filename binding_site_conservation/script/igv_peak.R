#!/usr/bin/env Rscript
suppressMessages(library(optparse))
suppressMessages(library(dplyr))

option_list <- list(
  make_option(c("-f", "--file1"), type="character",
              help="输入 file1 路径，文件名格式 RBP_humanTranscript.bed", metavar="FILE1"),
  make_option(c("-b", "--file2"), type="character",
              help="输入 file2 路径，mapping 区间 文件2，5 列：transcript, map_start, map_end, map_len2, some_score", metavar="FILE2"),
  make_option(c("-s", "--species"), type="character",
              help="输入 species 映射文件路径，格式: species<TAB>prefix", metavar="SPECIES"),
  make_option(c("-o", "--outdir"), type="character",
              help="输出目录，用于存放 narrowPeak 文件", metavar="OUTDIR")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$file1) || is.null(opt$file2) || is.null(opt$species) || is.null(opt$outdir)) {
  stop("请提供 --file1, --file2, --species, --outdir 参数。")
}
if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive=TRUE)

# 提取 RBP 名称和 human 转录本 ID
file1_base  <- basename(opt$file1)
file1_noext <- sub("\\.bed$", "", file1_base)
parts       <- strsplit(file1_noext, "_", fixed=TRUE)[[1]]
if (length(parts) < 2) stop("file1 文件名应为 'RBP_humanTranscript.bed' 格式。")
rbp_name <- parts[1]
human_tx <- parts[2]

# 读取文件1（HAR 数据）
file1 <- read.table(opt$file1, sep="\t", header=FALSE, stringsAsFactors=FALSE,
                    col.names=c("transcript","har_start","har_end","map_len","har_score"))
# 读取文件2（mapping 区间）
file2 <- read.table(opt$file2, sep="\t", header=FALSE, stringsAsFactors=FALSE,
                    col.names=c("transcript","map_start","map_end","map_len2","some_score"))
# 读取 species 映射表
spec  <- read.table(opt$species, sep="\t", header=FALSE, stringsAsFactors=FALSE,
                    col.names=c("species","prefix_raw"))

# 对每行 file1，匹配对应 file2 区间：har_start 和 har_end 都必须在 map_start-map_end 范围内
matched_list <- lapply(seq_len(nrow(file1)), function(i) {
  r1 <- file1[i, ]
  f2 <- subset(file2, transcript == r1$transcript)
  sel <- which(r1$har_start >= f2$map_start & r1$har_end <= f2$map_end)
  if (length(sel) == 0) {
    warning(sprintf("行 %d: 找不到匹配的 map 区间，跳过: %s %d-%d",
                    i, r1$transcript, r1$har_start, r1$har_end))
    return(NULL)
  }
  f2_sel <- f2[sel[1], ]
  cbind(r1, f2_sel)
})
merged <- do.call(rbind, matched_list)

# 按物种前缀拆分并输出 narrowPeak 文件
for (i in seq_len(nrow(spec))) {
  sp       <- spec$species[i]
  prefixes <- unlist(strsplit(spec$prefix_raw[i], "/", fixed=TRUE))
  sub      <- merged[sapply(merged[, "transcript"], function(t) any(startsWith(t, prefixes))), ]
  if (nrow(sub) == 0) next

  out_df <- data.frame(
    chrom       = human_tx,
    chromStart  = sub[, "map_start"],
    chromEnd    = sub[, "map_end"],
    name        = rbp_name,
    score       = round(sub[, "har_score"] * 1000),
    strand      = ".",
    signalValue = sub[, "har_score"],
    pValue      = sub[, "har_score"],
    qValue      = sub[, "har_score"],
    peak        = sub[, "har_start"] - sub[, "map_start"],  # 调整为 HAR start 相对 map_start 的偏移
    stringsAsFactors = FALSE
  )

  # MSA延长的fragment校正
  out_df <- out_df %>% 
    mutate (chromStart_old = chromStart,
            chromStart = pmax(chromStart, peak + chromStart - 100),
            chromEnd   = pmin(chromEnd, peak + chromStart_old + 100) -2,
            peak = peak - chromStart + chromStart_old) %>%
    select(-chromStart_old)


  # 输出文件名：RBP_物种.narrowPeak，使用追加模式
  out_file <- file.path(opt$outdir, paste0(rbp_name, "_", sp, ".narrowPeak"))
  write.table(out_df, file=out_file, sep="\t", quote=FALSE,
              row.names=FALSE, col.names=FALSE, append=TRUE)
}

message("处理完成，输出目录：", normalizePath(opt$outdir))
