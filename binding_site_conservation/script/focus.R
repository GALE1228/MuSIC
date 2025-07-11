#!/usr/bin/env Rscript

suppressMessages(library(optparse))

# 定义可接受的命令行选项：短参数和长参数
option_list <- list(
  make_option(c("-a", "--file1"), type = "character", help = "file1 路径（格式同示例）", metavar = "FILE1"),
  make_option(c("-b", "--file2"), type = "character", help = "file2 路径（格式同示例）", metavar = "FILE2"),
  make_option(c("-o", "--output"), type = "character", help = "输出文件路径", metavar = "OUTPUT")
)

opt_parser <- OptionParser(option_list = option_list, 
                           description = "用法示例：Rscript update_coordinates.R -a file1.txt -b file2.txt -o output.txt")
opt <- parse_args(opt_parser)

# 如果任意一个必需参数缺失，则打印帮助并退出
if (is.null(opt$file1) || is.null(opt$file2) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("请使用 --file1/-a、--file2/-b、--output/-o 参数。")
}

file1_path  <- opt$file1
file2_path  <- opt$file2
output_path <- opt$output

# 读入 file1：四列，第一列 idx，第三列 rel_start，第四列 rel_end
file1 <- read.table(file1_path, header = FALSE, stringsAsFactors = FALSE)

# 读入 file2：五列，用 '|' 分隔
file2 <- read.table(file2_path,
                    sep = "|",
                    header = FALSE,
                    stringsAsFactors = FALSE,
                    col.names = c("ID", "region_start", "region_end", "boundary_end", "score"))

# 检查行数是否一致
if (nrow(file1) != nrow(file2)) {
  stop("file1 和 file2 行数不一致，无法逐行对应。")
}

# 复制一份作为输出
output <- file2

for (i in seq_len(nrow(file1))) {
  idx          <- file1[i, 1]
  row_idx      <- idx + 1
  rel_start    <- as.numeric(file1[i, 3])
  rel_end      <- as.numeric(file1[i, 4])
  region_start <- as.numeric(file2[row_idx, "region_start"])
  boundary_end <- as.numeric(file2[row_idx, "boundary_end"])

  # 计算新的全局坐标
  new_start <- region_start + rel_start - 1
  new_end   <- region_start + rel_end   - 1
  if (new_end > boundary_end) {
    new_end <- boundary_end
  }

  # 更新输出表的 region_start 和 region_end
  output[row_idx, "region_start"] <- new_start
  output[row_idx, "region_end"]   <- new_end
}

write.table(output,
            file = output_path,
            sep = "|",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
