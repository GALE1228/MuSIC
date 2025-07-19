#!/usr/bin/env Rscript

suppressMessages(library(optparse))
library(stringr)

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input TSV file path, no header, contains id,f,kc,f,lenQ,lenR columns", metavar="file"),
  make_option(c("-o", "--output"), type="character", default="output_maxscore.tsv",
              help="Output file path [default= %default]", metavar="file"),
  make_option(c("-w", "--weight"), type="double", default=0.8,
              help="id*kc weight, length weight is 1-w, range [0,1] [default= %default]")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Must provide input file path --input", call.=FALSE)
}

# 1. Read original data (no header)
df <- read.table(opt$input, header=FALSE, sep="\t", stringsAsFactors=FALSE, fill=TRUE)

# 2. Adjust column numbers based on actual data (example adjusted based on your example)
# Assume:
# V1 = Query transcript ID
# V3 = lenQ (query sequence length)
# V8 = lenR (reference sequence length)
# Last two columns are id:f:0.xxxxxx and kc:f:0.xxxxxx

id_col <- ncol(df) - 1
kc_col <- ncol(df)

# Extract id and kc values
df$id <- as.numeric(str_extract(df[[id_col]], "(?<=id:f:)\\d*\\.?\\d+"))
df$kc <- as.numeric(str_extract(df[[kc_col]], "(?<=kc:f:)\\d*\\.?\\d+"))

# Extract length columns
df$lenQ <- as.numeric(df$V3)
df$lenR <- as.numeric(df$V8)

# 3. Calculate score function
calculate_score_weighted <- function(id, kc, lenQ, lenR, w1) {
  maxLenDiff <- pmax(pmin(2*lenR, 2000), 1)  # Ensure at least 1, avoid division by 0
  norm_len_diff <- pmin(abs(lenQ - lenR), maxLenDiff) / maxLenDiff
  length_factor <- 1 - norm_len_diff
  score <- w1 * (id * kc) + (1 - w1) * length_factor
  return(score)
}


df$score <- mapply(calculate_score_weighted, df$id, df$kc, df$lenQ, df$lenR,
                   MoreArgs = list(w1 = opt$weight))

# 4. Take maximum score record for each transcript pair (Query and Reference ID)
# Assume Query ID in V1, Reference ID in V6
library(dplyr)

df_max <- df %>%
  group_by(RefID = V6) %>%
  slice_max(score, n = 1, with_ties = FALSE) %>%
  ungroup()


# 5. Save results, include all original columns + score
write.table(df_max, file=opt$output, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

cat(sprintf("Processing complete, maximum score results saved to %s, record count %d\n", opt$output, nrow(df_max)))
