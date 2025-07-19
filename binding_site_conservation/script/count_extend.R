#!/usr/bin/env Rscript
suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(stringr))

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input TSV file path, no header, at least 5 columns: ID,start,end,length,score", metavar="file"),
  make_option(c("-e", "--extend"), type="integer", default=0,
              help="Range extension value, extend N bases left and right before checking overlap [default: %default]", metavar="N")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input)) {
  stop("Must specify input file path -i")
}

# Read data
df <- fread(opt$input, header=FALSE, sep="\t", data.table=FALSE,
            col.names = c("ID","start","end","length","score"))

# First extend intervals
df$start_ext <- pmax(0, df$start - opt$extend)
df$end_ext   <- df$end   + opt$extend

# 1. Count unique elements in first column
unique_ids <- unique(df$ID)
cat(sprintf("%s\t%d\t", opt$input, length(unique_ids)))

# 2. Group: ENST starting vs non-ENST starting
enst_rows    <- grepl("^ENST", df$ID)
df_enst      <- df[enst_rows, , drop=FALSE]
df_non_enst  <- df[!enst_rows, , drop=FALSE]

# Define function to check if intervals overlap
is_overlap <- function(start1, end1, start2, end2) {
  !(end1 < start2 || end2 < start1)
}

# Count non-ENST rows that overlap with any ENST row
overlap_non_enst <- vapply(seq_len(nrow(df_non_enst)), function(i) {
  any(mapply(is_overlap,
             df_non_enst$start_ext[i], df_non_enst$end_ext[i],
             df_enst$start_ext,    df_enst$end_ext))
}, logical(1))

# Count unique IDs of overlapping non-ENST rows
ids_overlap_non_enst <- unique(df_non_enst$ID[overlap_non_enst])
cat(sprintf("%d\t", length(ids_overlap_non_enst)))

# 3. For each ENST row, count unique non-ENST row IDs that overlap, take maximum
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
