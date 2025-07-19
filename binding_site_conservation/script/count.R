#!/usr/bin/env Rscript
suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(stringr))

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input TSV file path, no header, at least 5 columns: ID,start,end,length,score", metavar="file")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input)) {
  stop("Must specify input file path -i")
}

# Read data
df <- fread(opt$input, header=FALSE, sep="\t", data.table=FALSE,
            col.names = c("ID","start","end","length","score"))

# 1. Count unique elements in first column
unique_ids <- unique(df$ID)
cat(sprintf("%s\t%d\t", opt$input, length(unique_ids)))

# 2. Group: ENST starting vs non-ENST starting
enst_rows <- grepl("^ENST", df$ID)
df_enst <- df[enst_rows, , drop=FALSE]
df_non_enst <- df[!enst_rows, , drop=FALSE]

# Define function to check if intervals overlap, intervals [start1,end1] and [start2,end2]
is_overlap <- function(start1, end1, start2, end2) {
  return(!(end1 < start2 || end2 < start1))
}

# Count non-ENST rows that overlap with any ENST row
overlap_non_enst <- logical(nrow(df_non_enst))
for (i in seq_len(nrow(df_non_enst))) {
  non_st <- df_non_enst$start[i]
  non_en <- df_non_enst$end[i]
  # Check if overlaps with any ENST row
  overlap_any <- any(mapply(is_overlap,
                            non_st, non_en,
                            df_enst$start, df_enst$end))
  overlap_non_enst[i] <- overlap_any
}
# Count unique IDs of overlapping non-ENST rows
ids_overlap_non_enst <- unique(df_non_enst$ID[overlap_non_enst])
cat(sprintf("%d\t", length(ids_overlap_non_enst)))

# 3. For each ENST row, count unique non-ENST row IDs that overlap, take maximum
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
