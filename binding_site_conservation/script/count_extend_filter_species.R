#!/usr/bin/env Rscript
suppressMessages({
  library(optparse)
  library(data.table)
  library(stringr)
})

option_list <- list(
  make_option(c("-i", "--input"),   type="character", default=NULL,
              help="Input TSV file path, no header, at least 5 columns: ID,start,end,length,score", metavar="file"),
  make_option(c("-s", "--species"), type="character", default=NULL,
              help="Species list TSV, first column is species name, second column is ID prefix pattern for filtering [optional]", metavar="file"),
  make_option(c("-e", "--extend"),  type="integer", default=0,
              help="Range extension value, extend N bases left and right before checking overlap [default %default]", metavar="N")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input)) {
  stop("Must specify input file path -i")
}

# 1. Read data
df <- fread(opt$input, header=FALSE, sep="\t", data.table=FALSE,
            col.names = c("ID","start","end","length","score"))

# 2. Filter by species list (if -s is provided)
if (!is.null(opt$species)) {
  sp <- fread(opt$species, header=FALSE, sep="\t", data.table=FALSE,
              col.names = c("species","prefix"))
  # Construct regex: starts with any prefix
  pattern <- paste0("^(", paste0(sp$prefix, collapse="|"), ")")
  df <- df[str_detect(df$ID, pattern), , drop=FALSE]
}

# 3. First extend intervals
df$start_ext <- pmax(0, df$start - opt$extend)
df$end_ext   <- df$end   + opt$extend

# 4. Count unique elements in first column
unique_ids <- unique(df$ID)
cat(sprintf("%s\t%d\t", opt$input, length(unique_ids)))

# 5. Group: ENST starting vs non-ENST starting
enst_rows   <- grepl("^ENST", df$ID)
df_enst     <- df[enst_rows, , drop=FALSE]
df_non_enst <- df[!enst_rows, , drop=FALSE]

# Define function to check if intervals overlap
is_overlap <- function(s1, e1, s2, e2) !(e1 < s2 || e2 < s1)

# 6. Count non-ENST rows that overlap with any ENST row
overlap_non_enst <- vapply(seq_len(nrow(df_non_enst)), function(i) {
  any(mapply(is_overlap,
             df_non_enst$start_ext[i], df_non_enst$end_ext[i],
             df_enst$start_ext,    df_enst$end_ext))
}, logical(1))

# 7. Count unique IDs of overlapping non-ENST rows
ids_overlap_non_enst <- unique(df_non_enst$ID[overlap_non_enst])
cat(sprintf("%d\t", length(ids_overlap_non_enst)))

# 8. For each ENST row, count unique non-ENST row IDs that overlap, take maximum
max_overlap_count <- 0
for (j in seq_len(nrow(df_enst))) {
  st_j <- df_enst$start_ext[j]
  en_j <- df_enst$end_ext[j]
  flag <- mapply(is_overlap,
                 df_non_enst$start_ext, df_non_enst$end_ext,
                 MoreArgs = list(s2 = st_j, e2 = en_j))
  ids_ov <- unique(df_non_enst$ID[flag])
  max_overlap_count <- max(max_overlap_count, length(ids_ov))
}
cat(sprintf("%d\n", max_overlap_count))
