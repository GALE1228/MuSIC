#!/usr/bin/env Rscript
suppressMessages({
  library(optparse)
  library(data.table)
  library(stringr)
  library(Biostrings)
})

option_list <- list(
  make_option(c("-u", "--updated"), type="character", default=NULL,
              help="Updated BED/TSV file (based on multi-sequence alignment coordinates), at least 5 columns: ID, start_upd, end_upd, extra1, extra2", metavar="file"),
  make_option(c("-r", "--raw"), type="character", default=NULL,
              help="Original BED file (before convert.R), format 'ID|start_raw|end_raw|extra1|extra2', no header", metavar="file"),
  make_option(c("-f", "--fasta"), type="character", default=NULL,
              help="Multi-sequence alignment FASTA file path for coordinate restoration", metavar="file"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output file path, only includes original BED entries overlapping with best ENST entries (sorted by ID and coordinates)", metavar="file")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$updated) || is.null(opt$raw) || is.null(opt$fasta)) {
  stop("Must specify -u updated file, -r original BED file, -f FASTA")
}

# 1. Read updated BED/TSV, assume no header, tab separated
df_upd <- fread(opt$updated, header=FALSE, sep="\t", data.table=FALSE,
                col.names = c("ID", "start_upd", "end_upd", "extra1", "extra2"))
if (!all(c("ID","start_upd","end_upd") %in% colnames(df_upd))) {
  stop("Updated file must contain at least three columns: ID, start_upd, end_upd")
}

# 2. Read original BED, assume "|" separated, no header
raw_df <- fread(opt$raw, header=FALSE, sep="|", data.table=FALSE,
                col.names = c("ID", "start_raw", "end_raw", "extra1", "extra2"))
if (nrow(df_upd) != nrow(raw_df)) {
  stop("Updated file and original BED must have same number of rows, ensure one-to-one correspondence")
}

# 3. Load MSA for coordinate restoration
msa <- readDNAMultipleAlignment(opt$fasta, format="fasta")
seq_names <- rownames(msa)
if (length(seq_names) == 0) {
  stop("Cannot read any sequences from FASTA")
}

# Helper: MSA reverse mapping to original coordinates
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

# 4. Separate ENST and non-ENST
is_enst <- grepl("^ENST", df_upd$ID)
df_enst <- df_upd[is_enst, , drop=FALSE]
df_non  <- df_upd[!is_enst, , drop=FALSE]

# 5. Interval overlap check
is_overlap <- function(s1,e1,s2,e2) !(e1 < s2 || e2 < s1)

# 6. Calculate overlap count for each ENST with different non-ENST
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

# 7. Find maximum overlap count, select all ENST rows reaching that value
max_overlap <- if (n_enst > 0) max(overlap_count) else 0
if (max_overlap == 0) {
  stop("No ENST rows overlap with any non-ENST rows")
}
best_enst_rows <- enst_indices[overlap_count == max_overlap]

# 8. For each best ENST row, collect that row and all overlapping non-ENST rows
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

# 9. Restore original coordinates line by line and collect output
out_list <- list()
for (idx in all_indices) {
  row_upd <- df_upd[idx, ]
  id      <- row_upd$ID
  msa_id  <- find_msa_id(id)
  if (is.na(msa_id)) {
    warning(sprintf("MSA sequence %s not found, skipping", id))
    next
  }
  msa_idx <- which(seq_names == msa_id)
  if (length(msa_idx) == 0) {
    warning(sprintf("Sequence %s not found in MSA, skipping", msa_id))
    next
  }
  msa_seq <- msa@unmasked[msa_idx][[1]]
  orig_st <- get_orig_pos(msa_seq, row_upd$start_upd)
  orig_en <- get_orig_pos(msa_seq, row_upd$end_upd)
  if (is.na(orig_st) || is.na(orig_en)) {
    warning(sprintf("Coordinate mapping failed for ID %s, skipping", id))
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
  stop("No output entries generated")
}
out_df <- rbindlist(out_list)

# 10. Sort by ID alphabetically and start ascending
out_df <- out_df[order(out_df$ID, out_df$start), ]

# 11. Write results
if (!is.null(opt$output)) {
  fwrite(out_df, file = opt$output, sep="|", quote=FALSE, col.names=FALSE, row.names=FALSE)
  cat(sprintf("Output %d sorted original BED entries to %s\n", nrow(out_df), opt$output))
} else {
  fwrite(out_df, file="", sep="|", quote=FALSE, col.names=FALSE, row.names=FALSE)
}
