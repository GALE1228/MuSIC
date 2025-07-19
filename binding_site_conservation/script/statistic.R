#!/usr/bin/env Rscript
suppressMessages({
  library(optparse)
  library(Biostrings)
  library(data.table)
  library(stringr)
})

option_list <- list(
  make_option(c("-b", "--bed"), type="character", default=NULL,
              help="Input bed file, format ID|start|end|..., no header", metavar="file"),
  make_option(c("-f", "--fasta"), type="character", default=NULL,
              help="Multi-sequence alignment FASTA file path", metavar="file"),
  make_option(c("-o", "--outbed"), type="character", default="updated.bed",
              help="Output updated bed file path", metavar="file")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$bed) || is.null(opt$fasta)) {
  stop("Must specify both -b input bed file and -f multi-sequence alignment FASTA file")
}

# Read bed file
bed_raw <- fread(opt$bed, sep="|", header=FALSE, col.names=c("ID", "start", "end", "extra1", "extra2"))
cat(sprintf("Read %d BED entries\n", nrow(bed_raw)))

# Read multi-sequence alignment FASTA as DNAMultipleAlignment object
msa <- readDNAMultipleAlignment(opt$fasta, format="fasta")
seq_names <- rownames(msa)
cat(sprintf("Read %d multi-sequence alignment sequences\n", length(seq_names)))

# Gap-aware coordinate mapping function
get_msa_pos <- function(seq, orig_pos) {
  count <- 0
  for (i in seq_along(seq)) {
    if (substr(as.character(seq), i, i) != "-") count <- count + 1
    if (count == orig_pos) return(i)
  }
  return(NA)
}

# Simple prefix matching to find alignment sequence name
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
    warning(sprintf("Alignment sequence %s not found, keep original coordinates", id))
    next
  }
  idx <- which(rownames(msa) == msa_id)
  if (length(idx) == 0) {
    warning(sprintf("Alignment sequence index %s not found, keep original coordinates", msa_id))
    next
  }
  seq_set <- msa@unmasked[idx]
  seq <- seq_set[[1]]
  start_orig <- bed_raw$start[i]
  end_orig <- bed_raw$end[i]
  start_msa <- get_msa_pos(seq, start_orig)
  end_msa <- get_msa_pos(seq, end_orig)
  if (is.na(start_msa) || is.na(end_msa)) {
    warning(sprintf("Coordinate conversion failed for ID %s, keep original coordinates", id))
    next
  }
  bed_updated$start[i] <- start_msa
  bed_updated$end[i] <- end_msa
}

fwrite(bed_updated, file=opt$outbed, sep="\t", col.names=FALSE)
cat(sprintf("Updated bed file saved to: %s\n", opt$outbed))
