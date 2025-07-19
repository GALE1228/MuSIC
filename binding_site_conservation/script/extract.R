#!/usr/bin/env Rscript
suppressMessages(library(optparse))
suppressMessages(library(Biostrings))
suppressMessages(library(data.table))
suppressMessages(library(stringr))

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Filtered TSV file path, no header, contains Query ID and Ref ID columns", metavar="file"),
  make_option(c("-f", "--fasta"), type="character", default=NULL,
              help="FASTA file path containing all query and ref sequences", metavar="file"),
  make_option(c("-o", "--outdir"), type="character", default="output_fasta_by_ref",
              help="Output directory, default 'output_fasta_by_ref'", metavar="dir")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$input) || is.null(opt$fasta)) {
  cat("Must specify input TSV and FASTA file paths\n")
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

# Create a mapping: simplified ID (first segment after | or space) -> complete seq name
simplify_id <- function(full_id) {
  str_split(full_id, "\\||\\s")[[1]][1]
}

id_map <- setNames(seq_names, sapply(seq_names, simplify_id))

# Find complete seq name by simplified ID
find_seq_name <- function(id) {
  matched_names <- id_map[names(id_map) == id]
  if (length(matched_names) == 0) {
    # Try prefix matching: match seq_name starting with id
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

  # Find complete seq name for ref
  ref_seq_name <- find_seq_name(ref_id)
  if (!is.na(ref_seq_name)) {
    records[[ref_seq_name]] <- seqs[[ref_seq_name]]
  } else {
    warning(sprintf("Reference sequence %s not found, skipping this file", ref_id))
    next
  }

  query_ids <- unique(ref_to_queries[[ref_id]])
  for (q_id in query_ids) {
    q_seq_name <- find_seq_name(q_id)
    if (!is.na(q_seq_name)) {
      records[[q_seq_name]] <- seqs[[q_seq_name]]
    } else {
      warning(sprintf("Query sequence %s not found", q_id))
    }
  }


  if (length(records) > 0) {
    xstring_set <- DNAStringSet(records)
    writeXStringSet(xstring_set, filepath = out_fasta)
    cat(sprintf("Written %s, total %d sequences\n", out_fasta, length(records)))

    # Generate corresponding .list file, write simplified IDs
    simple_ids <- sapply(names(records), simplify_id)
    writeLines(simple_ids, con = out_list)
    cat(sprintf("Written %s, contains %d IDs\n", out_list, length(simple_ids)))
  } else {
    warning(sprintf("%s has no valid sequences to write", ref_id))
  }

}
