#!/usr/bin/env Rscript
suppressMessages({
  library(optparse)
  library(data.table)
  library(stringr)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="BED/TSV file path to be counted, first column is ID", metavar="file"),
  make_option(c("-m", "--map"), type="character", default="id_species.txt",
              help="Species ID prefix mapping file, format: <Species> <ID_prefixs>, space separated, no header, multiple prefixes separated by '/'", metavar="file"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output file path (if not specified, print to standard output)", metavar="file"),
  make_option(c("--header"), action="store_true", default=FALSE,
              help="Whether to output header [default: no output]")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input)) {
  stop("Must specify input file path: -i")
}

# 1. Read id_species mapping file, use read.table to read any whitespace separated
map_df <- read.table(opt$map, header=FALSE, sep="", stringsAsFactors=FALSE,
                     col.names = c("Species", "Prefixes"))

# 2. Read file to be counted, extract first column ID
bed_df <- fread(opt$input, header=FALSE, sep="\t", data.table=FALSE)
ids <- bed_df[[1]]

# 3. For each species, determine if there are matches based on possible multiple prefixes
results <- integer(nrow(map_df))
for (j in seq_len(nrow(map_df))) {
  # Split multiple prefixes
  raw_prefixes <- str_trim(map_df$Prefixes[j])
  prefix_list <- unlist(strsplit(raw_prefixes, "/", fixed=TRUE))
  # Build regex, e.g., "^(NM|NR|XM|XR)"
  pattern <- paste0("^(", paste(prefix_list, collapse="|"), ")")
  present <- any(str_detect(ids, pattern))
  results[j] <- ifelse(present, 1L, 0L)
}

# 4. Prepare output: first column is input filename, subsequent columns are 0/1 for each species
file_base <- basename(opt$input)
out_df <- data.frame(
  File = file_base,
  t(results),
  stringsAsFactors = FALSE
)
colnames(out_df) <- c("File", map_df$Species)

# 5. Write or print, decide whether to output column names based on --header parameter
if (!is.null(opt$output)) {
  fwrite(out_df, file = opt$output, sep = "\t", quote = FALSE,
         row.names = FALSE, col.names = opt$header)
} else {
  fwrite(out_df, file = "", sep = "\t", quote = FALSE,
         row.names = FALSE, col.names = opt$header)
}
