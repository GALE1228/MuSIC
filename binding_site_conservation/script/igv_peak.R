#!/usr/bin/env Rscript
suppressMessages(library(optparse))
suppressMessages(library(dplyr))

option_list <- list(
  make_option(c("-f", "--file1"), type="character",
              help="Input file1 path, filename format RBP_humanTranscript.bed", metavar="FILE1"),
  make_option(c("-b", "--file2"), type="character",
              help="Input file2 path, mapping interval file2, 5 columns: transcript, map_start, map_end, map_len2, some_score", metavar="FILE2"),
  make_option(c("-s", "--species"), type="character",
              help="Input species mapping file path, format: species<TAB>prefix", metavar="SPECIES"),
  make_option(c("-o", "--outdir"), type="character",
              help="Output directory for storing narrowPeak files", metavar="OUTDIR")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$file1) || is.null(opt$file2) || is.null(opt$species) || is.null(opt$outdir)) {
  stop("Please provide --file1, --file2, --species, --outdir parameters.")
}
if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive=TRUE)

# Extract RBP name and human transcript ID
file1_base  <- basename(opt$file1)
file1_noext <- sub("\\.bed$", "", file1_base)
parts       <- strsplit(file1_noext, "_", fixed=TRUE)[[1]]
if (length(parts) < 2) stop("file1 filename should be 'RBP_humanTranscript.bed' format.")
rbp_name <- parts[1]
human_tx <- parts[2]

# Read file1 (HAR data)
file1 <- read.table(opt$file1, sep="\t", header=FALSE, stringsAsFactors=FALSE,
                    col.names=c("transcript","har_start","har_end","map_len","har_score"))
# Read file2 (mapping intervals)
file2 <- read.table(opt$file2, sep="\t", header=FALSE, stringsAsFactors=FALSE,
                    col.names=c("transcript","map_start","map_end","map_len2","some_score"))
# Read species mapping table
spec  <- read.table(opt$species, sep="\t", header=FALSE, stringsAsFactors=FALSE,
                    col.names=c("species","prefix_raw"))

# For each row in file1, match corresponding file2 interval: har_start and har_end must both be within map_start-map_end range
matched_list <- lapply(seq_len(nrow(file1)), function(i) {
  r1 <- file1[i, ]
  f2 <- subset(file2, transcript == r1$transcript)
  sel <- which(r1$har_start >= f2$map_start & r1$har_end <= f2$map_end)
  if (length(sel) == 0) {
    warning(sprintf("Row %d: cannot find matching map interval, skipping: %s %d-%d",
                    i, r1$transcript, r1$har_start, r1$har_end))
    return(NULL)
  }
  f2_sel <- f2[sel[1], ]
  cbind(r1, f2_sel)
})
merged <- do.call(rbind, matched_list)

# Split by species prefix and output narrowPeak files
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
    peak        = sub[, "har_start"] - sub[, "map_start"],  # Adjust to HAR start relative to map_start offset
    stringsAsFactors = FALSE
  )

  # MSA extended fragment correction
  out_df <- out_df %>% 
    mutate (chromStart_old = chromStart,
            chromStart = pmax(chromStart, peak + chromStart - 100),
            chromEnd   = pmin(chromEnd, peak + chromStart_old + 100) -2,
            peak = peak - chromStart + chromStart_old) %>%
    select(-chromStart_old)


  # Output filename: RBP_species.narrowPeak, use append mode
  out_file <- file.path(opt$outdir, paste0(rbp_name, "_", sp, ".narrowPeak"))
  write.table(out_df, file=out_file, sep="\t", quote=FALSE,
              row.names=FALSE, col.names=FALSE, append=TRUE)
}

message("Processing complete, output directory: ", normalizePath(opt$outdir))
