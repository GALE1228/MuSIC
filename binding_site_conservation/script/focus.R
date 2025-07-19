#!/usr/bin/env Rscript

suppressMessages(library(optparse))

# Define acceptable command line options: short and long parameters
option_list <- list(
  make_option(c("-a", "--file1"), type = "character", help = "file1 path (format same as example)", metavar = "FILE1"),
  make_option(c("-b", "--file2"), type = "character", help = "file2 path (format same as example)", metavar = "FILE2"),
  make_option(c("-o", "--output"), type = "character", help = "Output file path", metavar = "OUTPUT")
)

opt_parser <- OptionParser(option_list = option_list, 
                           description = "Usage example: Rscript update_coordinates.R -a file1.txt -b file2.txt -o output.txt")
opt <- parse_args(opt_parser)

# If any required parameter is missing, print help and exit
if (is.null(opt$file1) || is.null(opt$file2) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("Please use --file1/-a, --file2/-b, --output/-o parameters.")
}

file1_path  <- opt$file1
file2_path  <- opt$file2
output_path <- opt$output

# Read file1: four columns, first column idx, third column rel_start, fourth column rel_end
file1 <- read.table(file1_path, header = FALSE, stringsAsFactors = FALSE)

# Read file2: five columns, separated by '|'
file2 <- read.table(file2_path,
                    sep = "|",
                    header = FALSE,
                    stringsAsFactors = FALSE,
                    col.names = c("ID", "region_start", "region_end", "boundary_end", "score"))

# Check if row counts match
if (nrow(file1) != nrow(file2)) {
  stop("file1 and file2 have different row counts, cannot correspond line by line.")
}

# Copy as output
output <- file2

for (i in seq_len(nrow(file1))) {
  idx          <- file1[i, 1]
  row_idx      <- idx + 1
  rel_start    <- as.numeric(file1[i, 3])
  rel_end      <- as.numeric(file1[i, 4])
  region_start <- as.numeric(file2[row_idx, "region_start"])
  boundary_end <- as.numeric(file2[row_idx, "boundary_end"])

  # Calculate new global coordinates
  new_start <- region_start + rel_start - 1
  new_end   <- region_start + rel_end   - 1
  if (new_end > boundary_end) {
    new_end <- boundary_end
  }

  # Update region_start and region_end in output table
  output[row_idx, "region_start"] <- new_start
  output[row_idx, "region_end"]   <- new_end
}

write.table(output,
            file = output_path,
            sep = "|",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
