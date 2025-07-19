#!/usr/bin/env Rscript

# Load required packages
library(dplyr)
library(ggplot2)
library(optparse)

# Set command line arguments
option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Input file paths, multiple files separated by commas"),
  make_option(c("-o", "--output"), type = "character", default = "forest_plot_output.pdf", help = "Output PDF filename [default: %default]"),
  make_option(c("-c", "--column"), type = "integer", default = 1, help = "Select column for y-axis (1, 2, or 3) [default: %default]")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

# Read input file paths
input_files <- strsplit(opts$input, ",")[[1]]

# Initialize an empty data frame
all_data <- data.frame()

# Define species_name mapping
species_map <- c(
  "species_1" = "PONAB",
  "species_2" = "MACFA",
  "species_3" = "MOUSE",
  "species_4" = "RAT",
  "species_5" = "CHICK",
  "species_6" = "XENLA",
  "species_7" = "DANRE",
  "species_8" = "DROME",
  "species_9" = "ARATH",
  "species_10" = "YEAST"
)

# Process each file
for (file_path in input_files) {
  # Extract filename and file path
  file_name <- basename(file_path)
  
  # Extract grouping information
  group_info <- sub(".*_(species_\\d+).*", "\\1", file_name)  # Get grouping info from filename
  species_number <- sub("species_(\\d+)", "\\1", group_info)  # Extract number from species_x

  # Read file content
  data <- read.table(file_path, header = FALSE)
  colnames(data) <- c("path", "col1", "col2", "col3")  # Name the columns
  
  # Process col1, col2, col3 based on species_number and convert to percentage
  data$col1 <- (data$col1 / (as.numeric(species_number) + 1)) * 100
  data$col2 <- (data$col2 / as.numeric(species_number)) * 100
  data$col3 <- (data$col3 / as.numeric(species_number)) * 100

  # Add RBP name and grouping information to data frame
  data$RBP <- sub("^.*/([A-Za-z0-9_]+)_.*$", "\\1", data$path)
  data$group <- group_info
  
  # Add species_name column
  data$species_name <- species_map[group_info]
  
  data <- data %>%
    group_by(RBP, species_name) %>%
    mutate(mean_y = mean(.data[[paste("col", opts$column, sep = "")]]),
           sorted_col = sort(.data[[paste("col", opts$column, sep = "")]]),
           head_25_mean = mean(head(sorted_col, floor(length(sorted_col) * 0.25))),
           tail_25_mean = mean(tail(sorted_col, floor(length(sorted_col) * 0.25)))) %>%
    ungroup()  # Cancel grouping

  # Bind data to all_data
  all_data <- rbind(all_data, data)
}

# Ensure selected column exists
y_col <- paste("col", opts$column, sep = "")
if (!y_col %in% colnames(all_data)) {
  stop(paste("Column", opts$column, "does not exist in the data."))
}

# Ensure x-axis follows species order
species_order <- c("PONAB", "MACFA", "MOUSE", "RAT", "CHICK", "XENLA", "DANRE", "DROME", "ARATH", "YEAST")
all_data$species_name <- factor(all_data$species_name, levels = species_order)

# Draw forest plot
p <- ggplot(all_data, aes(x = species_name, y = mean_y, color = mean_y)) +
  geom_point(position = position_dodge(width = 0.8), size = 3) +  # Draw mean points
  geom_errorbar(aes(ymin = head_25_mean, ymax = tail_25_mean), width = 0.3, position = position_dodge(width = 0.8)) +  # Error bars from head 25% mean to tail 25% mean
  facet_grid(RBP ~ ., scales = "free_y") +
  scale_color_gradient(low = "#60a9cf", high = "#e8b6ae", limits = c(0, 100)) +  # Fill color gradient based on y mean
  scale_y_continuous(
    limits = c(0, 100),
    breaks = c(0, 50, 100),               # Displayed tick labels
    minor_breaks = c(25, 75)              # Displayed but unlabeled minor tick lines
  ) +
  theme_minimal() +
  labs(
    title = "IGF2BP family of RBPs",  # Chart title in English
    x = "Species (Accumulated)",  # x-axis title in English
    y = "Extent of Conservation of Binding Sites",  # y-axis title updated
    color = "Average Conservancy"  # Legend title in English
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  # Keep four-sided border
    panel.grid.major = element_blank(),  # Cancel major reference lines
    panel.grid.minor = element_blank(),   # Cancel minor reference lines
    panel.spacing = unit(1.2, "lines")
  )

# Output as PDF file
ggsave(opts$output, plot = p, width = 8, height = 6)
