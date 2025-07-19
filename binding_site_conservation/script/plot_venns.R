# Load required packages
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite")
}
if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  install.packages("VennDiagram")
}
if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse")
}
if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(jsonlite)
library(VennDiagram)
library(optparse)
library(patchwork)
library(ggplot2)
library(grid)

# Define command line arguments
option_list <- list(
  make_option(c("-j", "--json_file"), type = "character", default = "/data1/RiboGen/data/filter_specie_rbp.json",
              help = "JSON file path, default is %default"),
  make_option(c("-o", "--output_file"), type = "character", default = "species_vs_human_venn.pdf",
              help = "Venn diagram output file path, default is %default")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Check if JSON file exists
if (file.exists(opt$json_file)) {
  # Read JSON file
  json_data <- read_json(opt$json_file)
  message("Successfully read JSON file")
  
  # Get human RBP data
  human_rbps <- json_data$human
  
  # Get other species names except human
  other_species <- names(json_data)[names(json_data) != "human"]
  
  # Initialize list to store Venn diagrams
  venn_plots <- list()
  
  # Iterate through each species except human
  for (i in seq_along(other_species)) {
    species <- other_species[i]
    species_rbps <- json_data[[species]]
    
    # Draw Venn diagram
    venn_grob <- draw.pairwise.venn(
      area1 = length(human_rbps),
      area2 = length(species_rbps),
      cross.area = length(intersect(human_rbps, species_rbps)),
      category = c("Human", species),
      fill = c("lightblue", "lightgreen"),
      alpha = 0.5,
      cat.pos = c(0, 0),
      cat.dist = c(0.025, 0.025),
      scaled = FALSE,
      cat.fontface = 4,
      cat.fontfamily = "serif",
      cat.col = c("darkblue", "darkgreen"),
      fontfamily = "serif",
      lty = "blank",
      cex = 1.5
    )
    
    # Convert gList to single grob object
    single_venn_grob <- grid::gTree(children = venn_grob)
    
    # Convert Venn diagram to ggplot object
    venn_ggplot <- ggplot() + 
      annotation_custom(grob = single_venn_grob) +
      theme_void() +
      theme(
        plot.margin = margin(10, 10, 10, 10),
        panel.spacing = unit(2, "cm")  # Increase spacing between subplots
      )
    
    venn_plots[[i]] <- venn_ggplot
  }
  
  # Dynamically calculate number of columns, can be adjusted as needed
  ncol_value <- 2
  # Calculate number of rows
  nrow_value <- ceiling(length(venn_plots) / ncol_value)
  
  # Use patchwork to combine plots
  combined_plot <- wrap_plots(venn_plots, ncol = ncol_value, nrow = nrow_value) + 
    plot_layout(
      guides = "collect",
      widths = rep(1, ncol_value)  # Adjust column widths
    ) +
    theme(
      panel.spacing.x = unit(3, "cm"),  # Increase spacing between columns
      panel.spacing.y = unit(1, "cm")   # Increase spacing between rows
    )
  
  # Set PDF output file
  pdf(opt$output_file)
  
  # Output combined plot
  print(combined_plot)
  
  # Close PDF device
  dev.off()
  message(sprintf("Venn diagram saved as %s", opt$output_file))
} else {
  stop(sprintf("JSON file %s does not exist", opt$json_file), call. = FALSE)
}
