# Load required libraries
library(ape)         # For reading Newick trees
library(ggtree)      # For tree visualization
library(ggplot2)     # For plotting
library(dplyr)       # For data manipulation
library(RColorBrewer) # For color palettes

# Function to visualize transfers on a phylogenetic tree
visualize_transfers <- function(tree_file, data_file, node_of_interest) {
  # Read the phylogenetic tree
  tree <- read.tree(tree_file)
  
  # Read the transfer data
  transfers_data <- read.csv(data_file)
  
  # Convert node IDs to character to ensure proper matching
  transfers_data$donor <- as.character(transfers_data$donor)
  transfers_data$receiver <- as.character(transfers_data$receiver)
  transfers_data$removed_node <- as.character(transfers_data$removed_node)
  transfers_data$induced_donor <- as.character(transfers_data$induced_donor)
  
  # Filter data for the specified removed_node
  filtered_data <- transfers_data %>%
    filter(removed_node == as.character(node_of_interest))
  
  if(nrow(filtered_data) == 0) {
    stop(paste("No transfers found for removed_node =", node_of_interest))
  }
  
  # For each topology_id, count how many transfers exist
  topology_counts <- filtered_data %>%
    group_by(topology_id) %>%
    summarize(count = n())
  
  # Merge the counts back to the filtered data
  filtered_data <- filtered_data %>%
    left_join(topology_counts, by = "topology_id")
  
  # Filter for transfers from the induced_donor
  induced_transfers <- filtered_data %>%
    filter(donor == induced_donor)
  
  # Print some debugging info
  cat("Tree nodes:", paste(tree$tip.label, collapse=", "), "\n")
  cat("Internal nodes:", paste(as.character((length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode)), collapse=", "), "\n")
  cat("Filtered data (first few rows):\n")
  print(head(filtered_data))
  cat("Induced transfers (first few rows):\n")
  print(head(induced_transfers))
  
  # Create a base tree plot
  p <- ggtree(tree)
  
  # Get the tree data and node positions
  tree_data <- p$data
  
  # Add tiplab to visualize the tips
  p <- p + geom_tiplab(size=3)
  
  # Add node labels for all nodes including internal ones
  p <- p + 
    geom_text(data = tree_data,
              aes(x = x, y = y, label = label),
              size = 3, vjust = -0.5, hjust = -0.2)
  
  # Create a color palette
  if(nrow(topology_counts) > 0) {
    max_count <- max(topology_counts$count)
    color_palette <- colorRampPalette(brewer.pal(min(9, max_count), "YlOrRd"))(max_count)
  } else {
    color_palette <- brewer.pal(3, "YlOrRd")
  }
  
  # Add arrows for each transfer if there are any induced transfers
  if(nrow(induced_transfers) > 0) {
    # Check if nodes exist in the tree
    valid_transfers <- induced_transfers %>%
      filter(donor %in% c(tree$tip.label, as.character((length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode))),
             receiver %in% c(tree$tip.label, as.character((length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode))))
    
    if(nrow(valid_transfers) > 0) {
      # Add a legend for the counts
      p <- p + 
        scale_color_gradientn(colors = color_palette,
                             name = "Transfer count\nin topology",
                             limits = c(1, max(valid_transfers$count)))
      
      # Add custom arrows
      for(i in 1:nrow(valid_transfers)) {
        transfer <- valid_transfers[i,]
        donor_node <- transfer$donor
        receiver_node <- transfer$receiver
        
        # Get positions for donor and receiver
        donor_pos <- tree_data[tree_data$label == donor_node, c("x", "y")]
        receiver_pos <- tree_data[tree_data$label == receiver_node, c("x", "y")]
        
        # Skip if positions are identical or missing
        if(nrow(donor_pos) == 0 || nrow(receiver_pos) == 0) {
          cat("Warning: Missing position for donor:", donor_node, "or receiver:", receiver_node, "\n")
          next
        }
        
        if(donor_pos$x == receiver_pos$x && donor_pos$y == receiver_pos$y) {
          cat("Warning: Identical positions for donor:", donor_node, "and receiver:", receiver_node, "\n")
          next
        }
        
        # Add a bit of curvature and offset to avoid overlapping with tree lines
        curve_offset <- 0.2 * (max(tree_data$x) - min(tree_data$x))
        
        # Add the arrow
        p <- p + 
          geom_curve(data = data.frame(x = donor_pos$x, y = donor_pos$y, 
                                       xend = receiver_pos$x, yend = receiver_pos$y),
                     aes(x = x, y = y, xend = xend, yend = yend),
                     arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
                     color = color_palette[transfer$count],
                     size = 0.8,
                     alpha = 0.7,
                     curvature = 0.3,
                     ncp = 10) # More control points for smoother curves
      }
    } else {
      cat("Warning: No valid transfers found. Check that donor and receiver nodes exist in the tree.\n")
    }
  } else {
    cat("No transfers from induced_donor found in the filtered data.\n")
  }
  
  # Adjust the plot layout
  p <- p + 
    theme_tree2() +
    labs(title = paste("Transfers with removed_node =", node_of_interest),
         subtitle = "Color intensity indicates number of transfers with same topology_id") +
    theme(plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 12))
  
  # Save the plot
  output_file <- paste0("transfers_removed_node_", node_of_interest, ".pdf")
  ggsave(output_file, p, width = 12, height = 10)
  
  # Return the plot
  return(p)
}

# Command line argument handling
args <- commandArgs(trailingOnly = TRUE)

tree_file <- "results/original_tree.nwk"
data_file <- "results/combined_spr_classifications.csv"

if(length(args) == 0) {
  cat("Usage: Rscript plot_transfers.R <removed_node>\n")
  cat("Example: Rscript plot_transfers.R 10\n")
  
  # Default value for testing
  node_of_interest <- "10"  
} else {
  node_of_interest <- args[1]
}

# Print input parameters
cat("Tree file:", tree_file, "\n")
cat("Data file:", data_file, "\n")
cat("Removed node:", node_of_interest, "\n")

# Run the visualization function
tryCatch({
  result_plot <- visualize_transfers(tree_file, data_file, node_of_interest)
  print(result_plot)  # Display the plot if running interactively
}, error = function(e) {
  cat("Error occurred:", conditionMessage(e), "\n")
  
  # Print more debugging information
  if(file.exists(tree_file)) {
    cat("Tree file exists. Checking contents...\n")
    tree <- tryCatch(read.tree(tree_file), error = function(e) NULL)
    if(!is.null(tree)) {
      cat("Tree successfully read with", length(tree$tip.label), "tips and", tree$Nnode, "internal nodes.\n")
    } else {
      cat("Failed to read tree file. Check if it's in valid Newick format.\n")
    }
  } else {
    cat("Tree file does not exist:", tree_file, "\n")
  }
  
  if(file.exists(data_file)) {
    cat("Data file exists. Checking contents...\n")
    data <- tryCatch(read.csv(data_file), error = function(e) NULL)
    if(!is.null(data)) {
      cat("Data successfully read with", nrow(data), "rows and", ncol(data), "columns.\n")
      cat("Column names:", paste(colnames(data), collapse=", "), "\n")
    } else {
      cat("Failed to read data file. Check if it's a valid CSV.\n")
    }
  } else {
    cat("Data file does not exist:", data_file, "\n")
  }
})