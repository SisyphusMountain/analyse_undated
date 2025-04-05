# --- Configuration ---
csv_file <- "results/combined_spr_classifications.csv" # Path to your CSV data
tree_file <- "results/original_tree.nwk"                      # Path to your Newick tree file
output_pdf_file <- "all_removed_node_plots_unit_branch.pdf" # Updated PDF name

# --- Load Libraries ---
library(readr)
library(dplyr)
library(ape)
library(ggtree)
library(ggplot2)
library(ggnewscale) # For multiple color scales

# --- Load Data ---
transfer_data <- read_csv(csv_file, col_types = cols(
  donor = col_character(),
  receiver = col_character(),
  topology_id = col_integer(),
  removed_node = col_character(),
  induced_donor = col_character()
))
tree <- read.tree(tree_file)

# <<<--- SET ALL BRANCH LENGTHS TO 1 ---<<<
if (!is.null(tree$edge.length)) {
    cat("Original branch lengths detected. Setting all to 1.\n")
    tree$edge.length <- rep(1, nrow(tree$edge))
} else {
    cat("Tree has no branch lengths. Layout will be cladogram-like by default.\n")
    # Optionally, add unit lengths if none exist, though ggtree handles this
    # tree$edge.length <- rep(1, nrow(tree$edge))
}
# >>>------------------------------------->>>

# --- Get Unique Removed Nodes ---
unique_removed_nodes <- unique(transfer_data$removed_node)
# Optional: Sort nodes if desired, ensure they are character for consistency later
unique_removed_nodes <- sort(as.character(unique_removed_nodes))
cat("Found unique removed_node values:", paste(unique_removed_nodes, collapse=", "), "\n")

# --- List to store plots ---
plot_list <- list()

# --- Loop Through Each Removed Node ---
for (target_removed_node in unique_removed_nodes) {

  cat("\n--- Processing for Removed Node:", target_removed_node, "---\n")

  # Use tryCatch to continue loop even if one plot fails
  tryCatch({

    # 1. Filter the dataframe for the current target removed_node
    sub_df <- transfer_data %>%
      filter(removed_node == as.character(target_removed_node)) # Ensure comparison is character

    if (nrow(sub_df) == 0) {
      cat("Skipping removed_node =", target_removed_node, "- No data found in initial filter.\n")
      next # Skip to the next iteration
    }
    # cat("Dimension of sub_df (filtered by removed_node):", dim(sub_df), "\n") # Less verbose

    # 2. Calculate counts based on BOTH topology_id AND receiver within the FULL filtered subset (sub_df)
    topology_receiver_counts <- sub_df %>%
      group_by(topology_id, receiver) %>%
      summarise(transfer_count = n(), .groups = 'drop')

    # 3. Prepare the data SPECIFICALLY for plotting:
    plot_transfers <- sub_df %>%
      filter(donor == induced_donor) %>%
      select(topology_id, induced_donor, receiver) %>%
      left_join(topology_receiver_counts, by = c("topology_id", "receiver")) %>%
      filter(!is.na(transfer_count)) %>%
      rename(from_label = induced_donor, to_label = receiver) %>%
      distinct(from_label, to_label, topology_id, transfer_count, .keep_all = TRUE)

    # cat("Dimension of final plot_transfers data:", dim(plot_transfers), "\n") # Less verbose

    # --- Plotting ---

    # 1. Create the basic tree plot object (using tree with modified lengths)
    p_base <- ggtree(tree) +
           geom_tiplab(size = 3, align = TRUE, linesize = 0.2) + # Added linesize
           geom_text2(aes(subset = !isTip, label = label), hjust = -0.2, size = 2.5)

    # 2. Calculate the desired x-limit (might need adjustment for unit lengths)
    max_x_coord <- max(p_base$data$x, na.rm = TRUE)
    # Extend slightly more to accommodate labels and arrows with unit branches
    upper_x_limit <- max_x_coord + (0.6 * max_x_coord) # Extend by 60% of max depth

    # 3. Add the x-axis limit modification to the base plot
    p <- p_base + xlim(NA, upper_x_limit)


    # --- Add Red Color for Removed Clade ---
    tree_plot_data <- p$data # Get data for modification

    # Find the node number for the target removed node label (tip or internal)
    removed_node_info <- tree_plot_data %>%
                         filter(label == as.character(target_removed_node)) %>% # Ensure comparison is character
                         select(node, isTip)

    nodes_to_color_red <- c()
    removed_node_is_tip <- FALSE

    if (nrow(removed_node_info) == 0) {
        cat("Warning: Target removed node '", target_removed_node, "' not found as any label in the tree.\n")
        removed_node_num <- NA
    } else {
        removed_node_num <- removed_node_info$node[1]
        removed_node_is_tip <- removed_node_info$isTip[1]

        parent_edge_row <- which(tree$edge[, 2] == removed_node_num)
        if (length(parent_edge_row) > 0) {
            parent_node_num <- tree$edge[parent_edge_row, 1]
            nodes_to_color_red <- c(nodes_to_color_red, removed_node_num)
        } else {
            parent_node_num <- NA
            cat("Warning: Could not find parent edge for removed node: ", removed_node_num, "\n")
        }

        if (!removed_node_is_tip) {
            direct_children <- tree$edge[tree$edge[, 1] == removed_node_num, 2]
            nodes_to_color_red <- c(nodes_to_color_red, direct_children) %>% unique()
        }
    }

    tree_plot_data <- tree_plot_data %>% mutate(clade_color = "black")

    if (length(nodes_to_color_red) > 0) {
        # cat("Attempting to color incoming edges red for nodes:", paste(nodes_to_color_red, collapse=", "), "\n") # Less verbose
        tree_plot_data <- tree_plot_data %>%
            mutate(clade_color = ifelse(node %in% nodes_to_color_red, "red", clade_color))
    } else {
        cat("No nodes identified for coloring red for removed_node =", target_removed_node, "\n")
    }

    p$data <- tree_plot_data # Update plot object's data

    # Add the tree layer using the modified data and mapping color
    p <- p + geom_tree(aes(color = clade_color), linetype = "solid", size=0.5) + # Added size
        scale_color_identity(guide = "none") # Use literal colors for TREE


    # Get node coordinates (use the updated p$data)
    node_coords <- p$data %>%
                   select(node, label, isTip, x, y) %>%
                   filter(!is.na(label))

    # Map transfer labels to node numbers and coordinates
    plot_data_coords <- plot_transfers %>%
      left_join(node_coords, by = c("from_label" = "label")) %>%
      rename(from_node = node, x_start = x, y_start = y) %>%
      select(-isTip) %>%
      left_join(node_coords, by = c("to_label" = "label")) %>%
      rename(to_node = node, x_end = x, y_end = y) %>%
      select(-isTip) %>%
      filter(!is.na(from_node), !is.na(to_node))

    # Check if there are any transfers left to plot for this specific node
    if (nrow(plot_data_coords) > 0) {
      # cat("Plotting", nrow(plot_data_coords), "transfers.\n") # Less verbose
      p <- p + ggnewscale::new_scale_color() # Introduce new color scale for arrows

      p <- p +
        geom_segment(data = plot_data_coords,
                     aes(x = x_start, y = y_start, xend = x_end, yend = y_end,
                         color = factor(transfer_count), # Arrow color based on count
                         linetype = ifelse(transfer_count == 1, "dashed", "solid")), # Arrow linetype based on count
                     arrow = arrow(length = unit(0.1, "inches"), type = "closed", angle = 20), # Smaller arrow
                     linewidth = 0.6) + # Thinner arrows
        # Use a discrete color scale for the ARROWS
        scale_color_viridis_d(option = "D", name = "Transfers") + # Simpler legend title
        # Use linetype values directly for arrows and hide its legend
        scale_linetype_identity(guide = "none")

    } else {
      cat("No valid transfers found to plot for removed_node =", target_removed_node, "\n")
      # Add dummy scale to avoid potential errors later if no arrows drawn
      p <- p + scale_linetype_identity(guide = "none")
    }

    # Add title and theme elements
    p <- p +
      ggtitle(paste("Removed Node:", target_removed_node)) + # Shorter title
      theme_tree2() + # Use a minimal tree theme
      theme(plot.title = element_text(hjust = 0.5, size=10), # Smaller title
            legend.position = "right",
            legend.key.size = unit(0.4, 'cm'), # Smaller legend keys
            legend.title=element_text(size=8),  # Smaller legend title
            legend.text=element_text(size=7))   # Smaller legend text

    # --- Store the plot ---
    plot_list[[paste0("node_", target_removed_node)]] <- p
    cat("Plot generated for removed_node =", target_removed_node, "\n")

  }, error = function(e) {
    # Handle errors during plot generation for a specific node
    cat("ERROR generating plot for removed_node =", target_removed_node, ":\n")
    print(e) # Print the error message
    cat("Skipping this node and continuing...\n")
  }) # End tryCatch

} # --- End loop through removed nodes ---


# --- Save all plots to a single PDF ---
if (length(plot_list) > 0) {
  cat("\n--- Saving", length(plot_list), "plots to PDF:", output_pdf_file, "---\n")
  # Use ggsave for multi-page PDF
  ggsave(
       filename = output_pdf_file,
       plot = gridExtra::marrangeGrob(grobs = plot_list, nrow=1, ncol=1),
       width = 7, height = 7 # Adjust dimensions as needed
    )

  # Alternative using base R pdf device (less flexible page size handling maybe)
  # pdf(output_pdf_file, width = 7, height = 7) # Adjust dimensions
  # for (plot_name in names(plot_list)) {
  #   # cat("Printing plot:", plot_name, "\n") # Less verbose
  #   print(plot_list[[plot_name]])
  # }
  # dev.off() # Close the PDF device

  cat("--- PDF saving complete ---\n")
} else {
  cat("\n--- No plots were generated successfully to save to PDF. ---\n")
}

cat("\n--- Script Finished ---\n")