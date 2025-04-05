import os
import pandas as pd
import multiprocessing
from tqdm import tqdm
from helper import (create_tree, get_induced_donor, make_all_sampled_trees,
                    compute_spr_for_all_pairs, classify_sprs)
from ete3 import Tree
import time

def process_ablated_tree(node_name, ablated_trees_dir, output_dir, induced_donor, compute_spr_bin):
    """
    Process a single ablated tree (for parallel execution).
    
    Args:
        node_name: Name of the node that was removed
        ablated_trees_dir: Directory containing ablated trees
        output_dir: Directory to store results
        induced_donor: Induced donor for this node
        compute_spr_bin: Path to the SPR computation binary
        
    Returns:
        DataFrame with SPR classifications for this ablated tree
    """
    start_time = time.time()
    ablated_tree_path = os.path.join(ablated_trees_dir, f"sampled_tree_{node_name}.newick")
    
    # Skip if the ablated tree does not exist
    if not os.path.exists(ablated_tree_path):
        return None
    
    # Create a directory for SPR results for this ablated tree
    spr_dir = os.path.join(output_dir, f"spr_results_{node_name}")
    os.makedirs(spr_dir, exist_ok=True)
    
    # Compute all possible SPRs for this ablated tree
    compute_spr_for_all_pairs(ablated_tree_path, spr_dir, compute_spr_bin)
    
    # Classify SPRs by topology
    classification_df = classify_sprs(spr_dir)
    
    # Add the removed node information
    classification_df["removed_node"] = node_name
    
    # Add the induced donor information
    classification_df["induced_donor"] = induced_donor
    
    # Save individual result to avoid memory issues with large datasets
    result_path = os.path.join(output_dir, f"spr_classification_{node_name}.csv")
    classification_df.to_csv(result_path, index=False)
    
    elapsed = time.time() - start_time
    return {
        "node_name": node_name,
        "elapsed_time": elapsed,
        "num_classifications": len(classification_df) if not classification_df.empty else 0
    }

def run_complete_analysis(output_dir, num_extant_tips=20, birth_rate=1.0, death_rate=1.0, 
                          seed=42, compute_spr_bin=None, max_workers=None):
    """
    Run a complete analysis pipeline with simple process-based parallelism
    """
    if compute_spr_bin is None:
        raise ValueError("Path to the SPR computation binary must be provided.")

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Step 1: Create a tree and save it
    tree_path = os.path.join(output_dir, "original_tree.nwk")
    create_tree(tree_path, num_extant_tips, birth_rate, death_rate, seed)
    
    # Step 2: Find induced donors for each removable clade
    tree = Tree(tree_path, format=1)
    induced_donors = {}
    list_nodes_to_test = []
    
    root_children = tree.get_children()
    for node in tree.traverse():
        if not node.is_root() and node not in root_children:
            list_nodes_to_test.append(node.name)
            induced_donors[node.name] = get_induced_donor(tree_path, node.name)
    
    # Step 3: Create ablated trees for all possible clades
    ablated_trees_dir = os.path.join(output_dir, "ablated_trees")
    os.makedirs(ablated_trees_dir, exist_ok=True)
    make_all_sampled_trees(tree_path, ablated_trees_dir)
    
    # Step 4 & 5: For each ablated tree, compute SPRs and classify them (in parallel)
    print(f"Processing {len(list_nodes_to_test)} ablated trees in parallel...")
    
    # If max_workers not specified, use CPU count - 1 (leave one core for system)
    if max_workers is None:
        max_workers = max(1, multiprocessing.cpu_count() - 1)
    
    print(f"Using {max_workers} worker processes")
    
    # Create tasks with correct induced donors
    tasks = [
        (node, ablated_trees_dir, output_dir, induced_donors.get(node, None), compute_spr_bin) 
        for node in list_nodes_to_test
    ]
    
    all_stats = []
    
    # Process all tasks at once with process pool
    with multiprocessing.Pool(processes=max_workers) as pool:
        # Execute in parallel with progress bar
        for result in tqdm(
            pool.starmap(process_ablated_tree, tasks),
            total=len(tasks),
            desc="Processing ablated trees"
        ):
            if result:
                all_stats.append(result)
    
    # Print statistics
    if all_stats:
        print("\nProcessing statistics:")
        total_time = sum(stat["elapsed_time"] for stat in all_stats)
        total_classifications = sum(stat["num_classifications"] for stat in all_stats)
        
        print(f"Total processing time: {total_time:.2f} seconds")
        print(f"Total classifications: {total_classifications}")
        print(f"Average time per node: {total_time/len(all_stats):.2f} seconds")
    
    # Step 6: Combine all classifications from individual CSV files
    print("Combining all classification results...")
    classification_files = [
        os.path.join(output_dir, f"spr_classification_{node}.csv") 
        for node in list_nodes_to_test
    ]
    
    # Only read files that exist
    existing_files = [f for f in classification_files if os.path.exists(f)]
    
    if existing_files:
        # Read and combine all CSV files
        dfs = []
        for file in tqdm(existing_files, desc="Reading classification files"):
            try:
                df = pd.read_csv(file)
                dfs.append(df)
            except Exception as e:
                print(f"Error reading {file}: {e}")
        
        if dfs:
            combined_df = pd.concat(dfs, ignore_index=True)
            
            # Save the combined results
            output_path = os.path.join(output_dir, "combined_spr_classifications.csv")
            combined_df.to_csv(output_path, index=False)
            
            print(f"Analysis complete. Results saved to {output_path}")
            return combined_df
    
    print("No classifications were found.")
    return pd.DataFrame()  # Return empty dataframe if no classifications were found

if __name__ == "__main__":
    # Example usage
    output_dir = "/home/enzo/Documents/git/ghost_experiments/analyse_undated/results"
    num_extant_tips = 10
    compute_spr_bin = "./make_spr"
    
    # Use process-based parallelism with automatic core detection
    run_complete_analysis(
        output_dir, 
        num_extant_tips=num_extant_tips,
        compute_spr_bin=compute_spr_bin,
    )