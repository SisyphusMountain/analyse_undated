import os
import pandas as pd
import concurrent.futures
from tqdm import tqdm
from helper import (create_tree, get_induced_donor, make_all_sampled_trees,
                    compute_spr_for_all_pairs, classify_sprs)
from ete3 import Tree

def process_ablated_tree(args):
    """
    Process a single ablated tree (for parallel execution).
    
    Args:
        args: Tuple containing (node_name, ablated_tree_path, output_dir, induced_donor, compute_spr_bin)
        
    Returns:
        DataFrame with SPR classifications for this ablated tree
    """
    print(f"currently processing {args}")
    node_name, ablated_tree_path, output_dir, induced_donor, compute_spr_bin = args
    
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
    return classification_df

def run_complete_analysis(output_dir, num_extant_tips=20, birth_rate=1.0, death_rate=1.0, 
                          seed=42, compute_spr_bin=None, max_workers=None):
    """
    Run a complete analysis pipeline:
    1. Create a tree and save it
    2. For each removable clade, find the induced donor
    3. Create ablated trees for all possible clades
    4. For each ablated tree, compute all possible SPRs (in parallel)
    5. Classify SPRs by topology
    6. Combine all results into a single dataframe
    
    Args:
        output_dir: Directory to store all results
        num_extant_tips: Number of tips for the tree
        birth_rate: Birth rate for the birth-death process
        death_rate: Death rate for the birth-death process
        seed: Random seed
        compute_spr_bin: Path to the SPR computation binary
        max_workers: Maximum number of worker threads (None uses all available cores)
        
    Returns:
        DataFrame with all SPR classifications and the corresponding removed nodes
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
    
    # Prepare arguments for parallel processing
    process_args = []
    for node_name in list_nodes_to_test:
        ablated_tree_path = os.path.join(ablated_trees_dir, f"sampled_tree_{node_name}.newick")
        process_args.append((
            node_name,
            ablated_tree_path,
            output_dir,
            induced_donors.get(node_name, None),
            compute_spr_bin
        ))
    
    # Process ablated trees in parallel with a progress bar
    all_classifications = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks and get future objects
        future_to_node = {executor.submit(process_ablated_tree, arg): arg[0] for arg in process_args}
        
        # Process results as they complete with a progress bar
        for future in tqdm(concurrent.futures.as_completed(future_to_node), 
                          total=len(future_to_node),
                          desc="Processing ablated trees"):
            node_name = future_to_node[future]
            try:
                result = future.result()
                if result is not None and not result.empty:
                    all_classifications.append(result)
            except Exception as e:
                print(f"Error processing ablated tree for node {node_name}: {e}")
    
    # Step 6: Combine all classifications into a single dataframe
    if all_classifications:
        combined_df = pd.concat(all_classifications, ignore_index=True)
        
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
    compute_spr_bin = "./make_spr"
    # Use 4 worker threads - adjust as needed
    run_complete_analysis(output_dir, compute_spr_bin=compute_spr_bin, max_workers=1)
