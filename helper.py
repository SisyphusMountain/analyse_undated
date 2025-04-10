from dendropy.model.birthdeath import birth_death_tree as bdt
from ete3 import Tree
import random
import os
import subprocess
from collections import defaultdict
import pandas as pd


def create_tree(tree_path, num_extant_tips, birth_rate=1.0, death_rate=1.0, seed=0):
    """
    Create a tree using the birth-death process and save it to a file.
    """
    rng = random.Random(seed)
    tree = bdt(birth_rate=birth_rate, death_rate=death_rate, num_extant_tips=num_extant_tips, rng=rng)
    tree.write(path=tree_path, schema="newick", suppress_rooting=True)
    # open the tree with ete3
    tree = Tree(tree_path, format=1)
    # Rename tree nodes to have unique names
    counter = 0
    for node in tree.traverse():
        node.name = str(counter)
        counter += 1
    # save the tree
    tree.write(format=1, outfile=tree_path, format_root_node=True)
    return tree

def detach_group(tree, node_name):
    """Detach a node and its descendants from the tree."""
    node = tree.search_nodes(name=node_name)[0]
    node.detach()
    for n in tree.traverse():
        if len(n.children) == 1:
            join_branch(n)
    return tree

def join_branch(node):
    """Substitute node for its only child."""
    if len(node.children) != 1:
        return
        
    child = node.children[0]
    
    if hasattr(node, "dist"):
        child.dist = (child.dist or 0) + node.dist
    
    parent = node.up
    if parent:
        idx = parent.children.index(node)
        parent.children[idx] = child
        child.up = parent

def get_sampled_tree(tree_path, node_name):
    """Create a tree with the specified node removed."""
    tree = Tree(tree_path, format=1)
    # if the node is a child of the root, the result is just the sister node
    left_root_child, right_root_child = tree.get_children()
    if left_root_child.name == node_name:
        return right_root_child
    elif right_root_child.name == node_name:
        return left_root_child
    # otherwise, detach the node and its descendants
    tree = detach_group(tree, node_name)
    for node in tree.traverse():
        if len(node.children) == 1:
            join_branch(node)
    return tree

def make_all_sampled_trees(tree_path, output_folder):
    tree = Tree(tree_path, format=1)
    list_nodes_to_test = []
    # root_children = tree.get_children()
    for node in tree.traverse():
        if not node.is_root():
            list_nodes_to_test.append(node.name)
    
    for node in list_nodes_to_test:
        sampled_tree = get_sampled_tree(tree_path, node)
        sampled_tree_path = os.path.join(output_folder, f"sampled_tree_{node}.newick")
        sampled_tree.write(format=1, outfile=sampled_tree_path, format_root_node=True)

def get_induced_donor(tree_path, clade_name):
    """Find the induced donor after removing a clade."""
    tree = Tree(tree_path, format=1)
    node = tree.search_nodes(name=clade_name)[0]
    return node.get_sisters()[0].name if node.get_sisters() else None

def compute_ghost_species(tree_path, node_name):
    list_species = []
    tree = Tree(tree_path, format=1)
    list_species.append(node_name)
    node = tree.search_nodes(name=node_name)[0]
    for descendant in node.get_descendants():
        list_species.append(descendant.name)
    return list(set(list_species))

def compute_contemporaneity(tree_path, output_dir, compute_contemporaneity_bin):
    """Compute contemporaneity using the external binary."""
    output_path = os.path.join(output_dir, "contemporaneity.csv")
    command = [compute_contemporaneity_bin, tree_path, output_path]
    result = subprocess.run(command, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Error computing contemporaneity: {result.stderr}")
        return None
        
    # Read the output file
    df = pd.read_csv(output_path, dtype={
        "IntervalStart": float, 
        "IntervalEnd": float, 
        "IntervalLength": float, 
        "SpeciesAlive": str
    })
    return df

def compute_transfers_expectancy(contemporaneity_df, output_dir):
    """Compute transfer expectancy based on contemporaneity."""
    pairs_transfers = defaultdict(float)
    
    for _, row in contemporaneity_df.iterrows():
        # Look at species alive
        species_alive = row["SpeciesAlive"].split(";")
        n_species = len(species_alive)
        
        if n_species <= 1:
            continue
            
        for i in range(n_species):
            for j in range(n_species):
                if i != j:
                    # Compute transfer expectancy
                    pairs_transfers[(species_alive[i], species_alive[j])] += row["IntervalLength"] / (n_species - 1)
    
    # Convert the dictionary to a DataFrame
    df = pd.DataFrame(
        [(k[0], k[1], v) for k, v in pairs_transfers.items()],
        columns=["Species1", "Species2", "Value"]
    )
    # Save the DataFrame to a CSV file
    output_path = os.path.join(output_dir, "transfers_expectancy.csv")
    df.to_csv(output_path, index=False)
    return df

def compute_spr(tree_path, donor, receiver, output_path, compute_spr_bin):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    command = [compute_spr_bin, tree_path, donor, receiver, output_path]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error executing SPR command for {donor}->{receiver}: {result.stderr}")
    return result.returncode == 0

def compute_spr_for_all_pairs(tree_path, output_dir, compute_spr_bin):
    """
    Compute SPR for all pairs of nodes in the tree.
    """
    from tqdm import tqdm
    
    tree = Tree(tree_path, format=1)
    tree_len = len(tree)
    
    pairs = []
    for node_1 in tree.traverse():
        if not node_1.is_root():
            node_1_ancestors = node_1.get_ancestors()
            for node_2 in tree.traverse():
                if not node_2 in node_1_ancestors and not node_2==node_1:
                    pairs.append((node_1.name, node_2.name))

    os.makedirs(output_dir, exist_ok=True)
    # Compute SPR for each pair
    for donor, receiver in tqdm(pairs, desc="Computing SPR operations"):
        output_path = os.path.join(output_dir, f"{donor}_{receiver}.newick")
        success = compute_spr(tree_path, donor, receiver, output_path, compute_spr_bin)
        
        # Check that the resulting tree has the same length as the original tree
        if success:
            try:
                output_tree = Tree(output_path, format=1)
                if len(output_tree) != tree_len:
                    print(f"Warning: Output tree for {donor}->{receiver} has different size ({len(output_tree)}) than original tree ({tree_len})")
            except Exception as e:
                print(f"Error reading output tree for {donor}->{receiver}: {str(e)}")

def are_isomorphic(tree1, tree2):
    """
    Check if two trees are isomorphic.
    """
    rf = tree1.robinson_foulds(tree2, unrooted_trees=True)[0]
    return rf == 0

def classify_sprs(spr_dir):
    # We classify the spr results according to their unrooted topology
    from tqdm import tqdm
    
    # Get list of SPR files first to show proper progress
    spr_files = [f for f in os.listdir(spr_dir) if "_" in f and f.endswith(".newick")]
    
    spr_classification = defaultdict(list)
    for spr_file in tqdm(spr_files, desc="Classifying SPR operations"):
        # match to the format "{donor}_{receiver}.newick"
        match = spr_file.split("_")
        if len(match) != 2:
            continue
        donor = match[0]
        receiver = match[1].split(".")[0]
        # Read the tree
        tree = Tree(os.path.join(spr_dir, spr_file), format=1)
        # compare to the currently stored trees. If the topology is the same, add it to the list of sprs corresponding
        # to this topology.
        # else, create a new entry in the dictionary
        found = False
        for key in spr_classification.keys():
            if are_isomorphic(tree, key):
                spr_classification[key].append((donor, receiver))
                found = True
                break
        if not found:
            spr_classification[tree] = [(donor, receiver)]
    
    # Now we need to write the classification to a table, containing the donor, the receiver, and an identifier for the
    # topology
    classification_table = []
    for i, (key, value) in enumerate(spr_classification.items()):
        topology_id = i + 1  # Use 1-based indexing for the topology identifier
        for donor, receiver in value:
            classification_table.append({
                "donor": donor,
                "receiver": receiver,
                "topology_id": topology_id
            })
    
    # Return the classification table as a DataFrame
    return pd.DataFrame(classification_table)


