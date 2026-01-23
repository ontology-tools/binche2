import os
import networkx as nx
# import matplotlib.pyplot as plt
from load_chebi import load_ontology, load_chebi
from matplotlib.patches import Ellipse
import xml.etree.ElementTree as ET
from tqdm import tqdm
from math import inf
import json
import time
from networkx.readwrite import json_graph

def id_to_name(class_id):
    id_to_name_map_file = "data/chebi_id_to_name_map.json" 

    with open(id_to_name_map_file, 'r') as f:
        id_to_name_map = json.load(f)

    prefix = "http://purl.obolibrary.org/obo/"
    if class_id.startswith(prefix):
        # remove prefix
        class_id = class_id.replace(prefix, "")
    name = id_to_name_map.get(class_id)
    return f"{name} ({class_id})" if name else class_id

def strip_prefix(class_id):
    prefix = "http://purl.obolibrary.org/obo/"
    if class_id.startswith(prefix):
        return class_id.replace(prefix, "")
    return class_id

#####################################
# Forming graph
#####################################

# def find_paths_to_root_old(ontology, start_class):
#     paths = []

#     def dfs(current_class, current_path):
#         superclasses = ontology.get_superclasses(current_class)
#         superclasses = [s for s in superclasses if s not in current_path] # Remove circular references

#         if not superclasses: # Reached root
#             paths.append(current_path)
#             return
        
#         for superclass in superclasses:
#             dfs(superclass, current_path + [superclass]) # Appends superclass to path

#     dfs(start_class, [start_class])
#     return paths

def find_paths_to_root_with_map(start_class, parents_map):
    paths = []

    def dfs(current_class, current_path):
        # if the class has no parents in the map, it's a root
        parents = parents_map.get(current_class, [])
        if not parents:
            paths.append(current_path)
            return
        
        for parent in parents:
            dfs(parent, current_path + [parent])

    dfs(start_class, [start_class])
    return paths
            
# Not used anymore
def find_paths_to_root_with_ontology(ontology, start_class, 
                       leaf_to_parents_json_file = 'data/removed_leaf_classes_to_direct_parents_map.json'): # should work for leaf classes too
    paths = []

    # Unified DFS (works for both ontology and leaf parents)
    def dfs(current_class, current_path):
        superclasses = ontology.get_superclasses(current_class)
        superclasses = [s for s in superclasses if s not in current_path]

        if not superclasses:  # reached root
            paths.append(current_path)
            return

        for superclass in superclasses:
            dfs(superclass, current_path + [superclass])


    # -----------------------------------------
    # CASE 1: class is a leaf → load parent map
    # -----------------------------------------
    with open(leaf_to_parents_json_file, "r") as f:
        leaf_to_parents = json.load(f)

    if start_class in leaf_to_parents:

        direct_parents = leaf_to_parents[start_class]
        if not isinstance(direct_parents, list):
            direct_parents = [direct_parents]

        # Start DFS from each parent, but include the leaf in the initial path
        for parent in direct_parents:
            dfs(parent, [start_class, parent])

        return paths
    
    else:
    # -----------------------------------------
    # CASE 2: class exists in filtered ontology (is not a leaf)
    # -----------------------------------------
    
        dfs(start_class, [start_class])
        return paths

def get_name(chebi_ontology, iri):
    tree = ET.parse(chebi_ontology)
    root = tree.getroot()
    ns = {
        'owl': 'http://www.w3.org/2002/07/owl#',
        'rdfs': 'http://www.w3.org/2000/01/rdf-schema#',
        'rdf': 'http://www.w3.org/1999/02/22-rdf-syntax-ns#'
    }
    
    # Find all OWL classes
    for cls in root.findall('owl:Class', ns):
        about = cls.attrib.get('{http://www.w3.org/1999/02/22-rdf-syntax-ns#}about')
        if about == iri:
            # print(f"Found class for IRI: {iri}")
            # Find rdfs:label element
            label_elem = cls.find('rdfs:label', ns)
            if label_elem is not None:
                return label_elem.text.strip()

    return None  # if not found

def create_graph_from_paths(paths):
    G = nx.DiGraph()
    for path in paths:
        for i in range(len(path) - 1):
            G.add_edge(path[i], path[i + 1])
    
    return G

# Not used anymore. If to be used, potentially remove color_map parameter
def create_graph_from_ontology(classes, classification, color_map = ['#FFB6C1', "#F44280", "#AA83A7", "#83163A", "#E63FE6", '#FFA07A', '#FF69B4'], max_n_leaf_classes=inf):
    G = nx.DiGraph()
    j = 0

    if classification == "structural":
        ontology = load_ontology("data/filtered_chebi_no_leaves_with_smiles_no_deprecated_structural.owl")
    elif classification == "functional":
        ontology = load_ontology("data/filtered_chebi_no_leaves_with_smiles_no_deprecated_functional.owl")
    else:
        ontology = load_chebi()

    for i, cls in enumerate(classes):
        print(f"Adding to graph... Starting node: {cls}")
        paths = find_paths_to_root(ontology, cls)
        H = create_graph_from_paths(paths)
        # Color the nodes of H 
        color = color_map[j % len(color_map)]
        nx.set_node_attributes(H, color, 'color')

        G = nx.compose(G, H)  # Combine graphs
        j += 1    
        if j >= max_n_leaf_classes:
            break

        print(f"Total number of starting leaf classes processed in graph: {j}")
    return G

def create_graph_from_map(classes, parent_map_json_file, max_n_leaf_classes=inf):

    with open(parent_map_json_file, "r") as f:
        parents_map = json.load(f)

    G = nx.DiGraph()
    j = 0

    for cls in classes:
        if j % 10 == 0:
            print(f"Processing class {j+1}/{len(classes)}")

        paths = find_paths_to_root_with_map(cls, parents_map)

        for path in paths:
            for i in range(len(path) - 1):
                u, v = path[i], path[i + 1]

                if not G.has_edge(u, v):
                    G.add_edge(u, v)

                    # add label once, when node first appears
                    if 'label' not in G.nodes[u]:
                        G.nodes[u]['label'] = id_to_name(u)
                    if 'label' not in G.nodes[v]:
                        G.nodes[v]['label'] = id_to_name(v)

        j += 1
        if j >= max_n_leaf_classes:
            break

    print(f"Total number of starting leaf classes processed in graph: {j}")
    return G

# Similar to the above but doesn't first create separate graphs gor each class
def create_graph_from_map_original(classes, parent_map_json_file, max_n_leaf_classes=inf):
    
    with open(parent_map_json_file, "r") as f:
        parents_map = json.load(f)

    G = nx.DiGraph()
    j = 0
    

    for cls in classes:
        if j % 10 == 0:
            print(f"Proceeesing class {j+1}/{len(classes)}")

        paths = find_paths_to_root_with_map(cls, parents_map)
        H = create_graph_from_paths(paths)

        # Add labels for Cytospace compatibility
        lable_dict = {node: id_to_name(node) for node in H.nodes()}
        nx.set_node_attributes(H, lable_dict, 'label')

        G = nx.compose(G, H)  # Combine graphs

        j += 1    
        if j >= max_n_leaf_classes:
            break

        print(f"Total number of starting leaf classes processed in graph: {j}")
    return G



# def draw_graph(G, graphing_layout, title):
#     if graphing_layout == "default":
#         pos = None  # Default layout
#     elif graphing_layout == "kamada_kawai":
#         pos = nx.kamada_kawai_layout(G)
#     elif graphing_layout == "spectral":
#         pos = nx.spectral_layout(G)
#     elif graphing_layout == "layer_based":
#         # Calculate depth of each node from root
#         roots = [n for n, d in G.in_degree() if d == 0]
#         if roots:
#             # Assign layer based on shortest path from root
#             layers = {}
#             for node in G.nodes():
#                 min_dist = float('inf')
#                 for root in roots:
#                     if nx.has_path(G, root, node):
#                         dist = nx.shortest_path_length(G, root, node)
#                         min_dist = min(min_dist, dist)
#                 layers[node] = min_dist if min_dist != float('inf') else 0
            
#             # Set subset attribute for multipartite layout
#             nx.set_node_attributes(G, layers, 'subset')
#             pos = nx.multipartite_layout(G, subset_key='subset', align='horizontal')
#     else:
#         print(f"Unknown graphing layout: {graphing_layout}. Using default.")
#         pos = None  # Default layout


#     plt.figure(figsize=(20, 10))

#     # Draw nodes with their assigned colors
#     node_colors = [G.nodes[n].get("color") for n in G.nodes()]

#     nx.draw(G, pos, with_labels=True, node_size=500, node_shape='s', font_size=8, font_weight='bold', node_color=node_colors, arrows=True,
#             arrowsize=12, edge_color='black', alpha=1)

#     plt.title(title, fontsize=12)
#     plt.show()

#####################################
# Pruning strategies
#####################################

""" RootChildrenPruner:
This pruner aims to delete the roots (the more general entities in the ontology) and its children up to a defined level.
Removes the defined level of children from the root class of the ontology."""
### 2 levels remove two levels: root and its direct children. 
### Different from Binche1 which only removes root and 2 of its direct children.

def delete_children(node, G, next_level, removed_nodes):

    # Traverse down to the specified level and remove nodes
    if next_level > 0:
        next_level -= 1
        children = list(G.predecessors(node)) # In DiGraph, predecessors are children in this case
        for child in children:
            delete_children(child, G, next_level, removed_nodes)

        # Record note removal    
        G.remove_node(node) 
        removed_nodes.add(node)

def root_children_pruner(G, levels, allow_re_execution = False, execution_count = 0):
    removed_nodes = set()
    if allow_re_execution or execution_count == 0:
        roots = [n for n, d in G.out_degree() if d == 0]
        for root in roots:
            delete_children(root, G, levels, removed_nodes)
        execution_count += 1
    return G, removed_nodes, execution_count
#TODO: look over if execution count and allow_re_execution is needed anywhere


"""Linear branch collapser pruner - remove fewer nodes:
This pruner collapses linear branches in the ontology graph by removing fewer intermediate nodes in branches 
where each node has exactly one child, effectively connecting every n:th node in such branches directly.
Keeps every n:th node in linear branches."""
# Removing fewer nodes in linear branches

# new version 
def process_branch_remove_less(head, node, G, n, removed_nodes):
    # Check if node still exists (might have been removed in another branch)
    if not G.has_node(node):
        return G
    
    branch_nodes = []
    current_node = node
    children = list(G.predecessors(current_node))
    
    while len(children) == 1:
        branch_nodes.append(current_node)
        current_node = children[0]
        # Check if current_node still exists before getting its predecessors
        if not G.has_node(current_node):
            break
        children = list(G.predecessors(current_node))
    
    # print(f"Current node: {current_node}, Children: {children}")
    last_node = current_node
    
    # Capture children BEFORE modifying the graph
    if G.has_node(last_node):
        children_before = list(G.predecessors(last_node))
    else:
        children_before = []
    
    if len(branch_nodes):
        # Determine nodes to keep
        if n == 0:
            # Remove all intermediate nodes
            nodes_to_keep = [head, last_node]
        else:
            # Keep every n:th node in the branch
            i = 0
            nodes_to_keep = [head]
            for node in branch_nodes:
                i += 1
                if i % n == 0:
                    nodes_to_keep.append(node)
            nodes_to_keep.append(last_node)
            
            # print(f"Nodes to keep in branch: {nodes_to_keep}")
            
        for node in branch_nodes:
            if node not in nodes_to_keep:
                removed_nodes.add(node)
                G.remove_node(node)
        
        for i in range(len(nodes_to_keep)-1):
            # Check that both nodes exist before adding edge
            if G.has_node(nodes_to_keep[i]) and G.has_node(nodes_to_keep[i+1]):
                if not G.has_edge(nodes_to_keep[i+1], nodes_to_keep[i]):
                    G.add_edge(nodes_to_keep[i+1], nodes_to_keep[i])
        
    # Recurse only over children that still exist in the graph
    for child in children_before:
        if G.has_node(child):  # Check before recursing
            process_branch_remove_less(last_node, child, G, n, removed_nodes)
    
    return G

def linear_branch_collapser_pruner_remove_less(G, n):
    removed_nodes = set()
    roots = [n for n, d in G.out_degree() if d == 0]
    for root in roots:
        # print(f"Processing root: {root}")
        direct_children = list(G.predecessors(root))
        # print(f"Direct children of root {root}: {direct_children}")
        for child in direct_children:
            process_branch_remove_less(root, child, G, n, removed_nodes)
    return G, removed_nodes


"""High P-Value Branch Pruner: 
Removes branches from the graph components that contain only vertices with a p-value greater than 0.05.
 * This pruner looks for branches of the ontology in the enrichment analysis result graph which
 * have only high p-values, and gets rid of them. If the branch inspected has at least one node with
 * a p-value below threshold set, then the branch is kept."""

def high_p_value_branch_pruner(G, p_value_dict, p_value_threshold = 0.05):
    removed_nodes = set()
    roots = [n for n, d in G.out_degree() if d == 0]
    
    for root in roots:

        size_before = G.number_of_nodes()
        process_node_for_p_value_pruner(root, G, p_value_dict, p_value_threshold, removed_nodes)
        size_after = G.number_of_nodes()
        # Keep pruning until nothing more is removed
        while size_after < size_before:
            size_before = size_after
            process_node_for_p_value_pruner(root, G, p_value_dict, p_value_threshold, removed_nodes)
            size_after = G.number_of_nodes()

    return G, removed_nodes

def process_node_for_p_value_pruner(node, G, p_value_dict, p_value_threshold, removed_nodes): # Returns a boolean. Tells whether node or any descendant has p-value below threshold
    
    # Check if node still exists (might have been removed in another branch)
    if not G.has_node(node):
        print(f"Node {node} no longer exists in graph.")
        return False
    
    has_good_decendant = False # whether any descendant has p-value below threshold
    children = list(G.predecessors(node))
    nodes_to_remove = []
    for child in children:
        if process_node_for_p_value_pruner(child, G, p_value_dict, p_value_threshold, removed_nodes):
            has_good_decendant = True
        else:
            nodes_to_remove.append(child)

    # Node's own p-value check
    # Get correct p-value for the node
    if "p_value_corrected" in p_value_dict.get(node, {}):
        node_p_value = p_value_dict.get(node, {}).get("p_value_corrected")
    else:
        node_p_value = p_value_dict.get(node, {}).get("p_value", None)

    if node_p_value is None:
        print(f"Warning: p-value for node {node} not found.")
        return True # Keep node if p-value not found since it is most likely a leaf node.

    if node_p_value <= p_value_threshold:
        has_good_decendant = True
    elif children == [] or (len(nodes_to_remove) == len(children)):
        # No good descendants and node itself is bad → mark for removal
        removed_nodes.add(node)
        if G.has_node(node): # Check before removing
            G.remove_node(node)
        return False # Whole branch is bad
    
    return has_good_decendant
    

def zero_degree_pruner(G):
    to_remove = []
    removed = []

    for node in G.nodes():
        if G.degree(node) == 0:
            to_remove.append(node)

    for node in to_remove:
        G.remove_node(node)
        removed.append(node)

    return G, removed







#####################################
# Converting NetworkX graph to Cytoscape compatible format
#####################################
def clean_label(label):
    """ Changes label from 'Name (CHEBI_ID)' to 'Name' """
    return label.split(' (')[0] if ' (' in label else label

def extract_chebi_id(label):
    """ Finds CHEBI_XXXX inside parentheses"""
    if ' (' in label and 'CHEBI_' in label:
        return label.split(' (')[1][:-1]  # Extract CHEBI_ID without parentheses
    return label

def graph_to_cytospace_json(G, output_file, enrichment_results=None):
    data = {"elements": []}

    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Extract the nested enrichment results
    enr_dict = enrichment_results.get('enrichment_results') if enrichment_results else None

    # Nodes
    for node, attrs in G.nodes(data=True):
        label = attrs.get("label", node)
        short_label = clean_label(label)
        node_data = {
            "id": node,
            "label": label,
            "short_label": short_label,
            "color": attrs.get("color", "#706C6C")
        }
        # Add enrichment results to nodes if provided
            
        if enr_dict and label in enr_dict: 
            enr = enr_dict[label]
            node_data["p_value"] = enr.get("p_value")
            node_data["p_value_corrected"] = enr.get("p_value_corrected")

        data["elements"].append({"data": node_data})
    
    """In enrichment_reults:
            results[class_to_check]={
            "class": id_to_name(class_to_check),
            "n_ss_annotated": n_ss_annotated,
            "n_ss_leaves": n_ss_leaves,
            "n_bg_annotated": n_bg_annotated,
            "n_bg_leaves": n_bg_leaves,
            "odds_ratio": odds,
            "p_value": p_value,
            "p_value_corrected": p_value_corrected (Added after correction)
        }"""

    # Edges
    for source, target in G.edges():
        data["elements"].append({
            "data": {
                "id": f"{source}_to_{target}",
                "source": source,
                "target": target
            }
        })

    with open(output_file, "w") as f:
        json.dump(data, f, indent=4)


"""
### Example usage
start_class = "http://purl.obolibrary.org/obo/CHEBI_33675"
start_time = time.time()
ontology = load_ontology("data/filtered_chebi_no_leaves_with_smiles_no_deprecated.owl")

paths = find_paths_to_root(ontology, start_class)
end_time = time.time()
print(f"Time taken using ontology: {end_time - start_time} seconds")

print("Using ontology:")
for path in paths:
    print(" -> ".join(path))

start_time = time.time()
paths = find_paths_to_root_with_map("data/chebi_parent_map.json", start_class)
end_time = time.time()
print(f"Time taken using parent map: {end_time - start_time} seconds")

print("Using parent map:")
for path in paths:
    print(" -> ".join(path))

"""


"""
# ---- Usage ----
# Variables
levels = 2 # Number of levels to prune from root. 1 only prunes root and it's direct neighbour, and so on.
allow_re_execution = False  # True or False. whether the pruner can be executed multiple times on a given graph.
execution_count = 0  # Counter for the number of executions

# OBS: first run cell to create G, then run the pruner function below. 
pruned_G = G.copy()
pruned_G, execution_count = root_children_pruner(pruned_G, levels, allow_re_execution, execution_count)

# graphing_layout = "kamada_kawai" # options: "default", "kamada_kawai", "spectral", "layer_based"
draw_graph(pruned_G, graphing_layout, f"Ontology graph pruned {levels} levels from root(s)")

"""