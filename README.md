# BiNChE2

BiNChE2 is an updated version of [BiNChE](https://github.com/pcm32/BiNCheWeb/wiki/BiNChE#graph-pruning-strategies). It is a tool for ontology-based chemical enrichment analysis. It is based on the [ChEBI](https://www.ebi.ac.uk/chebi/) ontology of chemical entities.

## The Web Application
The web application is available at https://binche2.hastingslab.org/

### Running The Analysis

The web application is hosted at https://binche2.hastingslab.org/. To run calculations locally, execute `website/app.py`. Note that all necessary data files must be generated beforehand for local execution (as described in the [Workflow](#workflow) section below).

### Study Set
On the home page, enter your study set as ChEBI IDs (one per line), as shown in the example. You can optionally provide weights for each compound (tab or space-separated). Then specify the target of enrichment:

- **Structure:** Enrichment based on ChEBI structural classification (molecular composition and connectivity)
- **Role:** Enrichment based on ChEBI role classification (biological context or intended use)
- **Both:** Union of structure and role classifications (note: structural classification is significantly larger)

### Correction Method
For multiple hypothesis testing correction, p-value correction methods are available. The options are Benjamini-Hochberg, Bonferroni, and None. Benjamini-Hochberg is recommended.

### Pruning Strategies
Pruning options are available to make the graph less cluttered. The following pruners are available:

* **Root Children Pruner:** Removes the roots and their children up to a defined level (number of levels being an adaptable parameter). This allows removal of more general, and less meaningful, entities in the ontology. For example, levels set to 2 will remove the roots and one level of their descendants.

* **Linear Branch Collapser Pruner:** Removes linear branches within the graph; only nodes with one parent and one child can be removed. Either a chosen number of nodes (n) in the linear branch will be kept, or all intermediate nodes in the branch can be removed (set n = 0). E.g., n = 3 will keep every third node in the branch.

* **High P-Value Branch Pruner:** Removes branches from the graph that only contain nodes with a p-value greater than 0.05 (this value can be changed). A node with a higher p-value will still be kept if it has at least one descendant with a p-value lower than the threshold. 
 
* **Zero-degree Pruner:** Removes nodes that have no connections with other nodes; that is, nodes with a total degree of zero. 

If you manually choose which pruning strategies to apply, they will be implemented once each. Alternatively, pruning strategies can be implemented in a looping manner. In this scenario, there is first a pre-loop phase where pruners are applied once, and then a loop phase where pruners are applied in a loop until no more changes are made. The looping option is:

* **Plain Enrichment Pruning Strategy:** The pre-loop phase applies the high p-value branch pruner (with a threshold of 0.05), the linear branch collapser pruner (with n = 0), and the root children pruner (levels = 2). The loop phase applies the high p-value branch pruner (with a threshold of 0.05), the linear branch collapser pruner (with n = 0), and the zero-degree vertex pruner. Benjamini-Hochberg is used as the p-value correction method.

### The Graph

The coloring based on the significance of the p-value is dependent on the values in that session; it is relative by default. Making the coloring absolute can currently only be done by changing the code (not available on the online webpage). To make this change in your local version, go to `website/templates/graph.html` and change the following line:

```const colourScaleMode = 'relative'; // 'absolute' or 'relative'```

The corrected p-value is used for the coloring if available. 

Hovering over a node displays more detailed information about it. Both raw and corrected p-values are shown, as well as its ChEBI identity number. 

Nodes can be selected by clicking on them. Right-clicking on a node provides options such as 'Select first neighbors' and 'Select descendants'. 

The graph will initially show only the most relevant branches. This means that all nodes with p-value under 0.05 will be shown, including all nodes in the paths from these nodes up to the root. If all nodes have a higher p-value, the same will be done but for nodes with p-value lower than 1.

There are options to choose the layout of the graph, which nodes are shown, and how to export it. Note that if you change the target of enrichment or pruning options under Settings, all calculations will run again. 

## Workflow

#### 0. UV environment
A uv environment with the following installations was used:

```uv pip install py-horned-owl rdkit networkx matplotlib pandas tqdm flask flask_sqlalchemy scipy```

#### 1. Load ChEBI
Download and load the ChEBI ontology by running `load_chebi.py`. In the script, the OWL file is downloaded from https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.owl and cached as `data/chebi.owl`. Re-running the script will not re-download the file if it already exists; delete `data/chebi.owl` if you want to attain a newer version. The version currently used in the webapplication was updated 2025.10.09.

#### 2. Remove leaf classes and save maps
To remove leaf classes with SMILES from the ontology, use `pruning_smiles.py`.

Run task *"remove_leaves_with_smiles"* to find and filter out deprecated classes and leaf classes with SMILES. The following files are created:

- A filtered OWL file with the remaining classes (from `save_filtered_owl`). The current file in this workspace is `data/data_owl/filtered_chebi_no_leaves_with_smiles_no_deprecated.owl`.
- A subclass map JSON file containing all classes with their direct subclasses (from `find_leaf_classes_with_smiles_and_deprecated`). The file also includes deprecated classes.
- A leaf-to-parents map JSON file mapping the removed leaf classes to all of their ancestors (from `find_leaf_classes_with_smiles_and_deprecated`). This file is used in later calculations.

Run task *"build_parent_map"* to create:

- `data/chebi_parent_map.json`, a map of all classes to their direct parents (deprecated classes are excluded). This file is later used when creating paths for graph construction.

Run task *"map_names_to_classes"* to build:

- `data/chebi_id_to_name_map.json`, which maps short CHEBI IDs (e.g., `CHEBI_111`) to their names.

#### 2.5 Save maps connected to the roles of the classes

Maps that include the roles of the classes are needed for some enrichment calculations. These are made in `prepare_role_calculations.py`.

First, run the task *"find has_role connections"*. This parses the OWL file directly and produces a map from all classes to their **direct** roles (not including any roles that ancestors have):

- `data/class_to_direct_roles_map.json`

This task also calls `create_leaves_to_all_roles_map`, which writes 

- `data/removed_leaf_classes_to_ALL_roles_map.json` (using `data/removed_leaf_classes_to_ALL_parents_map.json` and `data/chebi_parent_map.json`)

Second, run the task *"build leaf to all roles map"*. This builds:

- `data/removed_leaf_classes_to_ALL_roles_map.json`
- `data/roles_to_leaves_map.json`

The first file maps each removed leaf class to (a) its direct roles, (b) roles inherited from ancestor classes, and (c) **ancestors of those roles** in the role hierarchy. The second file is the inverse: it maps each role class to all leaf classes connected to that role.

Third, run the task *"build class to all roles map"* to create:

- `data/class_to_all_roles_map.json`

Here, each class is mapped to its direct roles, roles inherited from ancestor classes, and **descendants** of those roles in the role hierarchy.


#### 3. Split up the ontology based on structure
To split the ontology classes into functional vs structural, use `pruning_split_up_structure.py` and run task *"split_structural_functional"*.

- This creates three versions of the previously filtered ontology: `_structural.owl`, `_functional.owl`, and `_unknown.owl` (the last contains classes that are not classified by the two roots, and has just been used for trouble shooting). The current workspace includes:
	- `data/data_owl/filtered_chebi_no_leaves_with_smiles_no_deprecated_structural.owl`
	- `data/data_owl/filtered_chebi_no_leaves_with_smiles_no_deprecated_functional.owl`

#### 4. Save a file with the removed leaf classes
Go back to `pruning_smiles.py` and run task *"save_removed_leaf_classes"* to save the removed leaf classes in a CSV file. The current file in this workspace is 

- `data/removed_leaf_classes_with_smiles.csv`

The CSV contains `IRI`, `SMILES`, and `Classification`, where the classification is inferred from the class’ direct parents in the structural/functional split. In ChEBI, classes with SMILES are expected to fall under structural roots, so entries classified as **functional** are likely misclassified and are excluded from downstream calculations (they are kept in the file for troubleshooting).


#### 5. Fisher's Calculations

Use `fishers_calculations.py` for enrichment runs.

First (only needed once), run task *"build_class_to_leaf_map"* in `pre_fishers_calculations.py` to create `data/class_to_leaf_descendants_map.json`, which maps each class to all of its removed leaf descendants using `data/removed_leaf_classes_to_ALL_parents_map.json`. By default, leaf classes themselves are **not** included as keys in this map.

Enrichment calculations can be run in `fishers_calculations.py`, but this is easiest done via the web application. Either use the website link (easiest since no preparation steps to obtain all the necessary files are needed) or run `website/app.py` locally.













    



