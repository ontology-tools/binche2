# BiNChE2 (preliminary name)

BiNChE2 is an updated version of [BiNChE](https://github.com/pcm32/BiNCheWeb/wiki/BiNChE#graph-pruning-strategies). It is a tool for ontology-based chemical enrichment analysis. It is based on the [ChEBI](https://www.ebi.ac.uk/chebi/) ontology of chemical entities.

## BiNChE2 - Current workflow

#### 0. UV environment
A UV environment with the following installations was used:

```uv pip install py_horned_owl rdkit networkx matplotlib pandas tqdm flask flask_sqlalchemy scipy```

#### 1. Load ChEBI
Download and load ChEBI using `load_chebi.py` (TODO: double check how download was computed)

#### 2. Remove leaf classes and save maps
To remove all leaf classes with SMILES from the ontology, use file `pruning_smiles.py`. 

Run task *"remove_leaves_with_smiles"* to find and filter out depracated classes and leaf classes with SMILES. The following files are created:

- A filtered owl file with the remaining classes (from function *save_filtered_owl*)
- A subclass map json file containing all classes with their direct subclasses. Used in the function. Observe that the filealso contains dercated classes (from function *find_leaf_classes_with_smiles_and_deprecated*)
- Leaf to parents map json file mapping the removed leaf classes to all of their respective ancestor. This file will be used later (from function *find_leaf_classes_with_smiles_and_deprecated*)

Run task *"build_parent_map"* to attain:

- A map (*chebi_parent_map.json*) with all classes pointing to their respective direct parents. To be used when creating paths for the graph creations.
- The same corresponding map (*chebi_parent_map_shortened_id.json*) but not with the long iri:s, just eg "CHEBI_123"

Run task *"map_names_to_classes*" to build:

- A map (*"chebi_id_to_name_map.json*") that maps the Chebi IDs (the short IDs, eg Chebi_111) to their actual names.

#### 2.5 (New) Save maps connected to the roles of the classes

Maps that include the roles of the classes are needed for some the enrichment calculations. These are made in the script `prepare_role_calculations.py`.

First, run the task "find has_role connections". This will produce a map from all classes to its direct roles (not including any roles that its ancestors have). This file is named 

* "data/class_to_direct_roles_map.json", 

Second, run the task "build leaf to all roles map". This will produce the files 

* *"data/removed_leaf_classes_to_ALL_roles_map.json"* 
* *"data/roles_to_leaves_map.json"* 

The first file maps the leaf classes to all of its direct roles, all of its ancestors direct roles, and the ancestors to all of these roles. The second file does the opposite; it maps all of the role classes to all leaf classes that are connected to that role.

Third, run the task "build class to all roles map" to attain the file

* *"data/class_to_all_roles_map.json"*

Each class is mapped to its direct roles, its ancestors roles, as well as the decessors to all of those roles.


#### 3. Split up the ontology based on structure 
To split the ontology classes based on if they are functional or structural, use the file `pruning_split_up_structure.py`. Run task *"split_structural_functional"*. 

- Two versions of the previously filtered ontology are created; one with functional and one with structral classes.

#### 4. Save a file with the removed leaf classes
Go back to file `pruning_smiles.py` and run task *"save_removed_leaf_classes"* to save the removed leaf classes in a csv file. 

- CSV file with all removed leaf classes with their SMILES.
- Observe that their classification is mentioned as either functional or structural. There should however only exist classes with SMILES (that can stand for themselves) under the structural role. We are thus only only interested in these. The others will not be taken into consideration in calculations.


#### 5. Fishers calculations

Use file `fishers_calculations.py`. 

First (only needed once), run task *"build_class_to_leaf_map"* to create a json file with all classes mapped to all of their removed leaf classes (currently, unless updated, the file does not include the leaf classes as own entries).

Run the task *"count_total_removed_leaves"* to calculate the total number of removed leaf classes for either structural or functional classes. The csv with removed leaf classes is used for this.

Run the task *"count_removed_classes_for_class"* to calculate the removed leaf classes for a specific class. The new json file is used for this.

Run the task "enrichment_analysis_plain"
 ...

## Behind the webpage

### Get it working

Run file `website/app.py`. 
To be able to proceed with the next steps, all necessary files are needed but these should already be obtained by previous steps (or saved on a cloud somewhere......)

### Enrichment analysis
For the correction of multiple hypotehsis testing, P-value correction methods are available. The options are Benjamini-Hochberg, Bonferroni, and None, whereas the first one is recommended.

### Correction method
For the purpose 

### Pruning strategies
Pruning options are available for the possibility to make the graph less clustered. The following pruners can be chosen:

* **Root Children Pruner:** Removes the roots and its children up to a defined level (number of levels being an adaptable paramter). This allows removal of more general, and less meaningful, entities in the ontology. For example levels set to 2 will remove the roots and one level of their descendants.

* **Linear Branch Collapser Pruner:** Removes linear branches within the graph; only nodes with one parent and one child can be removed. Either all intermediate nodes in the branch can be removed (set n = 0) or a a chosen number of nodes (n) in the linear branch will be kept. E.g. n = 3 will keep every third node in the branch.

* **High P-Value Branch Pruner:** Removes branches from the graph components that only contain nodes with a p-value greater than 0.05 (this value can be changed). Therefore a node with a higher p-value can still be kept, if it has at least one decendant with a p-value that lower than the threshold. 
 
* **Zero-degree pruner:** Removes nodes that has no connections with other nodes; that is nodes with a total degree of zero. 

There is either the option to choose which pruning options that will be applied. The chosen srategies will then be implemented ones each. Alternatively, different pruning options can be implemented in a looping manner. There is first a pre-loop phase where pruners are applied ones, and thereafter a loop-phase where pruners are applied until in a loop until no more changes are made.

* **Plain Enrichment pruning strategy**: The pre-loop phase applies the high p-value branch pruner (with a threshold of 0.5), the linear branch collapser pruner (with n = 0), and the root children pruner (levels = 2). The loop-phase applies the high p-value branch pruner (with a threshold of 0.5), the linear branch collapser pruner (with n = 0), and the zero-degree vertex pruner. Benjamini-Hochberg is used at the P-value corrcetion method.


### The graph

The colouring based on the significance of the p-value is dependant on the values in that actual session. So it is simply relative. To make the the colouring absolute, go to `website/templates/graph.html`, and change the following line:


```const colourScaleMode = 'relative'; // 'absolute' or 'relative'```

 The corrected p-value is used if available. 

Hovering over a node gives more displays more detailed information about it. Bothe raw and correcte dp-values are shown, as well as its ChEBI identity number. 

Nodes can be selected by clicking on them. Right-clicking on a node provides options such as 'Select first neighbors' and 'Select descendants'. 












    



