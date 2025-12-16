# BiNChE2 - Current workflow

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

#### 3. Split up the ontology based on structure 
To split the ontology classes based on if they are functional or structural, use the file `pruning_split_up_structure.py`. Run task *"split_structural_functional"*. 

- Two versions of the previously filtered ontology are created; one with functional and one with structral classes.

#### 4. Save a file with the removed leaf classes
Go back to file `pruning_smiles.py` and run task *"save_removed_leaf_classes"* to save the removed leaf classes in a csv file. 

- CSV file with all removed leaf classes with their SMILES and if they are functional or structural.

#### 5. Fishers calculations

Use file `fishers_calculations.py`. 

First (only needed once), run task *"build_class_to_leaf_map"* to create a json file with all classes mapped to all of their removed leaf classes (currently, unless updated, the file does not include the leaf classes as own entries).

Run the task *"count_total_removed_leaves"* to calculate the total number of removed leaf classes for either structural or functional classes. The csv with removed leaf classes is used for this.

Run the task *"count_removed_classes_for_class"* to calculate the removed leaf classes for a specific class. The new json file is used for this.

Run the task "enrichment_analysis_plain"
 ...

### Using the webpage

#### Get it working

Run file `website/app.py`. 
To be bale to proceed with the next steps, all necessary files are needed but these should already be obtained by previous steps (or saved so on a cloud somewhere......)

#### Enrichment analysis

#### The graph

The colouring based on the significance of the p-value is dependant on the values in that actual session. So it is simply relative. The corrected p-value is used if available.

Nodes can be selected by clicking on them.

GOAL (as in not yet implemented).
Right-clicking on a node gives the user options such as 'Select first neighbors' and 'Select descendants'

Find other goals here (tools implemented in BiNChE):
https://github.com/pcm32/BiNCheWeb/wiki/BiNChE#graph-pruning-strategies

TODO: setg 1 att implementera så att högerklickning faktiskt får upp något. Would be lovely

...







    



