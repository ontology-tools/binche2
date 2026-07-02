# BiNChE2

BiNChE2 is an updated version of [BiNChE](https://github.com/pcm32/BiNCheWeb/wiki/BiNChE#graph-pruning-strategies). It is a tool for ontology-based chemical enrichment analysis. It is based on the [ChEBI](https://www.ebi.ac.uk/chebi/) ontology of chemical entities.

## The Web Application
The web application is available at https://binche2.hastingslab.org/.

### Running The Analysis

The web application is hosted at https://binche2.hastingslab.org/. To run calculations locally, execute `website/app.py`. Note that all necessary data files must be generated beforehand for local execution (as described in the [Workflow](#workflow) section below).

### Study Set
On the home page, enter your study set as ChEBI IDs (one per line), as shown in the example. You can optionally provide weights for each compound (tab or space-separated). Instead of ChEBI IDs, SMILES can be used to represent a molecular entity. Each SMILES is resolved to a ChEBI ID in this order: (1) an exact string match against the local table of ChEBI leaf classes, (2) a match via the InChIKey computed from the SMILES, (3) a direct lookup through the [Chebifier](https://chebifier.hastingslab.org/) API. If none of these resolve, its predicted direct parent classes (also from Chebifier) can optionally be used for enrichment calculations instead.

Then specify the target of enrichment:

- **Structure:** Enrichment based on ChEBI structural classification (molecular composition and connectivity)
- **Role:** Enrichment based on ChEBI role classification (biological context or intended use)
- **Both:** Union of structure and role classifications (note: structural classification is significantly larger)

### Correction Method
For multiple hypothesis testing correction, p-value correction methods are available. The options are Benjamini-Hochberg, Bonferroni, and None. Benjamini-Hochberg is recommended.

### Pruning Strategies
Pruning options are available to make the graph less cluttered. The following pruners are available and found in visualitations_and_pruning.py:

* **Root Children Pruner:** Removes the roots and their children up to a defined level (number of levels being an adaptable parameter). This allows removal of more general, and less meaningful, entities in the ontology. For example, levels set to 2 will remove the roots and one level of their descendants.

* **Linear Branch Collapser Pruner:** Removes linear branches within the graph; only nodes with one parent and one child can be removed. Either a chosen number of nodes (n) in the linear branch will be kept, or all intermediate nodes in the branch can be removed (set n = 0). E.g., n = 3 will keep every third node in the branch.

* **High P-Value Branch Pruner:** Removes branches from the graph that only contain nodes with a p-value greater than 0.05 (this value can be changed). A node with a higher p-value will still be kept if it has at least one descendant with a p-value lower than the threshold. 
 
* **Zero-degree Pruner:** Removes nodes that have no connections with other nodes; that is, nodes with a total degree of zero. 

If you manually choose which pruning strategies to apply, they will be implemented once each. Alternatively, pruning strategies can be implemented in a looping manner. In this scenario, there is first a pre-loop phase where pruners are applied once, and then a loop phase where pruners are applied in a loop until no more changes are made. The looping option is:

* **Plain Enrichment Pruning Strategy:** The pre-loop phase applies the high p-value branch pruner (with a threshold of 0.05), the linear branch collapser pruner (with n = 0), and the root children pruner (levels = 2). The loop phase applies the high p-value branch pruner (with a threshold of 0.05), the linear branch collapser pruner (with n = 0), and the zero-degree vertex pruner. Benjamini-Hochberg is used as the p-value correction method.

### Calculations
For multiple hypothesis testing correction, p-value correction methods are available.
To perform enrichment analysis calculations, Fisher's exact test is used for p-value calculations. For weighted enrichment analysis however, a SaddleSum method is used. All weights must be real and positive numbers.

Have used code (translated from c to python) from https://ftp.ncbi.nlm.nih.gov/pub/qmbpmn/SaddleSum/src/, version [SaddleSum-standalone-1.2.2.tar.gz](https://ftp.ncbi.nlm.nih.gov/pub/qmbpmn/SaddleSum/src/SaddleSum-standalone-1.2.2.tar.gz) 2010-08-11 17:55  1.3M

### The Graph

The coloring based on the significance of the p-value is dependent on the values in that session; it is relative by default. Making the coloring absolute can currently only be done by changing the code (not available on the online webpage). To make this change in your local version, go to `website/templates/graph.html` and change the following line:

```const colourScaleMode = 'relative'; // 'absolute' or 'relative'```

The corrected p-value is used for the coloring if available. 

Hovering over a node displays more detailed information about it. Both raw and corrected p-values are shown, as well as its ChEBI identity number. 

Nodes can be selected by clicking on them. Right-clicking on a node provides options such as 'Select first neighbors' and 'Select descendants'. 

The graph will initially show only the most relevant branches. This means that all nodes with p-value under 0.05 will be shown, including all nodes in the paths from these nodes up to the root. If all nodes have a higher p-value, the same will be done but for nodes with p-value lower than 1. If there only exists nodes where all p-values are 1 or N/A, then all nodes will be shown.

There are options to choose the layout of the graph, which nodes are shown, and how to export it. Note that if you change the target of enrichment or pruning options under Settings, all calculations will run again. 

## Human background

Compounds from HMDB and Wikidata (LOTUS) were mapped to ChEBI leaf classes to serve as a background for enrichment. Matching to a ChEBI ID was attempted in this order: (1) a ChEBI ID already present in the source data, (2) an exact SMILES match against the local table of ChEBI leaf classes (Wikidata only), (3) an InChIKey lookup against the same local table, (4) the Chebifier API, which performs both a direct lookup and parent-class classification — for parent-class matches, the deepest class in the ChEBI hierarchy was kept to avoid overly broad annotations. Where a matched ChEBI ID corresponded to a non-leaf class in the ontology, it was expanded to its leaf descendants; classes with more than 150 leaf descendants were excluded to prevent high-level classes from disproportionately inflating the background. The resulting set of leaf classes were used to form the narrow background used in the enrichment analysis.

## Endogenous human background (Recon3D)

A second, narrower human background was built from [Recon3D](http://bigg.ucsd.edu/models/Recon3D), a genome-scale reconstruction of human metabolism, downloaded as JSON from [BiGG Models](http://bigg.ucsd.edu/). Unlike the Human background above, this one is restricted to metabolites that participate in modelled human metabolic reactions, so it excludes externally-sourced human-associated compounds (e.g. drugs, diet).

Recon3D represents each metabolite once per cellular compartment it appears in (e.g. `10fthf_c`, `10fthf_m` for the cytosolic and mitochondrial pools of the same compound), so compartment-specific entries sharing a base BiGG ID were first collapsed into a single compound record — these always carry identical formula, charge, and database cross-references, confirming they are the same chemical species. This reduced Recon3D's 5,835 metabolite entries to 2,797 unique compounds.

Each compound's listed ChEBI ID(s) were then resolved to leaf classes as follows:

1. If any listed ChEBI ID is already a leaf, **all** such leaf candidates were kept. BiGG often lists several ChEBI IDs for one compound (e.g. different protonation or tautomer states), and these are typically genuinely distinct structures rather than duplicates, so none were discarded in favour of a single "primary" one.
2. If none of the listed IDs is a leaf, each was expanded to its leaf descendants, excluding any class with more than 150 leaf descendants (the same cutoff used for the Human background, to avoid over-generic classes).
3. For compounds with no ChEBI annotation at all, a ChEBI cross-reference was attempted via [UniChem](https://www.ebi.ac.uk/unichem/), first by InChIKey, then by HMDB ID (Recon3D stores HMDB IDs in an older 5-digit format, which was zero-padded to UniChem's expected 7-digit format before lookup). Any ChEBI IDs found this way were resolved to leaves using rules 1–2 above.

Compounds for which none of the above resolved to a leaf were left out of the background. In practice these are largely abstract macromolecule/generic placeholders (e.g. `Rtotal` fatty-acyl chains, cytochromes, thioredoxin, procollagen) rather than discrete chemical structures, so they would not have been usable in a ChEBI structure-based background regardless of the matching method.

## Workflow

#### 0. UV environment
A uv environment with the following installations was used:

```uv pip install py-horned-owl rdkit networkx matplotlib pandas tqdm flask flask_sqlalchemy scipy requests```

#### 1. Load ChEBI
Download and load the ChEBI ontology by running `load_chebi.py`. In the script, the OWL file is downloaded from https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.owl and cached as `data/chebi.owl`. Re-running the script will not re-download the file if it already exists; delete `data/chebi.owl` if you want to attain a newer version. The version used in the webapplication is automatically updated on the 1st of every month.

#### 2. Remove leaf classes and save maps
To identify leaf classes and flatten the hierarchy, use `pruning_smiles.py`.

**What counts as a leaf.** A class is a leaf if it has its own valid SMILES string that is *not* a wildcard/R-group placeholder (SMILES containing the dummy atom `*`, e.g. `*C(N)C(=O)O`, are rejected via RDKit). This holds **regardless of whether the class has subclasses** — a parent class with a proper SMILES (e.g. `proline`) is a leaf in its own right, alongside its more specific children (e.g. `L-proline`, `D-proline`).

**Flattening the hierarchy (the "splice").** Because a leaf must be terminal, whenever a class sits under a leaf it is reconnected to that leaf's nearest *non-leaf* ancestors, climbing past chains of stacked leaves. So `L-proline` stops pointing at the leaf `proline` and instead points directly at the real category above it (e.g. `alpha-amino acid`), while keeping every other ancestor it already had (e.g. `D-proline` keeps `D-alpha-amino acid`). After the splice, no leaf is any class's parent, so all SMILES-bearing classes become siblings under the genuine (non-leaf) category terms. A verification step asserts the splice preserved every class's set of reachable non-leaf ancestors before any file is written.

Run task *"remove_leaves_with_smiles"* to find leaf classes, splice the hierarchy, and filter out deprecated classes. The following files are created:

- A filtered OWL file with the remaining classes (from `save_filtered_owl`). The current file in this workspace is `data/data_owl/filtered_chebi_no_leaves_with_smiles_no_deprecated.owl`.
- A **flattened** subclass map JSON file (`data/chebi_subclass_map.json`) mapping all classes to their direct subclasses after the splice (from `find_leaf_classes_with_smiles_and_deprecated`). Leaves have no subclasses here; the file also includes deprecated classes.
- A leaf-to-parents map JSON file mapping each leaf class to all of its **non-leaf** ancestors (from `find_leaf_classes_with_smiles_and_deprecated`). This file is used in later calculations.

Run task *"build_parent_map"* to create:

- `data/chebi_parent_map.json`, a map of all classes to their direct parents after the splice (deprecated classes are excluded). It is derived from the flattened subclass map, so the graph built from it shows the same sibling structure.

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

Note: because the hierarchy was flattened in step 2, a leaf inherits roles only from its **non-leaf** ancestors. It no longer inherits roles asserted directly on a leaf-ancestor — e.g. `D-proline` no longer inherits roles (such as *human metabolite*) that are asserted on the generic `proline`, which is now itself a leaf. A leaf's own direct roles are always kept.

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

The CSV contains `IRI`, `SMILES`, and `Classification`, where the classification is inferred from the class’ direct parents in the structural/functional split. Every leaf has a row here, including classes that have subclasses but carry their own valid SMILES (e.g. `proline`); classes whose SMILES is a wildcard/R-group placeholder are not leaves and do not appear. In ChEBI, classes with SMILES are expected to fall under structural roots, so entries classified as **functional** are likely misclassified and are excluded from downstream calculations (they are kept in the file for troubleshooting).


#### 5. Fisher's Calculations

First (only needed once), run task *"build_class_to_leaf_map"* in `pre_fishers_calculations.py` to create `data/class_to_leaf_descendants_map.json`, which maps each class to all of its removed leaf descendants using `data/removed_leaf_classes_to_ALL_parents_map.json`. Leaf classes are **not** keys in this map: after the splice no leaf appears as another class's ancestor, so only the genuine (non-leaf) category terms become keys. A leaf is therefore only ever counted as a member of its categories, never tested as a category itself.

Enrichment calculations can be run in `fishers_calculations.py`, but this is easiest done via the web application. Either use the website link (easiest since no preparation steps to obtain all the necessary files are needed) or run `website/app.py` locally.

#### 6. Needed for human dataset

1. Download LOTUS compound–taxon data from Wikidata via the QLever SPARQL endpoint using `wikidata/get_lotus.py`. This is run automatically by `create_files.py`, but can also be run standalone:

```bash
python wikidata/get_lotus.py
```

Output:
- `data/lotus_homo_sapiens.csv`
- `data/lotus_arabidopsis_thaliana.csv`

2. Connect the LOTUS CSVs to ChEBI IDs using `connect_lotus_csv_to_chebi_ids()` in `wikidata/get_wikidata_lotus.py`.

Output:
- `data/wikidata/created/lotus_homo_sapiens_with_chebi_ids.tsv`
- `data/wikidata/created/lotus_arabidopsis_thaliana_with_chebi_ids.tsv`

3. Extract HMDB compounds using `extract_hmdb_to_file()` in `hmdb/extract_hmdb.py`.

`data/hmdb_metabolites.xml` is required and must be downloaded manually from https://hmdb.ca/downloads (use the 'All Metabolites' XML).

Output:
- `data/hmdb_metabolites_extract.tsv`

4. Filter HMDB to only keep compounds with status "quantified" or "detected" using `filter_hmdb_statuses_main()` in `hmdb/filter_hmdb_statuses.py`.

Output: `data/hmdb_metabolites_extract_quantified_detected.tsv`

5. Find missing ChEBI IDs using `run_find_missing_chebis(source)` in `wikidata/find_missing_chebis.py` (also runnable via `jobs/run_find_missing_chebis.sh [source]`). The `source` argument must be one of the presets in `SOURCE_PRESETS`: `"lotus_hs"`, `"lotus_at"`, or `"hmdb"`.

ChEBI ID matching is attempted in this order:
- Direct ChEBI matches (LOTUS)
- Exact SMILES match against ChEBI leaf classes
- InChIKey match against ChEBI leaf classes
- Chebifier API

Output (depending on source):
- `data/wikidata/created/lotus_homo_sapiens_with_chebi_ids_updatedchebis.tsv`
- `data/wikidata/created/lotus_arabidopsis_thaliana_with_chebi_ids_updatedchebis.tsv`
- `data/hmdb_metabolites_extract_quantified_detected_updatedchebis.tsv`

6. Combine HMDB and LOTUS Homo sapiens sources using `combine_datasets()` in `wikidata/combine_human_datasets.py`. Rows with no ChEBI ID are dropped.

Output: `data/combined_hmdb_wikidata.tsv`

7. Create a file with the human leaf classes using `gather_narrow_leaves()` in `wikidata/narrow_background_fishers.py`.

Files needed:
- `compounds_tsv = "data/combined_hmdb_wikidata.tsv"`
- `leaves_csv = "data/removed_leaf_classes_with_smiles.csv"`
- `class_to_leaf_map = "data/class_to_leaf_descendants_map.json"`
- `taxon_label = "homo_sapiens"` (recorded in the output JSON for traceability)

Output: `data/human_entities_leaves.json`

#### Narrow background for a single Wikidata taxon (e.g. Arabidopsis thaliana)

This is the same workflow as above, but since there is only one source (Wikidata), steps 3, 4, and 6 (HMDB extraction/filtering and combining datasets) are skipped entirely.

1. The Arabidopsis thaliana LOTUS CSV (`data/lotus_arabidopsis_thaliana.csv`) and its ChEBI-matched TSV (`data/wikidata/created/lotus_arabidopsis_thaliana_with_chebi_ids.tsv`) are already produced in steps 1–2 above.

2. Fill in any still-missing ChEBI IDs using the `"lotus_at"` preset in `wikidata/find_missing_chebis.py`.

Output: `data/wikidata/created/lotus_arabidopsis_thaliana_with_chebi_ids_updatedchebis.tsv`

3. Build the leaf classes with `gather_narrow_leaves()` in `wikidata/narrow_background_fishers.py`, passing the file from step 2 as `compounds_tsv` and `taxon_label="arabidopsis_thaliana"`.

Output: `data/arabidopsis_thaliana_leaves.json`

#### Endogenous human background (Recon3D)

This path is independent of steps 6 and the Wikidata/HMDB workflow above; it only needs the files from steps 1–5 (`data/removed_leaf_classes_with_smiles.csv` and `data/class_to_leaf_descendants_map.json`).

Run `BiGG/get_model.py`. This:

1. Downloads the Recon3D model JSON from BiGG (`http://bigg.ucsd.edu/static/models/Recon3D.json`) to `data/Recon3D.json`.
2. Calls `gather_recon3d_leaves()`, which collapses compartment-specific metabolite entries into unique compounds and resolves each to leaf ChEBI classes (directly, via parent expansion, or via UniChem InChIKey/HMDB cross-reference, as described above). Running the script prints a breakdown of how many compounds were resolved by each method, and how many were left unresolved.

Output: `data/recon3d_leaves.json` (same `narrow_leaves` JSON shape as the other narrow backgrounds above, so it plugs into the website as `NARROW_BACKGROUND_LEAVES_JSON['endogenous_human']` without further changes).










    



