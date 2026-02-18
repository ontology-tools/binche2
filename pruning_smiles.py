from platform import node
from load_chebi import load_chebi, load_ontology
import xml.etree.ElementTree as ET
import os
import csv
import pandas as pd
import json
from collections import defaultdict

"""This script filters the ontology to remove classes that have SMILES and that are leaves"""

# def find_leaf_classes_with_smiles_and_deprecated(chebi_ontology, smiles_property, deprecated_property, subclass_map_file):

#     print("\nFinding leaf classes with SMILES...")

#     # Load or build subclass mapping for faster lookup
#     try:
#         with open(subclass_map_file, 'r') as f:
#             subclass_map = json.load(f)
#         print(f"Loaded subclass mapping from {subclass_map_file} with {len(subclass_map)} entries.")
#     except FileNotFoundError:
#         print("Building subclass mapping (may take a while)...")
#         subclass_map = defaultdict(list)
#         for idx, cls in enumerate(chebi_ontology.get_classes()):
#             subclasses = chebi_ontology.get_subclasses(cls)
#             if subclasses:
#                 subclass_map[str(cls)] = [str(sub) for sub in subclasses]
#             if (idx + 1) % 10000 == 0:
#                 print(f"Processed {idx + 1} classes for subclass mapping...")
#         with open(subclass_map_file, 'w') as f:
#             json.dump(subclass_map, f)
#         print(f"Saved subclass mapping to {subclass_map_file}.")

#     # Get all classes
#     all_classes = chebi_ontology.get_classes()
#     print(f"Total classes: {len(all_classes)}")

#     # Precompute property strings
#     smiles_prop_str = f"<{smiles_property}>"
#     deprecated_prop_str = f"<{deprecated_property}>"

#     # Find classes with SMILES
#     classes_with_smiles = []
#     leaf_classes_with_smiles = []
#     deprecated_classes = []

#     # Counters
#     i = 0
#     j = 0
#     k = 0


#     for cls in all_classes:
#         axioms = chebi_ontology.get_axioms_for_iri(cls) # Get all axioms for the class
#         cls_str = str(cls)
#         is_deprecated = False
#         has_smiles = False

#         for axiom in axioms: 
#             component = axiom.component # The component of the axiom can eg be SubClassOf, AnnotationAssertion, etc.
#             if type(component).__name__ == 'AnnotationAssertion': # Check if axiom is an annotation, e.g. a SMILES or deprecated tag
#                 ann = component.ann # Get the annotation
#                 if hasattr(ann, 'ap'): # Check if annotation has an annotation property (ap). This can eg tell us if it is a SMILES string
#                     prop_str = str(ann.ap) # Converts property IRI to string for easier comparison
                    
#                     if prop_str == deprecated_prop_str: # Check if deprecated
#                         is_deprecated = True
#                         k += 1
#                         if k % 10000 == 0:
#                             print(f"Found {k} deprecated classes so far...")
                    
#                     if prop_str == smiles_prop_str: # Check if SMILES property
#                         has_smiles = True
#                         i += 1
#                         if i % 10000 == 0:
#                             print(f"Found {i} classes with SMILES so far...")
        

#         # Add to correct lists
#         if is_deprecated:
#             deprecated_classes.append(cls_str)
#             continue  # Skip further checks for deprecated classes
#         if has_smiles:
#             classes_with_smiles.append(cls_str)
#             if cls_str not in subclass_map or len(subclass_map[cls_str]) == 0: #  # Check if leaf class (no subclasses)

#                 leaf_classes_with_smiles.append(cls_str)    
#                 j += 1
#                 if j % 10000 == 0:
#                     print(f"Found {j} leaf classes with SMILES so far...")
                
    
#     # Summary
#     print(f"\nDeprecated classes: {len(deprecated_classes)}")
#     print(f"Classes WITH SMILES: {len(classes_with_smiles)}")
#     print(f"Leaf classes with SMILES: {len(leaf_classes_with_smiles)}")

#     return set(leaf_classes_with_smiles), set(deprecated_classes), all_classes

"""
def find_leaf_classes_with_smiles_and_deprecated_old(chebi_ontology, smiles_property, deprecated_property, subclass_map_file, leaf_parents_map_file):

    print("\nFinding leaf classes with SMILES...")

    # Load or build subclass mapping for faster lookup
    try:
        with open(subclass_map_file, 'r') as f:
            subclass
            _map = json.load(f)
        print(f"Loaded subclass mapping from {subclass_map_file} with {len(subclass_map)} entries.")
    except FileNotFoundError:
        print("Building subclass mapping (may take a while)...")
        subclass_map = defaultdict(list)
        for idx, cls in enumerate(chebi_ontology.get_classes()):
            subclasses = chebi_ontology.get_subclasses(cls)
            if subclasses:
                subclass_map[str(cls)] = [str(sub) for sub in subclasses]
            if (idx + 1) % 10000 == 0:
                print(f"Processed {idx + 1} classes for subclass mapping...")
        with open(subclass_map_file, 'w') as f:
            json.dump(subclass_map, f)
        print(f"Saved subclass mapping to {subclass_map_file}.")

    # Get all classes
    all_classes = chebi_ontology.get_classes()
    print(f"Total classes: {len(all_classes)}")

    # Precompute property strings
    smiles_prop_str = f"<{smiles_property}>"
    deprecated_prop_str = f"<{deprecated_property}>"

    # Find classes with SMILES
    classes_with_smiles = []
    leaf_classes_with_smiles = []
    deprecated_classes = []

    leaf_to_parents = {}

    # Counters
    i = 0
    j = 0
    k = 0


    for cls in all_classes:
        axioms = chebi_ontology.get_axioms_for_iri(cls) # Get all axioms for the class
        cls_str = str(cls)
        is_deprecated = False
        has_smiles = False

        for axiom in axioms: 
            component = axiom.component # The component of the axiom can eg be SubClassOf, AnnotationAssertion, etc.
            if type(component).__name__ == 'AnnotationAssertion': # Check if axiom is an annotation, e.g. a SMILES or deprecated tag
                ann = component.ann # Get the annotation
                if hasattr(ann, 'ap'): # Check if annotation has an annotation property (ap). This can eg tell us if it is a SMILES string
                    prop_str = str(ann.ap) # Converts property IRI to string for easier comparison
                    
                    if prop_str == deprecated_prop_str: # Check if deprecated
                        is_deprecated = True
                        k += 1
                        if k % 1000 == 0:
                            print(f"Found {k} deprecated classes so far...")
                    
                    if prop_str == smiles_prop_str: # Check if SMILES property
                        has_smiles = True
                        i += 1
                        if i % 1000 == 0:
                            print(f"Found {i} classes with SMILES so far...")
        

        # Add to correct lists
        if is_deprecated:
            deprecated_classes.append(cls_str)
            continue  # Skip further checks for deprecated classes
        if has_smiles:
            classes_with_smiles.append(cls_str)
            if cls_str not in subclass_map or len(subclass_map[cls_str]) == 0: #  # Check if leaf class (no subclasses)

                leaf_classes_with_smiles.append(cls_str)

                # Collect parent classes for this leaf
                parents = chebi_ontology.get_superclasses(cls)
                leaf_to_parents[cls_str] = [str(p) for p in parents]   

                j += 1
                if j % 10000 == 0:
                    print(f"Found {j} leaf classes with SMILES so far...")
                
    # Save leaf to parent map
    with open(leaf_parents_map_file, 'w') as f:
        json.dump(leaf_to_parents, f, indent=2)
    print(f"Saved leaf to parents map to {leaf_parents_map_file} with {len(leaf_to_parents)} entries.")
                    
    # Summary
    print(f"\nDeprecated classes: {len(deprecated_classes)}")
    print(f"Classes WITH SMILES: {len(classes_with_smiles)}")
    print(f"Leaf classes with SMILES: {len(leaf_classes_with_smiles)}")

    return set(leaf_classes_with_smiles), set(deprecated_classes)
"""


#NEW
ancestor_cache = {}
def get_all_ancestors(cls, visited=None):
    """Recursively collect all ancestor classes, using cache to avoid recomputation."""
    cls_str = str(cls)
    if cls_str in ancestor_cache:
        return ancestor_cache[cls_str]
    
    if visited is None:
        visited = set()
    
    ancestors = []
    direct_parents = chebi_ontology.get_superclasses(cls)
    
    for parent in direct_parents:
        parent_str = str(parent)
        if parent_str not in visited:  # Avoid cycles
            visited.add(parent_str)
            ancestors.append(parent_str)
            # Recursively get ancestors of this parent
            ancestors.extend(get_all_ancestors(parent, visited))

    ancestors = list(set(ancestors))  # Remove duplicates
    ancestor_cache[cls_str] = ancestors  # Cache the result
    
    return ancestors

def find_leaf_classes_with_smiles_and_deprecated(
        chebi_ontology, 
        smiles_property,
        deprecated_property, 
        subclass_map_file, 
        leaf_parents_map_file,
        use_found_leaf_classes = False,
        removed_leaf_classes_file = None):

    print("\nFinding leaf classes with SMILES...")

    # Load or build subclass mapping for faster lookup
    try:
        with open(subclass_map_file, 'r') as f:
            subclass_map = json.load(f)
        print(f"Loaded subclass mapping from {subclass_map_file} with {len(subclass_map)} entries.")
    except FileNotFoundError:
        print("Building subclass mapping (may take a while)...")
        subclass_map = defaultdict(list)
        for idx, cls in enumerate(chebi_ontology.get_classes()):
            subclasses = chebi_ontology.get_subclasses(cls)
            if subclasses:
                subclass_map[str(cls)] = [str(sub) for sub in subclasses]
            if (idx + 1) % 10000 == 0:
                print(f"Processed {idx + 1} classes for subclass mapping...")
        with open(subclass_map_file, 'w') as f:
            json.dump(subclass_map, f)
        print(f"Saved subclass mapping to {subclass_map_file}.")

    # Get all classes
    all_classes = chebi_ontology.get_classes()
    print(f"Total classes: {len(all_classes)}")

    # Precompute property strings
    smiles_prop_str = f"<{smiles_property}>"
    deprecated_prop_str = f"<{deprecated_property}>"

    classes_with_smiles = []
    leaf_classes_with_smiles = []
    deprecated_classes = []

    if use_found_leaf_classes and removed_leaf_classes_file:
        print(f"Loading previously found leaf classes with SMILES from {removed_leaf_classes_file}...")
        df = pd.read_csv(removed_leaf_classes_file)
        leaf_classes_with_smiles = df["IRI"].tolist()
        print(f"Loaded {len(leaf_classes_with_smiles)} leaf classes with SMILES from file.")

    else: # Scan ontology to find leaf classes with SMILES

        # Counters
        i = 0
        j = 0
        k = 0

        for cls in all_classes:
            axioms = chebi_ontology.get_axioms_for_iri(cls) # Get all axioms for the class
            cls_str = str(cls)
            is_deprecated = False
            has_smiles = False

            for axiom in axioms: 
                component = axiom.component # The component of the axiom can eg be SubClassOf, AnnotationAssertion, etc.
                if type(component).__name__ == 'AnnotationAssertion': # Check if axiom is an annotation, e.g. a SMILES or deprecated tag
                    ann = component.ann # Get the annotation
                    if hasattr(ann, 'ap'): # Check if annotation has an annotation property (ap). This can eg tell us if it is a SMILES string
                        prop_str = str(ann.ap) # Converts property IRI to string for easier comparison
                        
                        if prop_str == deprecated_prop_str: # Check if deprecated
                            is_deprecated = True
                            k += 1
                            if k % 1000 == 0:
                                print(f"Found {k} deprecated classes so far...")
                        
                        if prop_str == smiles_prop_str: # Check if SMILES property
                            has_smiles = True
                            i += 1
                            if i % 1000 == 0:
                                print(f"Found {i} classes with SMILES so far...")
            

            # Add to correct lists
            if is_deprecated:
                deprecated_classes.append(cls_str)
                continue  # Skip further checks for deprecated classes
            if has_smiles:
                classes_with_smiles.append(cls_str)
                if cls_str not in subclass_map or len(subclass_map[cls_str]) == 0: # Check if leaf class (no subclasses)

                    leaf_classes_with_smiles.append(cls_str)

                    j += 1
                    if j % 10000 == 0:
                        print(f"Found {j} leaf classes with SMILES so far...")
            # Helper function to get all ancestors recursively
    
    # Build leaf to all ancestors map
    leaf_to_parents = {}
    i = 0
    for leaf in leaf_classes_with_smiles:
        all_ancestors = get_all_ancestors(leaf)
        leaf_to_parents[leaf] = all_ancestors
        i += 1
        if i % 1000 == 0:
            print(f"Processed {i} leaf classes for ancestor mapping...")
                
    # Save leaf to ALL ancestors map
    with open(leaf_parents_map_file, 'w') as f:
        json.dump(leaf_to_parents, f, indent=2)
    print(f"Saved leaf to all ancestors map to {leaf_parents_map_file} with {len(leaf_to_parents)} entries.")
                    
    # Summary
    print(f"\nDeprecated classes: {len(deprecated_classes)}")
    print(f"Classes WITH SMILES: {len(classes_with_smiles)}")
    print(f"Leaf classes with SMILES: {len(leaf_classes_with_smiles)}")

    return set(leaf_classes_with_smiles), set(deprecated_classes)

def save_filtered_owl(chebi_file, classes_with_smiles_to_remove, deprecated_classes_to_remove, output_file):

    # Check if output file already exists to avoid overwriting
    try:
        with open(output_file, 'r') as f:
            print(f"Output file {output_file} already exists. Please remove it or change name before running this function.")
            return
    except FileNotFoundError:
        pass  # File does not exist, proceed


    print(f"\nSaving filtered ontology to {output_file}...")
    
    tree = ET.parse(chebi_file) # Parsing the original OWL file to an ElementTree to be able tto remove elements etc
    root = tree.getroot()

    classes_to_remove = classes_with_smiles_to_remove.union(deprecated_classes_to_remove)
    
    # Define namespaces
    ns = {
        'owl': 'http://www.w3.org/2002/07/owl#',
        'rdf': 'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
    }
    
    # Collect elements to remove
    elements_to_remove = []
    i = 0
    j = 0

    for elem in list(root):
        # Check if element is a Class
        if elem.tag == f"{{{ns['owl']}}}Class":
            class_iri = elem.attrib.get(f"{{{ns['rdf']}}}about") # Get the class IRI
            
            # If class NOT in keep list, mark for removal
            if class_iri and class_iri in classes_to_remove:
                elements_to_remove.append(elem)
                i += 1
                if i % 1000 == 0:
                    print(f"Marked {i} classes for removal so far...")
    
    print(f"Total elements marked for removal: {len(elements_to_remove)}")

    # Remove all marked elements
    # Could combine marking elements and removing them -- check for later
    for elem in elements_to_remove:
        root.remove(elem)
        j += 1
        if j % 5000 == 0:
            print(f"Removed {j} elements so far...")
    
    # Save filtered OWL
    tree.write(output_file, encoding="utf-8", xml_declaration=True)
    print(f"✓ Done! Saved {output_file}")

def save_leaf_classes_with_smiles(leaf_classes, chebi_ontology, smiles_property, output_file, structural_classes, functional_classes):
    """Save leaf classes with SMILES to a CSV file."""

    print(f"\nSaving {len(leaf_classes)} leaf classes with SMILES to {output_file}...")

    seen = set()
    rows = []

    for cls in leaf_classes:
        if cls in seen:
            print(f"Skipping duplicate class: {cls}")
            continue  # Skip duplicates
        seen.add(cls)

        axioms = chebi_ontology.get_axioms_for_iri(cls)
        smiles = None

        for axiom in axioms:
            component = axiom.component
            if type(component).__name__ == 'AnnotationAssertion':
                ann = component.ann
                if hasattr(ann, 'ap'):
                    prop_str = str(ann.ap)
                    if prop_str == f"<{smiles_property}>":
                        smiles = str(ann.av)  # Get the SMILES value
                        break  # No need to check more axioms for this class

        # Classficiation from parents
        classification = "neither"
        parents = set()
        try: 
            parents = chebi_ontology.get_superclasses(cls)
            if any(p in structural_classes for p in parents):
                classification = "structural"
            elif any(p in functional_classes for p in parents):
                classification = "functional"
        except Exception as e:
            print(f"⚠️ Could not get superclasses for {cls}: {e}")
            classification = "unknown"
            pass


        rows.append([cls, smiles, classification])

    # Write to CSV
    with open(output_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["IRI", "SMILES", "Classification"])
        writer.writerows(rows)

    length = len(rows)
    print(f"✓ Done! Saved {length} leaf classes with SMILES to {output_file}")

# Only used when deprecated classes were removed in a separate step. Should not be needed now. # Used in build_parent_map
def find_deprecated_classes(ontology, deprecated_property):
    """Find all deprecated classes in the ontology."""
    print("Finding deprecated classes...")
    all_classes = ontology.get_classes()
    deprecated_classes = set()
    i = 0

    for cls in all_classes:
        axioms = ontology.get_axioms_for_iri(cls)
        for axiom in axioms:
            component = axiom.component
            if type(component).__name__ == 'AnnotationAssertion':
                ann = component.ann
                if hasattr(ann, 'ap'):
                    prop_str = str(ann.ap)
                    if prop_str == f"<{deprecated_property}>":
                        deprecated_classes.add(cls)
                        i += 1
                        if i % 1000 == 0:
                            print(f"Found {i} deprecated classes so far...")
                        break  # No need to check more axioms for this class
    print(f"Total deprecated classes found: {len(deprecated_classes)}")
    return deprecated_classes

# Only used when deprecated classes were removed in a separate step. Should not be needed now.
def remove_classes_from_owl(input_file, classes_to_remove, output_file):
    """Remove given classes from OWL file and save a new version."""
    print(f"\nRemoving {len(classes_to_remove)} classes from {input_file}...")

    tree = ET.parse(input_file)
    root = tree.getroot()

    ns = {
        'owl': 'http://www.w3.org/2002/07/owl#',
        'rdf': 'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
    }

    removed = 0
    for elem in list(root):
        if elem.tag == f"{{{ns['owl']}}}Class":
            iri = elem.attrib.get(f"{{{ns['rdf']}}}about")
            if iri and iri in classes_to_remove:
                root.remove(elem)
                removed += 1
                if removed % 5000 == 0:
                    print(f"Removed {removed} deprecated classes so far...")

    tree.write(output_file, encoding="utf-8", xml_declaration=True)
    print(f"✓ Done! Removed {removed} classes.")
    print(f"Saved to {output_file}")

"""For building parent map""" #maybe move to separate file later

def build_parent_map(ontology, output_json, deprecated_property):
    """Build parent map excluding deprecated classes"""
    
    # Find deprecated classes first
    raw_deprecated = find_deprecated_classes(ontology, deprecated_property)
    deprecated_classes = {str(c) for c in raw_deprecated}

    print(f"Excluding {len(deprecated_classes)} deprecated classes")
    
    # Convert deprecated classes to strings for consistent comparison
    
    all_classes = ontology.get_classes()
    string_ids = {cls: str(cls) for cls in all_classes}  # fast lookup table
    parent_map = {}
    
    i = 0
    for cls in all_classes:
        if i % 1000 == 0:
            print(f"Processed {i} classes for parent map...")
        i += 1
        
        cls_id = string_ids[cls]
        
        # Skip if deprecated (now comparing strings)
        if cls_id in deprecated_classes:
            continue
            
        # Get parents, excluding deprecated ones
        direct_parents = ontology.get_superclasses(cls)
        parents = [string_ids[parent] 
                   for parent in direct_parents 
                   if string_ids[parent] not in deprecated_classes]
        
        parent_map[cls_id] = parents
    
    with open(output_json, 'w') as f:
        json.dump(parent_map, f, indent=2)

    # Print statistics
    root_classes = [cls for cls, parents in parent_map.items() if not parents]
    print(f"Saved {len(parent_map)} non-deprecated classes")
    print(f"  - {len(root_classes)} root classes (no parents)")
    print(f"  - {len(parent_map) - len(root_classes)} classes with parents")
    print(f"Saved parent map to {output_json}")

def get_name(root, iri, ns):
    """Return the rdfs:label for a class IRI."""
    for cls in root.findall('owl:Class', ns):
        about = cls.attrib.get('{http://www.w3.org/1999/02/22-rdf-syntax-ns#}about')
        if about == iri:
            label_elem = cls.find('rdfs:label', ns)
            if label_elem is not None:
                return label_elem.text.strip()
    return None

def shorten_parent_map(input_json, output_json):
    """ Shorten the IRIs in the parent map by removing the common prefix. """
    prefix = "http://purl.obolibrary.org/obo/"

    with open(input_json, "r") as f:
        parent_map = json.load(f)

    short_map = {}

    for iri, parents in parent_map.items():
        # Convert key
        short_key = iri.replace(prefix, "")

        # Convert parent IRIs
        short_parents = [p.replace(prefix, "") for p in parents]

        short_map[short_key] = short_parents

    with open(output_json, "w") as f:
        json.dump(short_map, f, indent=2)

    print(f"Saved shortened parent map → {output_json}")

def map_names_to_classes(chebi_ontology, output_json):
    ns = {
        'owl': 'http://www.w3.org/2002/07/owl#',
        'rdfs': 'http://www.w3.org/2000/01/rdf-schema#',
        'rdf': 'http://www.w3.org/1999/02/22-rdf-syntax-ns#'
    }

    tree = ET.parse(chebi_ontology)
    root = tree.getroot()

    prefix = "http://purl.obolibrary.org/obo/"

    iri_to_name = {}
    i = 0

    for cls in root.findall('owl:Class', ns):
        iri = cls.attrib.get('{http://www.w3.org/1999/02/22-rdf-syntax-ns#}about')
        if iri is None:
            continue

        # keep only the last part, e.g. CHEBI_12345
        if iri.startswith(prefix):
            short_id = iri[len(prefix):]
        else:
            short_id = iri  # fallback

        label_elem = cls.find('rdfs:label', ns)
        name = label_elem.text.strip() if label_elem is not None else None

        iri_to_name[short_id] = name

        i += 1
        if i % 1000 == 0:
            print(f"Processed {i} classes for name mapping...")

    with open(output_json, "w") as f:
        json.dump(iri_to_name, f, indent=2)

    print(f"Saved short-ID → name mapping to {output_json}")

    
if __name__ == "__main__":

    task = "build_parent_map" # Options: "remove_leaves_with_smiles", "save_removed_leaf_classes", "build_parent_map", "map_names_to_classes"
    print(f"Selected task: {task}")
    # property IRI
    smiles_property = "https://w3id.org/chemrof/smiles_string"
    deprecated_property = "http://www.w3.org/2002/07/owl#deprecated"

    if task == "remove_leaves_with_smiles": # Remove leaf classes with SMILES and deprecated classes from OWL 
                                            # (a new filtered OWL file will be saved)
        filtered_output_file = "data/filtered_chebi_no_leaves_with_smiles_no_deprecated_newww.owl"

        chebi_file = "data/chebi.owl"
        subclass_map_file = "data/chebi_subclass_map.json" 
        leaf_parents_map_file = "data/removed_leaf_classes_to_ALL_parents_map.json"

        use_found_leaf_classes = True # Set to True to use previously found leaf classes with SMILES from CSV file
        removed_leaf_classes_file = "data/removed_leaf_classes_with_smiles.csv" # Only needed if use_found_leaf_classes is True

        if os.path.exists(filtered_output_file):
            print(f"Output file {filtered_output_file} already exists. Are you sure you want to overwrite it? If so, please remove it before running this script.")
        else:
            print("Running code to remove leaf classes with SMILES and deprecated classes...")
            chebi_ontology = load_chebi()
            classes_with_smiles, deprecated_classes = find_leaf_classes_with_smiles_and_deprecated(
                chebi_ontology,
                smiles_property, 
                deprecated_property, 
                subclass_map_file, 
                leaf_parents_map_file, 
                use_found_leaf_classes, 
                removed_leaf_classes_file)
    
            # Comment out the next two lines if you do not want to save the filtered OWL in a new file
            save_filtered_owl(chebi_file, classes_with_smiles, deprecated_classes, filtered_output_file)
            filtered_ontology = load_ontology(filtered_output_file) # just to confirm it loads

    elif task == "save_removed_leaf_classes":
        print("Running code to save removed leaf classes with SMILES...")

        output_file = "data/removed_leaf_classes_with_smiles_new.csv"
        subclass_map_file = "data/chebi_subclass_map.json"

        structural_ontology = load_ontology("data/filtered_chebi_no_leaves_with_smiles_no_deprecated_structural.owl")
        functional_ontology = load_ontology("data/filtered_chebi_no_leaves_with_smiles_no_deprecated_functional.owl")

        # Extract sets of class IRIs for quick membership checking
        structural_classes = set(structural_ontology.get_classes())
        functional_classes = set(functional_ontology.get_classes())

        print(f"Structural classes loaded: {len(structural_classes)}")
        print(f"Functional classes loaded: {len(functional_classes)}")

        chebi_ontology = load_chebi()
        classes_with_smiles, _, _ = find_leaf_classes_with_smiles_and_deprecated(chebi_ontology, smiles_property, deprecated_property, subclass_map_file)

        # Save CSV with classification
        save_leaf_classes_with_smiles(
            classes_with_smiles,
            chebi_ontology,
            smiles_property,
            output_file,
            structural_classes,
            functional_classes,

        )  

        df = pd.read_csv(output_file)
        print(len(df), "rows total")
        print(df["IRI"].nunique(), "unique IRIs")
        print(df[df["SMILES"].isna()])

    elif task == "build_parent_map":
        print("Running code to build parent map excluding deprecated classes...")
        ontology = load_chebi()
        output_json = "data/chebi_parent_map.json" 
        # shortened_output_json = "data/chebi_parent_map_shortened_id.json"
        if os.path.exists(output_json):
            print(f"Output file {output_json} already exists. Are you sure you want to overwrite it? If so, please remove it before running this script.")
        else:
            build_parent_map(ontology, output_json, deprecated_property) # Depracated classes are removed inside the function

        # shortened_output_json = "data/chebi_parent_map_shortened_id.json"
        # if os.path.exists(shortened_output_json):
        #     print(f"Shortened output file {shortened_output_json} already exists. Are you sure you want to overwrite it? If so, please remove it before running this script.")
        # else:
        #     shorten_parent_map(output_json, shortened_output_json)
        #     print(f"Saved shortened parent map to {shortened_output_json}")
        
    
    elif task == "map_names_to_classes":
        print("Running code to map names to classes...")
        chebi_ontology = "data/chebi.owl"
        output_json = "data/chebi_id_to_name_map.json"
        if os.path.exists(output_json):
            print(f"Output file {output_json} already exists. Are you sure you want to overwrite it? If so, please remove it before running this script.")
        else:
            map_names_to_classes(chebi_ontology, output_json)

    else: 
        print("No valid task selected. Please choose a valid task.")