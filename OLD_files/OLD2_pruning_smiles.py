from load_chebi import load_chebi, load_ontology
import xml.etree.ElementTree as ET
import os
import csv
import pandas as pd

"""This script filters the ontology to remove classes that have SMILES and that are leaves"""

def find_leaf_classes_with_smiles(chebi_ontology, smiles_property, deprecated_property):
    print("\nFinding leaf classes with SMILES...")

    # Get all classes
    all_classes = chebi_ontology.get_classes()
    print(f"Total classes: {len(all_classes)}")

    # Find classes with SMILES
    classes_with_smiles = []
    leaf_classes_with_smiles = []
    classes_without_smiles = []
    deprecated_classes = []

    # Counters
    i = 0
    j = 0
    k = 0
    l = 0

    for cls in all_classes:
        axioms = chebi_ontology.get_axioms_for_iri(cls) # Get all axioms for the class

        is_deprecated = False 
        has_smiles = False

        for axiom in axioms: 
            component = axiom.component # The component of the axiom can eg be SubClassOf, AnnotationAssertion, etc.
            if type(component).__name__ == 'AnnotationAssertion': # Check if axiom is an annotation, e.g. a SMILES or deprecated tag
                ann = component.ann # Get the annotation
                if hasattr(ann, 'ap'): # Check if annotation has an annotation property (ap). This can eg tell us if it is a SMILES string
                    prop_str = str(ann.ap) # Converts property IRI to string for easier comparison
                    
                    if prop_str == f"<{deprecated_property}>": # Check if deprecated
                        is_deprecated = True
                    
                    if prop_str == f"<{smiles_property}>": # Check if SMILES property
                        has_smiles = True

        # Put class in appropriate list
        if is_deprecated:
            deprecated_classes.append(cls) 
            k += 1
            if k % 10000 == 0:
                print(f"Found {k} deprecated classes so far...")

        elif has_smiles:
            classes_with_smiles.append(cls)
            i += 1
            if i % 10000 == 0:
                print(f"Found {i} classes with SMILES so far...")
            
            # Check if leaf class (no subclasses)

            # This is extremely slow -- consider optimising later
            if not chebi_ontology.get_subclasses(cls):
                leaf_classes_with_smiles.append(cls)
                l += 1
                if l % 1000 == 0:
                    print(f"Found {l} leaf classes with SMILES so far...")
                

        else: # Not necessary, including for counting
            classes_without_smiles.append(cls)
            j += 1
            if j % 10000 == 0:
                print(f"Found {j} classes without SMILES so far...")

    # Summary
    print(f"\nDeprecated classes: {len(deprecated_classes)}")
    print(f"Active classes WITH SMILES: {len(classes_with_smiles)}")
    print(f"Active classes WITHOUT SMILES: {len(classes_without_smiles)}")
    print(f"Active leaf classes with SMILES: {len(leaf_classes_with_smiles)}")

    # # Show examples of classes from each category
    # print("\n--- Examples with SMILES ---")
    # for cls in list(classes_with_smiles)[:3]:
    #     print(f"{cls}")

    # print("\n--- Examples without SMILES ---")
    # for cls in list(classes_without_smiles)[:3]:
    #     print(f"{cls}")

    # print("\n--- Examples deprecated ---")
    # for cls in list(deprecated_classes)[:3]:
    #     print(f"{cls}")

    return set(leaf_classes_with_smiles), all_classes

def save_filtered_owl(chebi_file, classes_to_remove, output_file):

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


if __name__ == "__main__":

    task = "save_removed_leaf_classes" # Options: "remove_leaves_with_smiles", "remove_all_deprecated", "save_removed_leaf_classes"
    print(f"Selected task: {task}")
    # property IRI
    smiles_property = "https://w3id.org/chemrof/smiles_string"
    deprecated_property = "http://www.w3.org/2002/07/owl#deprecated"

    if task == "remove_leaves_with_smiles":
        filtered_output_file = "data/filtered_chebi_no_leaves_with_smiles.owl"

        chebi_file = "data/chebi.owl"

        if os.path.exists(filtered_output_file):
            print(f"Output file {filtered_output_file} already exists. Are you sure you want to overwrite it? If so, please remove it before running this script.")
        else:
            print("Running code")
            chebi_ontology = load_chebi()
            classes_with_smiles, all_classes = find_leaf_classes_with_smiles(chebi_ontology, smiles_property, deprecated_property)
            save_filtered_owl(chebi_file, classes_with_smiles, filtered_output_file)
            filtered_ontology = load_ontology(filtered_output_file) # just to confirm it loads

    elif task == "remove_all_deprecated":
        print("Running code to remove all deprecated classes...")
        
        input_file = "data/filtered_chebi_no_leaves_with_smiles.owl"
        output_file = "data/filtered_chebi_no_leaves_with_smiles_no_deprecated.owl"

        chebi_ontology = load_ontology(input_file)
        deprecated_classes = find_deprecated_classes(chebi_ontology, deprecated_property)
        remove_classes_from_owl(input_file, deprecated_classes, output_file)

    elif task == "save_removed_leaf_classes":
        print("Running code to save removed leaf classes with SMILES...")

        output_file = "data/removed_leaf_classes_with_smiles.csv"

        structural_ontology = load_ontology("data/filtered_chebi_no_leaves_with_smiles_no_deprecated_structural.owl")
        functional_ontology = load_ontology("data/filtered_chebi_no_leaves_with_smiles_no_deprecated_functional.owl")

        # Extract sets of class IRIs for quick membership checking
        structural_classes = set(structural_ontology.get_classes())
        functional_classes = set(functional_ontology.get_classes())

        print(f"Structural classes loaded: {len(structural_classes)}")
        print(f"Functional classes loaded: {len(functional_classes)}")

        chebi_ontology = load_chebi()
        classes_with_smiles, _ = find_leaf_classes_with_smiles(chebi_ontology, smiles_property, deprecated_property)

        # Save CSV with classification
        save_leaf_classes_with_smiles(
            classes_with_smiles,
            chebi_ontology,
            smiles_property,
            output_file,
            structural_classes,
            functional_classes
        )  

        df = pd.read_csv("data/removed_leaf_classes_with_smiles.csv")
        print(len(df), "rows total")
        print(df["IRI"].nunique(), "unique IRIs")
        print(df[df["SMILES"].isna()])

    else: 
        print("No valid task selected. Please choose a valid task.")