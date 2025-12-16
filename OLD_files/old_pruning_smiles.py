from load_chebi import load_chebi, load_ontology
import xml.etree.ElementTree as ET
import os

""" This script saves a file with only classes that have SMILES strings (and their ancestors) from the ChEBI ontology. """

def find_classes_with_smiles(chebi_ontology, smiles_property, deprecated_property):

    # Get all classes
    all_classes = chebi_ontology.get_classes()
    print(f"Total classes: {len(all_classes)}")

    # Find classes with SMILES
    classes_with_smiles = []
    classes_without_smiles = []
    deprecated_classes = []

    # Counters
    i = 0
    j = 0
    k = 0

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
        # Can remove the elif and else but interesting to see counts separately
        elif has_smiles:
            classes_with_smiles.append(cls)
            i += 1
            if i % 10000 == 0:
                print(f"Found {i} classes with SMILES so far...")
        else: 
            classes_without_smiles.append(cls)
            j += 1
            if j % 10000 == 0:
                print(f"Found {j} classes without SMILES so far...")

    # Summary
    print(f"\nDeprecated classes: {len(deprecated_classes)}")
    print(f"Active classes WITH SMILES: {len(classes_with_smiles)}")
    print(f"Active classes WITHOUT SMILES: {len(classes_without_smiles)}")
    print(f"\nCoverage (excluding deprecated): {len(classes_with_smiles) / (len(classes_with_smiles) + len(classes_without_smiles)) * 100:.1f}%")

    # Show examples of classes from each category
    print("\n--- Examples with SMILES ---")
    for cls in list(classes_with_smiles)[:3]:
        print(f"{cls}")

    print("\n--- Examples without SMILES ---")
    for cls in list(classes_without_smiles)[:3]:
        print(f"{cls}")

    print("\n--- Examples deprecated ---")
    for cls in list(deprecated_classes)[:3]:
        print(f"{cls}")

    return set(classes_with_smiles), all_classes

# def get_parents(ontology, class_iri):
#     """Return immediate parents (superclasses) of a given class IRI."""
#     parents = set()
#     for ax in ontology.get_axioms_for_iri(class_iri):
#         comp = ax.component
#         if comp.__class__.__name__ == "SubClassOf":
#             parents.add(str(comp.sup))
#     return parents

def get_parents(ontology, class_iri):
    """Return immediate named superclasses (exclude restrictions)."""
    parents = set()
    for ax in ontology.get_axioms_for_iri(class_iri):
        comp = ax.component
        if comp.__class__.__name__ == "SubClassOf":
            sup = str(comp.sup)
            # Keep only named class IRIs, ignore ObjectSomeValuesFrom. This is so that we only keep ancestors and not other restrictions
            if "ObjectSomeValuesFrom" not in sup:
                parents.add(sup)
    return parents

def get_all_ancestors(ontology, class_iri):
    """Return all transitive superclasses (ancestors) of a class."""
    to_visit = [class_iri] # Start with the given class
    visited = set() # Keep track of visited classes

    while to_visit: # While there are classes left to visit
        current = to_visit.pop() # Get a class to visit. Will be last element
        if current in visited:
            continue # Already visited
        visited.add(current) # Mark current as visited
        parents = get_parents(ontology, current) # Get immediate parents
        for p in parents:
            if p not in visited:
                to_visit.append(p) # Add unvisited parents to visit list
    
    visited.remove(class_iri)  # Remove the class itself if present as we only want ancestors
    return visited





def save_filtered_owl(chebi_ontology,classes_with_smiles, output_file):

    # Check if output file already exists to avoid overwriting
    try:
        with open(output_file, 'r') as f:
            print(f"Output file {output_file} already exists. Please remove it or change name before running this function.")
            return
    except FileNotFoundError:
        pass  # File does not exist, proceed

    print(f"\nSaving filtered ontology to {output_file}...")
    
    tree = ET.parse('data/chebi.owl') # Parsing the original OWL file to an ElementTree to be able tto remove elements etc
    root = tree.getroot()
    
    # Define namespaces
    ns = {
        'owl': 'http://www.w3.org/2002/07/owl#',
        'rdf': 'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
    }
    
    classes_to_keep = classes_with_smiles
    for cls in set(classes_with_smiles):
        classes_to_keep.update(get_all_ancestors(chebi_ontology, cls)) # Also keep all ancestors of classes with SMILES
    
    # Collect elements to remove
    elements_to_remove = []
    i = 0
    j = 0
    k = 0
    for elem in list(root):
        # Check if element is a Class
        if elem.tag == f"{{{ns['owl']}}}Class":
            class_iri = elem.attrib.get(f"{{{ns['rdf']}}}about") # Get the class IRI
            
            # If class NOT in keep list, mark for removal
            if class_iri and class_iri not in classes_to_keep:
                elements_to_remove.append(elem)
                i += 1
                if i % 1000 == 0:
                    print(f"Marked {i} classes for removal so far...")
        

        # # Also remove Axioms that reference classes not in keep list
        # elif elem.tag == f"{{{ns['owl']}}}Axiom":
        #     # Check annotatedSource
        #     annotated_source = elem.find(f"{{{ns['owl']}}}annotatedSource")
        #     if annotated_source is not None:
        #         source_iri = annotated_source.attrib.get(f"{{{ns['rdf']}}}resource")
        #         if source_iri and source_iri not in classes_to_keep:
        #             elements_to_remove.append(elem)
        #             j += 1
        #             if j % 1000 == 0:
        #                 print(f"Marked {j} axioms for removal so far...")

        # Removing axioms gives a very small dataset, do not think we want to do this(?)

    
    print(f"Total elements marked for removal: {len(elements_to_remove)}")

    # Remove all marked elements
    for elem in elements_to_remove:
        root.remove(elem)
        j += 1
        if j % 1000 == 0:
            print(f"Removed {j} elements so far...")
    
    # Save filtered OWL
    tree.write(output_file, encoding="utf-8", xml_declaration=True)
    print(f"✓ Done! Saved {output_file}")


if __name__ == "__main__":

    filtered_output_file = "data/filtered_chebi_neeew.owl"

    # SMILES property IRI
    smiles_property = "https://w3id.org/chemrof/smiles_string"
    deprecated_property = "http://www.w3.org/2002/07/owl#deprecated"

    if os.path.exists(filtered_output_file):
        print(f"Output file {filtered_output_file} already exists. Are you sure you want to overwrite it? If so, please remove it before running this script.")
    else:
        print("Running code")
        chebi_ontology = load_chebi()
        classes_with_smiles, all_classes = find_classes_with_smiles(chebi_ontology, smiles_property, deprecated_property)
        save_filtered_owl(chebi_ontology, classes_with_smiles, filtered_output_file)
        filtered_ontology = load_ontology(filtered_output_file) # just to confirm it loads

    