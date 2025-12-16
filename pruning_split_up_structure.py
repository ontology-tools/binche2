from load_chebi import load_chebi, load_ontology
import xml.etree.ElementTree as ET
import os

## Root IRIs for functional and structural hierarchies

# Functional
ROLE_ROOT_IRI = "http://purl.obolibrary.org/obo/CHEBI_50906" # role 

# Structural
MOLECULAR_ROOT_IRI = "http://purl.obolibrary.org/obo/CHEBI_23367" # molecular entity  
CHEMICAL_ROOT_IRI = "http://purl.obolibrary.org/obo/CHEBI_59999" # chemical substance
"""
The root for structural is "chemical entity" which has four subclasses: 
'atom', 'chemical substance', 'group' and 'molecular entity'. 
Want to include both 'chemical substance' and 'molecular entity' but remove 'atom' and 'group'.
"""

def get_descendants(ontology, root_iri):
    """Recursively collect all descendants (subclasses) of a given class."""
    descendants = set()
    to_visit = [root_iri] # Start with the root class
    while to_visit:
        current = to_visit.pop()
        try:
            subclasses = ontology.get_subclasses(current) # Returns all direct subclasses/children
        except Exception:
            continue
        for sub in subclasses:
            if sub not in descendants:
                descendants.add(sub)
                to_visit.append(sub)
    return descendants


def identify_structural_vs_functional(chebi_ontology):
    """Classify classes as structural or functional based on ChEBI hierarchy."""
    all_classes = set(chebi_ontology.get_classes())
    print(f"Total classes in ontology: {len(all_classes)}")

    # Identify functional hierarchies
    print("Collecting descendants of 'role' (functional)...")
    functional_classes = get_descendants(chebi_ontology, ROLE_ROOT_IRI)
    # Add the root itself
    functional_classes.add(ROLE_ROOT_IRI)
    print(f"Functional classes: {len(functional_classes)}")

    # Identify structural hierarchies
    print("Collecting descendants of 'molecular entity' (structural)...")
    structural_classes_molecular = get_descendants(chebi_ontology, MOLECULAR_ROOT_IRI)
    # Add the root itself
    structural_classes_molecular.add(MOLECULAR_ROOT_IRI)
    print(f"Structural classes with molecular entity root: {len(structural_classes_molecular)}")
    print("Collecting descendants of 'chemical substance' (structural)...")
    structural_classes_chemical = get_descendants(chebi_ontology, CHEMICAL_ROOT_IRI)
    # Add the root itself
    structural_classes_chemical.add(CHEMICAL_ROOT_IRI)
    print(f"Structural classes with chemical substance root: {len(structural_classes_chemical)}")
    # Combine both structural sets
    structural_classes = structural_classes_molecular.union(structural_classes_chemical)
    print(f"Total structural classes: {len(structural_classes)}")

    # Identify any classes not belonging to either category
    known = functional_classes.union(structural_classes)
    unknown_classes = all_classes - known

    print(f"Unknown/unclassified: {len(unknown_classes)}")

    # Check if overlaps exist
    overlap = functional_classes.intersection(structural_classes)
    if overlap:
        print(f"!! Warning: Overlapping classes found between structural and functional: {len(overlap)}")

    return structural_classes, functional_classes, unknown_classes


def split_owl_by_type(structural_classes, functional_classes, unknown_classes, input_file):
    """Write two new OWL files: one for structural and one for functional classes."""
    ns = {
        'owl': 'http://www.w3.org/2002/07/owl#',
        'rdf': 'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
    }

    print(f"Parsing {input_file} ...")
    tree = ET.parse(input_file)
    root = tree.getroot()

    # Function to create and save filtered tree
    def create_filtered_tree(classes_to_keep, output_file):
        new_root = ET.Element(root.tag, root.attrib)
        removed = 0

        for elem in list(root):
            if elem.tag == f"{{{ns['owl']}}}Class":
                iri = elem.attrib.get(f"{{{ns['rdf']}}}about")
                if iri and iri in classes_to_keep:
                    new_root.append(elem)
                else:
                    removed += 1
            else:
                new_root.append(elem)

        ET.ElementTree(new_root).write(output_file, encoding="utf-8", xml_declaration=True)
        print(f"✓ Saved {output_file} ({removed} elements removed)")

    # Save all splits
    
    base = os.path.basename(input_file)
    name, ext = os.path.splitext(base)
    create_filtered_tree(structural_classes, f"data/{name}_structural{ext}")
    create_filtered_tree(functional_classes, f"data/{name}_functional{ext}")
    create_filtered_tree(unknown_classes, f"data/{name}_unknown{ext}")

def check_roots_of_unknown(unknown_classes, chebi_ontology):
    """For each unknown class, find its top-most ancestors (roots)
    and give detailed stats for those under 'chemical entity'."""
    print(f"\nChecking root ancestors of {len(unknown_classes)} unknown classes...")

    CHEMICAL_ENTITY_IRI = "http://purl.obolibrary.org/obo/CHEBI_24431"
    root_summary = {}
    chemical_entity_sub_summary = {}
    i = 0

    for cls in unknown_classes:
        i += 1
        if i % 1000 == 0:
            print(f"Processed {i} classes...")

        try:
            ancestors = chebi_ontology.get_ancestors(cls)
        except Exception as e:
            print(f"⚠️ Could not get ancestors for {cls}: {e}")
            continue

        ancestors = set(a for a in ancestors if a != cls)

        if not ancestors:
            root_summary.setdefault("NO_ANCESTOR", []).append(cls)
            continue

        # Find top-most ancestors
        root_ancestors = []
        for anc in ancestors:
            supers = chebi_ontology.get_superclasses(anc)
            if not supers or all(s not in ancestors for s in supers):
                root_ancestors.append(anc)

        # Store results per root
        for root in root_ancestors:
            root_summary.setdefault(root, []).append(cls)

            # Special handling if root is 'chemical entity'
            if root == CHEMICAL_ENTITY_IRI:
                try:
                    subs = chebi_ontology.get_subclasses(root)
                    relevant_subs = [s for s in subs if s in ancestors]
                    if relevant_subs:
                        for sub in relevant_subs:
                            chemical_entity_sub_summary.setdefault(sub, []).append(cls)
                except Exception as e:
                    print(f"⚠️ Could not get subclasses of chemical entity: {e}")

    # --- PRIMARY SUMMARY ---
    print("\n===============================")
    print(" ROOT CLASS SUMMARY FOR UNKNOWNS")
    print("===============================")
    for root, members in sorted(root_summary.items(), key=lambda x: len(x[1]), reverse=True):
        label = chebi_ontology.get_label(root) if hasattr(chebi_ontology, "get_label") else root
        print(f"• {label} ({root}): {len(members)} unknowns under it")

    if "NO_ANCESTOR" in root_summary:
        print("\nExamples of isolated classes (no ancestors):")
        for example in root_summary["NO_ANCESTOR"][:5]:
            print(f" - {example}")

    # --- SECONDARY SUMMARY (under chemical entity) ---
    if chemical_entity_sub_summary:
        print("\n----------------------------------")
        print(" DETAIL: Classes under 'chemical entity'")
        print("----------------------------------")
        for subroot, members in sorted(chemical_entity_sub_summary.items(),
                                       key=lambda x: len(x[1]), reverse=True):
            label = chebi_ontology.get_label(subroot) if hasattr(chebi_ontology, "get_label") else subroot
            print(f"  • {label} ({subroot}): {len(members)} unknowns")

    return root_summary, chemical_entity_sub_summary





if __name__ == "__main__":

    task = "check_descendants" # Options: "split_structural_functional", "check_unknown", "check_descendants"

    if task == "split_structural_functional":

        print("Loading filtered ChEBI ontology...")
        # file = "data/filtered_chebi_no_leaves_with_smiles.owl"
        file = "data/filtered_chebi_no_leaves_with_smiles_no_deprecated.owl"
        print("Using file:", file, " Check if correct")
        chebi_ontology = load_ontology(file) 
        # chebi_ontology = load_chebi()

        print("\nIdentifying structural vs functional classes...")
        structural_classes, functional_classes, unknown_classes = identify_structural_vs_functional(chebi_ontology)

        print("\nSplitting ontology into structural and functional OWL files...")
        split_owl_by_type(structural_classes, functional_classes, unknown_classes, file)

        print("\n✓ Done!")

    elif task == "check_unknown": # Should not be needed anymore
        print("Checking unknown/unclassified classes...")

        # unknown_classes_file = "data/filtered_chebi_no_leaves_with_smiles_unknown.owl"
        unknown_classes_file = "data/filtered_chebi_no_leaves_with_smiles_no_deprecated_unknown.owl"
        unknown_classes = load_ontology(unknown_classes_file)

        file = "data/filtered_chebi_no_leaves_with_smiles.owl"
        chebi_ontology = load_ontology(file) 

        root_summary = check_roots_of_unknown(unknown_classes.get_classes(), chebi_ontology)

    elif task == "check_descendants": # Should not be needed anymore
        print("Checking descendants of a specific class...")
        class_iri = "http://purl.obolibrary.org/obo/CHEBI_52217"
        
        chebi_ontology = load_chebi()

        descendants = get_descendants(chebi_ontology, class_iri)
        print(f"Class {class_iri} has {len(descendants)} descendants in full ontology.")
        # print descendants
        for d in descendants:
            print(f" - {d}")

    else:
        print("Unknown task. Please choose a valid task.")
