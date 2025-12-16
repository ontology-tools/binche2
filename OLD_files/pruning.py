from load_chebi import load_chebi

# Load the ChEBI ontology
chebi_ontology = load_chebi()

# Convert CURIE like "CHEBI:50906" to full IRI
def curie_to_iri(curie):
    prefix, num = curie.split(":")
    return f"http://purl.obolibrary.org/obo/{prefix}_{num}"

# Get all descendants of a root class (including root)
def get_descendants(ontology, root_curie):
    root_iri = curie_to_iri(root_curie)
    descendants = set()
    to_visit = [root_iri]

    while to_visit:
        current = to_visit.pop()
        # print(f"Visiting: {current}")
        if current not in descendants:
            descendants.add(current)
            subclasses = list(ontology.get_subclasses(current))
            to_visit.extend(subclasses)
    return descendants

# Generate sets of IDs to remove
role_ids = get_descendants(chebi_ontology, "CHEBI:50906")
particle_ids = get_descendants(chebi_ontology, "CHEBI:36342")
ids_to_remove = role_ids.union(particle_ids)
print(f"Total IDs to remove: {len(ids_to_remove)}")

# Write cleaned ontology to a new .obo file
with open("chebi_clean??.obo", "w") as out_file:
    for cls in chebi_ontology.get_classes():
        print(f"Processing class: {cls}")
        # Skip classes in remove list
        if cls in ids_to_remove:
            continue

        axioms = chebi_ontology.get_axioms_for_iri(cls)
        term_lines = ["[Term]\n"]  # start a new term


###### works until here
        for ax in axioms:
            pred_iri = ax.annotation_property.iri #'pyhornedowl.model.AnnotatedComponent' object has no attribute 'annotation_property'
            for obj in ax.objects:
                if pred_iri == "http://www.w3.org/2000/01/rdf-schema#label":
                    term_lines.append(f"id: {cls}\n")
                    term_lines.append(f"name: {obj}\n")
                elif pred_iri == "http://purl.obolibrary.org/obo/chebi#smiles":
                    term_lines.append(f"property_value: SMILES \"{obj}\"\n")
                else:
                    # Optional: include other annotations
                    pass

        if any(line.startswith("id:") for line in term_lines):
            out_file.writelines(term_lines)

print("chebi_clean.obo created successfully.")
