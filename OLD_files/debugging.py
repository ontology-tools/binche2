from load_chebi import load_chebi

chebi_ontology = load_chebi()

# Classes that should have SMILES but weren't found
test_classes = [
    "http://purl.obolibrary.org/obo/CHEBI_59",
    "http://purl.obolibrary.org/obo/CHEBI_36694"
]

smiles_property = "https://w3id.org/chemrof/smiles_string"

for test_class in test_classes:
    print(f"\n=== Checking {test_class} ===")
    
    # Get label
    label_prop = "http://www.w3.org/2000/01/rdf-schema#label"
    label = chebi_ontology.get_annotation(test_class, label_prop)
    print(f"Label: {label}")
    
    # Get all axioms
    axioms = chebi_ontology.get_axioms_for_iri(test_class)
    print(f"Total axioms: {len(axioms)}")
    
    # Check each annotation
    found_smiles = False
    for i, axiom in enumerate(axioms):
        component = axiom.component
        if type(component).__name__ == 'AnnotationAssertion':
            ann = component.ann
            if hasattr(ann, 'ap'):
                prop = str(ann.ap)
                
                # Print all properties to see what's there
                if 'chemrof' in prop or 'smiles' in prop.lower():
                    print(f"  Axiom {i}: {prop}")
                    if hasattr(ann, 'av'):
                        print(f"    Value: {ann.av}")
                    
                    if 'smiles' in prop.lower():
                        found_smiles = True
                        print(f"    ✓ Found SMILES!")
    
    if not found_smiles:
        print("  ❌ No SMILES found in axioms")
        
        # Show ALL properties for debugging
        print("\n  All annotation properties:")
        for axiom in axioms:
            component = axiom.component
            if type(component).__name__ == 'AnnotationAssertion':
                ann = component.ann
                if hasattr(ann, 'ap'):
                    print(f"    - {ann.ap}")