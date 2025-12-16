from load_chebi import load_chebi

chebi_ontology = load_chebi()

label_pred = "http://www.w3.org/2000/01/rdf-schema#label"
smiles_pred = "http://purl.obolibrary.org/obo/chebi#smiles"

print("Testing on first classes...\n")

all_classes = list(chebi_ontology.get_classes())
count = 0
for cls in all_classes[1000:2000]:
    label = chebi_ontology.get_annotation(cls, label_pred)
    smiles = chebi_ontology.get_annotation(cls, smiles_pred)

    if label:
        label_text = label[0] if isinstance(label, list) else label
        #print(f"Class: {cls}")
        #print(f"  Label: {label_text}")
        if smiles:
            smiles_text = smiles[0] if isinstance(smiles, list) else smiles
            print(f"  SMILES: {smiles_text}")
            print("SMILES found!")
        print()
        count += 1

print(f"✅ Printed {count} labeled classes.")

#no smiles found in first 1000 classes