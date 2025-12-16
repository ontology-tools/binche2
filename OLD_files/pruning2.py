# pruning.py
from load_chebi import load_chebi
import time

chebi_ontology = load_chebi()

# ---- Helpers ----

def curie_to_iri(curie: str) -> str:
    """Convert CURIE 'CHEBI:12345' -> full IRI used by the ontology."""
    prefix, num = curie.split(":")
    return f"http://purl.obolibrary.org/obo/{prefix}_{num}"

def get_descendants(ontology, root_curie: str):
    """Return a set of descendant IRIs (including the root IRI)."""
    root_iri = curie_to_iri(root_curie)
    descendants = set()
    to_visit = [root_iri]

    while to_visit:
        current = to_visit.pop()
        if current in descendants:
            continue
        descendants.add(current)
        # get_subclasses expects a full IRI
        try:
            subs = list(ontology.get_subclasses(current))
        except Exception:
            subs = []
        to_visit.extend(subs)
    return descendants

def safe_get_annotation(ontology, iri: str, pred_iri: str):
    """Return first annotation value for predicate IRI or None."""
    try:
        vals = ontology.get_annotation(iri, pred_iri)
    except Exception:
        vals = None
    if not vals:
        return None
    # ontology.get_annotation often returns a list-like; take first
    return vals[0] if isinstance(vals, (list, tuple)) else vals

# ---- Build remove set ----

print("Computing descendants (this may take a little while)...")
t0 = time.time()
role_desc = get_descendants(chebi_ontology, "CHEBI:50906")      # role
particle_desc = get_descendants(chebi_ontology, "CHEBI:36342")  # subatomic particle
ids_to_remove = role_desc.union(particle_desc)
t1 = time.time()
print(f"Total IDs to remove: {len(ids_to_remove)}  (computed in {t1-t0:.1f}s)")

# ---- Write cleaned OBO ----

OUTPATH = "chebi_clean.obo"
label_pred = "http://www.w3.org/2000/01/rdf-schema#label"
smiles_pred = "http://purl.obolibrary.org/obo/chebi#smiles"

all_classes = list(chebi_ontology.get_classes())
total = len(all_classes)
print(f"Total classes in ontology: {total}")

written = 0
processed = 0
progress_step = max(1, total // 20)  # ~5% steps

t0 = time.time()
with open(OUTPATH, "w") as out_file:
    for cls in all_classes:
        processed += 1
        # progress
        if processed % progress_step == 0 or processed <= 5:
            elapsed = time.time() - t0
            print(f"Processed {processed}/{total} classes (elapsed {elapsed:.1f}s)")

        if cls in ids_to_remove:
            continue

        # get canonical label and SMILES (if any) via safe API
        label = safe_get_annotation(chebi_ontology, cls, label_pred)
        smiles = safe_get_annotation(chebi_ontology, cls, smiles_pred)

        # if no label, skip (mirrors Perl which writes only terms with an id/name)
        if not label:
            continue

        term_lines = []
        term_lines.append("[Term]\n")
        term_lines.append(f"id: {cls}\n")
        term_lines.append(f"name: {label}\n")
        if smiles:
            term_lines.append(f"property_value: SMILES \"{smiles}\"\n")

        out_file.writelines(term_lines)
        written += 1

t1 = time.time()
print(f"Finished. Wrote {written} terms to {OUTPATH} in {t1-t0:.1f}s.")
