import requests
import time

BASE = "http://bigg.ucsd.edu/api/v2"

def get_metabolite_bigg_ids(model_id):
    """Fetch all metabolite BiGG IDs for a model, stripping compartment suffixes."""
    url = f"{BASE}/models/{model_id}/metabolites"
    resp = requests.get(url)
    resp.raise_for_status()
    mets = resp.json()["results"]
    # Each has a bigg_id (e.g. "pyr") and compartment_bigg_id (e.g. "c")
    # Strip compartment to get the universal metabolite ID
    return set(m["bigg_id"] for m in mets)

# Pick your models — recommended combination:
models = {
    "E. coli": "iML1515",     
    "Yeast":   "iMM904",       
    "Human":   "Recon3D",      
}

metabolites = {}
for organism, model_id in models.items():
    print(f"Fetching {model_id} ({organism})...")
    metabolites[organism] = get_metabolite_bigg_ids(model_id)
    print(f"  → {len(metabolites[organism])} unique metabolites")
    time.sleep(0.2)

# Compute pairwise and triple intersections
ec  = metabolites["E. coli"]
ye  = metabolites["Yeast"]
hu  = metabolites["Human"]

print(f"\nE. coli ∩ Yeast:         {len(ec & ye)}")
print(f"E. coli ∩ Human:         {len(ec & hu)}")
print(f"Yeast   ∩ Human:         {len(ye & hu)}")
print(f"\nE. coli ∩ Yeast ∩ Human: {len(ec & ye & hu)}  ← your core metabolism set")

core = ec & ye & hu
print(f"\nCore metabolites (BiGG IDs):")
for mid in sorted(core):
    print(f"  {mid}")