import requests
import pandas as pd

ENDPOINT = "https://query.wikidata.org/sparql"

QUERY = """
SELECT ?compound ?compoundLabel ?inchikey ?smiles ?chebi_id WHERE {
  ?compound wdt:P703 wd:Q15978631 .
  OPTIONAL { ?compound wdt:P235 ?inchikey . }
  OPTIONAL { ?compound wdt:P233 ?smiles . }
  OPTIONAL { ?compound wdt:P683 ?chebi_id . }
  SERVICE wikibase:label { bd:serviceParam wikibase:language "en" . }
}
"""

r = requests.get(
    ENDPOINT,
    params={"query": QUERY, "format": "json"},
    headers={"User-Agent": "Miranda/ChEBI-N research (EPFL)"}
)
r.raise_for_status()

bindings = r.json()["results"]["bindings"]
rows = []
for b in bindings:
    rows.append({
        "wikidata_id": b["compound"]["value"].split("/")[-1],
        "name":        b.get("compoundLabel", {}).get("value", ""),
        "inchikey":    b.get("inchikey", {}).get("value", ""),
        "smiles":      b.get("smiles", {}).get("value", ""),
        "chebi_id":    b.get("chebi_id", {}).get("value", ""),
    })

df = pd.DataFrame(rows)
print(f"{len(df)} rows retrieved")

# Save to file
df.to_csv("data/new_lotus_homo_sapiens.csv", index=False)
print("Saved to data/new_lotus_homo_sapiens.csv")