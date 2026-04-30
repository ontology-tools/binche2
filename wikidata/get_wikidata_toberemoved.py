import requests
import pandas as pd
import time
import os
from pathlib import Path

# --- Configuration ---
ENDPOINT = "https://query.wikidata.org/sparql"
USER_AGENT = "MetabolomicsEnrichment/1.0 (miranda.carlsson@idiap.ch)"

QUERY = """
SELECT DISTINCT ?compound ?compoundLabel ?chebi ?smiles_canonical ?taxon_name WHERE {
  ?compound wdt:P31/wdt:P279* wd:Q407595. 
  
  OPTIONAL { ?compound wdt:P683 ?chebi . }
  OPTIONAL { ?compound wdt:P233 ?smiles_canonical . }
  OPTIONAL { 
    { ?compound wdt:P703 ?taxon . }
    UNION
    { ?compound wdt:P1582 ?taxon . }
    ?taxon wdt:P225 ?taxon_name .
  }
  SERVICE wikibase:label { bd:serviceParam wikibase:language "en" . }
}
LIMIT 1
"""

def run_sparql_query(query, endpoint=ENDPOINT, user_agent=USER_AGENT):
    headers = {
        "User-Agent": user_agent,
        "Accept": "application/sparql-results+json"
    }
    response = requests.get(
        endpoint,
        params={"query": query, "format": "json"},
        headers=headers,
        timeout=60  # seconds before your script gives up
    )
    response.raise_for_status()  # raises an error for 4xx/5xx responses
    return response.json()

def results_to_dataframe(results):
    rows = []
    for binding in results["results"]["bindings"]:
        row = {key: binding[key]["value"] if key in binding else None
               for key in ["compound", "compoundLabel", "chebi", 
                           "smiles_canonical", "taxon_name"]}
        rows.append(row)
    return pd.DataFrame(rows)

# --- Run ---
try:
    print("Querying Wikidata...")
    raw = run_sparql_query(QUERY)
    df = results_to_dataframe(raw)
    print(f"Got {len(df)} results")
    print(df.head())

    output_dir = Path(__file__).resolve().parent / "data" / "wikidata"
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_filename = "wikidata_compounds.csv"
    # If the file already exists, add a number suffix to avoid overwriting
    if os.path.exists(output_dir / csv_filename):
        base, ext = os.path.splitext(csv_filename)
        i = 1
        while os.path.exists(output_dir / f"{base}_{i}{ext}"):
            i += 1
        csv_filename = f"{base}_{i}{ext}"

    output_path = output_dir / csv_filename
    df.to_csv(output_path, index=False)
    print(f"Saved CSV to: {output_path}")
except requests.exceptions.Timeout:
    print("Query timed out — try narrowing with a filter or LIMIT")
except requests.exceptions.HTTPError as e:
    print(f"HTTP error: {e.response.status_code} — {e.response.text}")

