import urllib.request
import urllib.parse


SPARQL_TEMPLATE = """\
PREFIX xsd:    <http://www.w3.org/2001/XMLSchema#>
PREFIX rdfs:   <http://www.w3.org/2000/01/rdf-schema#>
PREFIX prov:   <http://www.w3.org/ns/prov#>
PREFIX wd:     <http://www.wikidata.org/entity/>
PREFIX wdt:    <http://www.wikidata.org/prop/direct/>
PREFIX p:      <http://www.wikidata.org/prop/>
PREFIX ps:     <http://www.wikidata.org/prop/statement/>
PREFIX pq:     <http://www.wikidata.org/prop/qualifier/>
PREFIX pr:     <http://www.wikidata.org/prop/reference/>
PREFIX wikibase: <http://wikiba.se/ontology#>
PREFIX schema: <http://schema.org/>


SELECT
  (xsd:integer(STRAFTER(STR(?c), "Q")) AS ?compound)
  ?compoundLabel
  ?compound_inchikey
  ?compound_smiles_conn
  ?compound_smiles_iso
  ?compound_mass
  (REPLACE(REPLACE(REPLACE(REPLACE(REPLACE(REPLACE(REPLACE(REPLACE(REPLACE(REPLACE(STR(?compound_formula_raw), "\\u2080", "0"), "\\u2081", "1"), "\\u2082", "2"), "\\u2083", "3"), "\\u2084", "4"), "\\u2085", "5"), "\\u2086", "6"), "\\u2087", "7"), "\\u2088", "8"), "\\u2089", "9") AS ?compound_formula)
  (xsd:integer(STRAFTER(STR(?t), "Q")) AS ?taxon)
  ?taxon_name
  (xsd:integer(STRAFTER(STR(?r), "Q")) AS ?ref_qid)
  ?ref
  ?ref_title
  ?ref_doi
  ?ref_date
  ?statement

WHERE {{
  {{
    SELECT
      ?c ?compound_inchikey ?compound_smiles_conn
      ?compound_smiles_iso ?compound_mass ?compound_formula_raw
      ?compoundLabel
      ?t ?taxon_name
      ?r ?ref
      ?ref_title ?ref_doi ?ref_date
      ?statement
    WHERE {{
      {{
        SELECT ?c ?compound_inchikey ?compound_smiles_conn ?t ?taxon_name ?r ?ref ?statement
        WHERE {{

  ?c wdt:P235 ?compound_inchikey;
     wdt:P233 ?compound_smiles_conn.


  ?c p:P703 ?statement.
  ?statement ps:P703 ?t;
             prov:wasDerivedFrom ?ref.
  ?ref pr:P248 ?r.
  ?t wdt:P225 ?taxon_name.

          ?t (wdt:P171*) wd:{taxon_qid}.
        }}
      }}

  OPTIONAL {{ ?r wdt:P1476 ?ref_title. }}
  OPTIONAL {{ ?r wdt:P356 ?ref_doi. }}
  OPTIONAL {{ ?r wdt:P577 ?ref_date. }}


  OPTIONAL {{ ?c wdt:P2017 ?compound_smiles_iso. }}
  OPTIONAL {{ ?c wdt:P2067 ?compound_mass. }}
  OPTIONAL {{ ?c wdt:P274 ?compound_formula_raw. }}
  OPTIONAL {{ ?c rdfs:label ?compoundLabelMul. FILTER(LANG(?compoundLabelMul) = "mul") }}
  OPTIONAL {{ ?c rdfs:label ?compoundLabelEn. FILTER(LANG(?compoundLabelEn) = "en") }}
  BIND(COALESCE(?compoundLabelMul, ?compoundLabelEn) AS ?compoundLabel)

    }}
  }}
}}"""

QLEVER_API = "https://qlever.dev/api/wikidata"


def download_lotus_csv(taxon_qid: str, output_path: str) -> None:
    """Download a LOTUS compound-taxon CSV for the given Wikidata taxon QID.

    Args:
        taxon_qid: Wikidata entity ID without prefix, e.g. "Q15978631".
        output_path: Path where the CSV file will be written.
    """
    query = SPARQL_TEMPLATE.format(taxon_qid=taxon_qid)
    filename = output_path.split("/")[-1]
    body = urllib.parse.urlencode({
        "query": query,
        "action": "csv_export",
        "filename": filename,
    }).encode("utf-8")

    req = urllib.request.Request(
        QLEVER_API,
        data=body,
        headers={"Content-Type": "application/x-www-form-urlencoded"},
        method="POST",
    )
    print(f"Downloading LOTUS data for taxon {taxon_qid} -> {output_path}")
    with urllib.request.urlopen(req) as response, open(output_path, "wb") as out:
        out.write(response.read())
    print(f"Saved {output_path}")


def download_lotus_homo_sapiens(output_path: str) -> None:
    download_lotus_csv("Q15978631", output_path)


def download_lotus_arabidopsis_thaliana(output_path: str) -> None:
    download_lotus_csv("Q158695", output_path)


if __name__ == "__main__":
    import os
    os.makedirs("data", exist_ok=True)
    download_lotus_homo_sapiens("data/lotus_homo_sapiens.csv")
    download_lotus_arabidopsis_thaliana("data/lotus_arabidopsis_thaliana.csv")
