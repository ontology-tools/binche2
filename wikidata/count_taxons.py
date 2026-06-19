# Count which 5 taxons are the most common in Wikidata

import pandas as pd

taxon_file = 'data/wikidata/compound_reference_taxon.tsv'
taxon_names_file = 'data/wikidata/taxa.tsv'


def top_n_taxons(n=5):
    pairs = pd.read_csv(taxon_file, sep="\t")
    taxa = pd.read_csv(taxon_names_file, sep="\t")

    taxon_to_name = dict(zip(taxa["wikidataId"], taxa["names_pipe_separated"]))

    top_counts = pairs["taxon"].value_counts().head(n)

    print(f"Top {n} most common taxons:")
    for taxon_id, count in top_counts.items():
        name = taxon_to_name.get(taxon_id, "Unknown")
        print(f"{name}\t{taxon_id}\t{count}")

    return top_counts

# Count the number of classes in class_file.tsv
def count_classes_in_file(class_file):
    df = pd.read_csv(class_file, sep="\t")
    num_classes = len(df)
    print(f"Number of classes in {class_file}: {num_classes}")
    return num_classes


if __name__ == "__main__":
    top_n_taxons(5)

    class_file_hs = 'data/wikidata/created/compounds_with_chebi_ids_homo_sapiens.tsv'
    count_classes_in_file(class_file_hs)

    class_file_at = 'data/wikidata/created/compounds_with_chebi_ids_arabidopsis_thaliana.tsv'
    count_classes_in_file(class_file_at)