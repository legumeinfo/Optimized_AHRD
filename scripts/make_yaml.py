import yaml
import argparse
import os
import json
import argparse

def create_yaml(proteins_fasta, interpro_result, db_tuples_json):
    out_yaml = os.path.join("ahrd_config.yml")

    db_tuples = json.loads(db_tuples_json)

    # Generalized regex, may need changing based on .fasta headers
    fasta_regex = r'^>(?<accession>\S+)\s+(?<description>.+?)(?=\s(?:n=|OS=)|$)'

    blast_dbs = {}
    for db_tuple in db_tuples:
        db_name = db_tuple["dbName"]
        diamond_result = db_tuple["blastFile"]
        db_path = db_tuple["dbPath"]

        blast_dbs[db_name] = {
            'weight': 653,
            'description_score_bit_score_weight': 2.717061,
            'file': os.path.abspath(diamond_result),
            'database': os.path.abspath(db_path),
            'blacklist': 'blank1.txt',
            'filter': 'blank2.txt',
            'token_blacklist': 'blank3.txt',
        }

        if db_name.lower() == "uniref90":
            db_section["fasta_header_regex"] = r"^>(?<accession>UniRef90_\S+)\s+(?<description>.+?)(?=\s(?:n=|OS=)|$)"
            db_section["short_accession_regex"] = r"UniRef90_(?<shortAccession>[A-Z0-9]+)"


    data = {
        'interpro_result': os.path.abspath(interpro_result),
        'gene_ontology_result': '/data/elavelle/databases/goa_uniprot_all.gaf',
        'reference_go_regex': '^UniProtKB\s+(?<shortAccession>\S+)\s+\S+\s+\S+\s+(?<goTerm>GO:\d{7})',
        'proteins_fasta': os.path.abspath(proteins_fasta),
        'token_score_bit_score_weight': 0.468,
        'token_score_database_score_weight': 0.2098,
        'token_score_overlap_score_weight': 0.3221,
        'output': './ahrd_interpro_output.csv',
        'blast_dbs': blast_dbs,
    }

    with open(out_yaml, 'w') as outfile:
        yaml.dump(data, outfile, default_flow_style=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create a YAML file for ahrd to run")
    parser.add_argument('proteins_fasta', type=str, help="Path to the proteins FASTA file")
    parser.add_argument('interpro_result', type=str, help="Path to the concatenated InterProScan result")
    parser.add_argument('db_tuples_json', type=str, help="JSON string: [(dbName, diamond_output_file, dbPath), ...]")

    args = parser.parse_args()

    create_yaml(args.proteins_fasta, args.interpro_result, args.db_tuples_json)
