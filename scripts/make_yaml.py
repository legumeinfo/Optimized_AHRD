import yaml
import argparse
import os
import json


def populate_blast_dbs(db_tuples):
    """Populate blast_dbs section of AHRD config."""
    dbs = {}
    for db_tuple in db_tuples:
        db_name = db_tuple["dbName"]
        diamond_result = db_tuple["blast_file"]
        db_path = db_tuple["dbPath"]

        db_section = {
            "weight": 653,
            "description_score_bit_score_weight": 2.717061,
            "file": os.path.abspath(diamond_result),
            "database": os.path.abspath(db_path),
            "blacklist": "blank1.txt",
            "filter": "blank2.txt",
            "token_blacklist": "blank3.txt",
        }

        if db_name.lower() == "uniref90":
            db_section["fasta_header_regex"] = (
                r"^>(?<accession>UniRef90_\S+)\s+(?<description>.+?)(?=\s(?:n=|OS=)|$)"
            )
            db_section["short_accession_regex"] = r"UniRef90_(?<shortAccession>[A-Z0-9]+)"
        elif db_name.lower() in ["trembl", "swissprot"]:
            db_section["fasta_header_regex"] = r"^>(?<accession>\S+)\s+(?<description>.+)"
            db_section["short_accession_regex"] = r"^\w+\|(?<shortAccession>\w+)\|"
        dbs[db_name] = db_section
    return dbs


def create_yaml(proteins_fasta, out_dir, db_tuples_json):
    """Workflow to creata an AHRD config YML file."""
    out_yaml = os.path.join("ahrd_config.yml")
    with open(db_tuples_json) as f:
        db_tuples = json.load(f)
    blast_dbs = populate_blast_dbs(db_tuples)
    raw = f"{out_dir}/interproscan_concatenated.raw"
    data = {
        "interpro_database": "/data/elavelle/databases/interpro.xml",
        "interpro_result": os.path.abspath(raw),
        "gene_ontology_result": r"/data/elavelle/databases/goa_uniprot_all.gaf",
        "reference_go_regex": r"^UniProtKB\s+(?<shortAccession>\S+)\s+\S+\s+\S+\s+(?<goTerm>GO:\d{7})",
        "proteins_fasta": os.path.abspath(proteins_fasta),
        "token_score_bit_score_weight": 0.468,
        "token_score_database_score_weight": 0.2098,
        "token_score_overlap_score_weight": 0.3221,
        "output": "./ahrd_interpro_output.csv",
        "blast_dbs": blast_dbs,
    }

    with open(out_yaml, "w") as outfile:
        yaml.dump(data, outfile, default_flow_style=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Write an AHRD YML config in the provided output directory."
    )
    parser.add_argument("proteins_fasta", type=str, help="Path to the proteins FASTA file")
    parser.add_argument("out_dir", type=str, help="Out directory")
    parser.add_argument(
        "db_tuples_json", type=str, help="JSON string: [(dbName, diamond_output_file, dbPath), ...]"
    )
    args = parser.parse_args()
    create_yaml(args.proteins_fasta, args.out_dir, args.db_tuples_json)
