import yaml
import argparse
import os
import json
import argparse
from pathlib import Path

def create_yaml(proteins_fasta, out_dir, db_tuples_json, gaf_path, ips_xml_path, desc_blacklist, token_blacklist):
    out_yaml = os.path.join("ahrd_config.yml")
    with open(db_tuples_json) as f:
        db_tuples = json.load(f)
    #db_tuples = json.loads(db_tuples_json)
    
    # Generalized regex, may need changing based on .fasta headers
    # fasta_regex = r'^>(?<accession>\S+)\s+(?<description>.+?)(?=\s(?:n=|OS=)|$)'

    blast_dbs = {}
    for db_tuple in db_tuples:
        db_name = db_tuple["dbName"]
        diamond_result = db_tuple["blast_file"]
        db_path = db_tuple["dbPath"]
         
        db_section = {
            'weight': 653,
            'description_score_bit_score_weight': 2.717061,
            'file': os.path.abspath(diamond_result),
            'database': os.path.abspath(db_path),
            'blacklist': os.path.abspath(desc_blacklist),
            'token_blacklist': os.path.abspath(token_blacklist),
        }
        print(db_name)
        if db_name.lower() == "uniref90":
            db_section["fasta_header_regex"] = r"^>(?<accession>UniRef90_\S+)\s+(?<description>.+?)(?=\s(?:n=|OS=)|$)"
            db_section["short_accession_regex"] = r"UniRef90_(?<shortAccession>[A-Z0-9]+)"
        
        elif db_name.lower() in ["trembl", "swissprot"]:
            #db_section["fasta_header_regex"] = r"^>(?<accession>\S+)\s+(?<description>.+?)(?=\s(?:n=|OS=)|$)"
            #db_section["fasta_header_regex"] = r"^>(?<accession>\S+)\s+(?<description>.+?)(?=\s+(?:OS=|n=)|$)"
            #db_section["fasta_header_regex"] = r"^>(?<accession>\\S+)\\s+(?<description>Uncharacterized protein TEST)"
            #works best without regex
            db_section["short_accession_regex"] = r"^\w+\|(?<shortAccession>\w+)\|"
        
        elif db_name.lower() == "glyma_refseq":
            db_section["fasta_header_regex"] = r"^>(?<accession>\S+)\s+(?<description>.+?)(?=\s*\[|$)"
            db_section["short_accession_regex"] = r"^\s*(?<shortAccession>\S+)"

        elif db_name.lower() == "medtr_lis":
            db_section["fasta_header_regex"] = r"^>(?<accession>\S+).*?def=(?<description>.+)$"
            db_section["short_accession_regex"] = r"^(?<shortAccession>\S+)"
        
        blast_dbs[db_name] = db_section

    data = {
        #'interpro_database': '/data/elavelle/databases/interpro.xml',
        'interpro_database' : os.path.abspath(ips_xml_path), 
        'interpro_result' : os.path.abspath(os.path.join(out_dir, 'interproscan_concatenated.raw')),
        #'gene_ontology_result': '/data/elavelle/databases/goa_uniprot_all.gaf',
        'gene_ontology_result' : os.path.abspath(gaf_path),
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
    parser.add_argument('out_dir', type=str, help="Out directory")
    parser.add_argument('db_tuples_json', type=str, help="JSON string: [(dbName, diamond_output_file, dbPath), ...]")
    parser.add_argument('gaf', type=str, help="Path to ontology .gaf")
    parser.add_argument('ips_xml', type=str, help="Path to interproscan .xml")
    parser.add_argument('desc_blacklist', type=str, help="Blacklist provided by ahrd")
    parser.add_argument('token_blacklist', type=str, help="Blacklist provided by ahrd")

    args = parser.parse_args()

    create_yaml(args.proteins_fasta, args.out_dir, args.db_tuples_json, args.gaf, args.ips_xml, args.desc_blacklist, args.token_blacklist)
