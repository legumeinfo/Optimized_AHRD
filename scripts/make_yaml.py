import yaml
import argparse
import os

def create_yaml(proteins_fasta, diamond_cat, interpro_result, database, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    out_dir = os.path.abspath(out_dir)
    out_yaml = os.path.join(out_dir, "ahrd_config.yml")

    data = {
        'proteins_fasta': proteins_fasta,
        'token_score_bit_score_weight': 0.468,
        'token_score_database_score_weight': 0.2098,
        'token_score_overlap_score_weight': 0.3221,
        'output': './ahrd_interpro_output.csv',
        'blast_dbs': {
            'diamond': {
                'weight': 653,
                'description_score_bit_score_weight': 2.717061,
                'file': diamond_cat,
                'database': database,  # Dynamically filled with the 'database' argument
                'blacklist': 'blank1.txt',
                'filter': 'blank2.txt',
                'token_blacklist': 'blank3.txt',
                'fasta_header_regex': '^>(?<accession>\\S+)\\s+(?<description>.+?)\\s+\\[(?<organism>.+?)\\]$'
            }
        },
        'interpro_result': interpro_result
    }

    with open(out_yaml, 'w') as outfile:
        yaml.dump(data, outfile, default_flow_style=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create a YAML file for ahrd to run")
    parser.add_argument('proteins_fasta', type=str, help="Path to the proteins FASTA file")
    parser.add_argument('diamond_file', type=str, help="Path to the concatenated Diamond output file")
    parser.add_argument('interpro_result', type=str, help="Path to the concatenated InterProScan result")
    parser.add_argument('database', type=str, help="Path to the Diamond database file")
    parser.add_argument('out_dir', type=str, help="Directory to output the .yml file")

    args = parser.parse_args()

    create_yaml(args.proteins_fasta, args.diamond_file, args.interpro_result, args.database, args.out_dir)
