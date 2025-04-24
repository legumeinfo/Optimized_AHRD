import os
import argparse
from Bio import SeqIO


def chunk_fasta(fasta, out_dir, chunk_size):
    """Takes a FASTA file as input and outputs chunks of chunk_size"""
    count = 1
    chunk_number = 1
    out_dir = os.path.abspath(out_dir)
    os.makedirs(f"{out_dir}/chunks", exist_ok=True)

    chunk = open(f"{out_dir}/chunks/{chunk_number}.fasta", "w")
    with open(fasta) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            # print(record.id)
            if count > chunk_size:
                chunk.close()
                count = 0
                chunk_number += 1
                chunk = open(f"{out_dir}/chunks/{chunk_number}.fasta", "w")
            record = f">{record.id}\n{record.seq}"
            chunk.write(record + "\n")
            count += 1
    chunk.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse a FASTA file and chunk at chunk size.")
    parser.add_argument("fasta", type=str, help="Path to the proteins FASTA file")
    parser.add_argument("out_dir", type=str, help="Out directory")
    parser.add_argument(
        "chunk_size", type=str, help="JSON string: [(dbName, diamond_output_file, dbPath), ...]"
    )
    args = parser.parse_args()
    chunk_fasta(args.fasta, args.out_dir, args.chunk_size)
