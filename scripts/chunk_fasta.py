import os
import sys
from Bio import SeqIO


def chunk_fasta(fasta, out_dir, chunk_size):
    """takes fasta file as input and outputs chunks of chunk_size"""
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
    if len(sys.argv) < 4:
        print(f"{sys.argv[0]}: <fasta file> <output directory> <chunk size>")
        sys.exit(1)
    fasta = sys.argv[1]
    out_dir = sys.argv[2]
    chunk_size = int(sys.argv[3])
    chunk_fasta(fasta, out_dir, chunk_size)
