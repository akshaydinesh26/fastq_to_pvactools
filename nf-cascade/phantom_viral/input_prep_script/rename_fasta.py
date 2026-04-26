#!/usr/bin/env python3
import sys
import argparse

def load_mapping(tsv_file):
    mapping = {}
    with open(tsv_file) as f:
        for line in f:
            if line.strip() == "":
                continue
            old, new = line.strip().split("\t")[:2]
            mapping[old] = new
    return mapping

def rename_fasta(args):

    fasta_file = args.genome_fasta
    mapping = load_mapping(args.rename_tsv)
    output_file = args.output_fasta

    with open(fasta_file) as fin, open(output_file, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                header = line[1:].strip().split()[0]  # first word after ">"
                rest = line.strip()[len(header)+1:]   # keep rest of header
                new_header = mapping.get(header, header)
                fout.write(f">{new_header}\n")
            else:
                fout.write(line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert fasta header to workflow compatible")
    parser.add_argument("genome_fasta", type=str, help="<Required> input fasta")
    parser.add_argument("rename_tsv", type=str, help="<Required> input table to rename as tsv with id new id")
    parser.add_argument("output_fasta", type=str, help="<Required> output rename fasta file name")
    args = parser.parse_args()
    rename_fasta(args)