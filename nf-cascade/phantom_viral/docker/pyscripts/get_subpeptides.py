#!/usr/bin/env python3
# coding: utf8

import sys
from Bio import SeqIO
import argparse
import logging
import os

# logging function
def setup_logging(output_dir=None, log_level=logging.INFO):
    logger = logging.getLogger()
    logger.setLevel(log_level)

    # Prevent duplicate handlers if called multiple times
    if logger.hasHandlers():
        logger.handlers.clear()

    # Formatter
    formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s',
                                  datefmt='%Y-%m-%d %H:%M:%S')

    # Console logging
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(formatter)
    logger.addHandler(sh)

    # File logging
    if output_dir:
        output_dir = os.path.abspath(output_dir)
        os.makedirs(output_dir, exist_ok=True)
        logfile = os.path.join(output_dir, 'get_subpeptides.log')   # <- removed leading dot
        fh = logging.FileHandler(logfile)
        fh.setFormatter(formatter)
        logger.addHandler(fh)
        logger.info(f"Logging to file: {logfile}")
    else:
        logger.info("No output directory provided; logging only to console.")

    return logger

# check weather the peptide is in human protome by string match    
def peptide_in_fasta(human_sequences,peptide) -> bool:
    return any(peptide in seq for seq in human_sequences)


def parse_lengths(lengths_str):
    """Convert comma-separated string of lengths into a list of ints."""
    return [int(x.strip()) for x in lengths_str.split(",") if x.strip().isdigit()]

def parse_problematic_residues(residues_str):
    """Convert comma-separated residues into a set."""
    return set(residues_str.replace(" ", "").split(",")) if residues_str else set()

def generate_subpeptides(seq, record_id, lengths, human_sequences, problematic_n=set(), problematic_c=set()):
    """Yield subpeptides, skipping those with problematic residues at N- or C-terminal."""

    for L in lengths:
        if len(seq) >= L:
            for i in range(len(seq) - L + 1):
                pep = seq[i:i+L]
                #skip if terminal aa is problematic
                if problematic_n and pep[0] in problematic_n:
                    continue
                if problematic_c and pep[-1] in problematic_c:
                    continue
                #skip if exist in human proteome
                if peptide_in_fasta(human_sequences, pep):
                    continue          
                start = i + 1
                end = i + L
                index = f"{record_id}_{start}_{end}_{L}"
                yield record_id, pep, start, end, L, index

#this function creates peptide fasta file and table for MHC class I and Class II lengths
def fasta_to_subpeptides(args):
    
    input_fasta=args.input_fasta
    output_fasta1=args.output_fasta1
    output_fasta2=args.output_fasta2
    output_table=args.output_table
    lengths1 = parse_lengths(args.lengths_str1)
    lengths2 = parse_lengths(args.lengths_str2)
    all_lengths = set(lengths1 + lengths2)
    problematic_n = parse_problematic_residues(args.problematic_n_str)
    problematic_c = parse_problematic_residues(args.problematic_c_str)

    #load human proteome 
    human_sequences = [str(record.seq) for record in SeqIO.parse(args.human_proteome, "fasta")]

    with open(output_fasta1, "w") as fasta1, open(output_fasta2, "w") as fasta2, open(output_table, "w") as table_out:
        table_out.write("OriginalID\tSubpeptide\tStart\tEnd\tLength\tIndex\n")
        
                
        for record in SeqIO.parse(input_fasta, "fasta"):
            seq = str(record.seq)
            
            # mhc1 peptides
            for orig_id, pep, start, end, L, idx in generate_subpeptides(
                    seq, record.id, lengths1, human_sequences, problematic_n, problematic_c):
                fasta1.write(f">{idx}\n{pep}\n")
                table_out.write(f"{orig_id}\t{pep}\t{start}\t{end}\t{L}\t{idx}\n")
            
            # mhc2 peptides
            for orig_id, pep, start, end, L, idx in generate_subpeptides(
                    seq, record.id, lengths2, human_sequences, problematic_n, problematic_c):
                fasta2.write(f">{idx}\n{pep}\n")
                table_out.write(f"{orig_id}\t{pep}\t{start}\t{end}\t{L}\t{idx}\n")

# Instantiate the parser
parser = argparse.ArgumentParser(prog='get_subpeptides.py', formatter_class=argparse.RawTextHelpFormatter)
subparsers = parser.add_subparsers(help='sub-command help')

fasta_to_subpeptides_parser = subparsers.add_parser("fasta_to_subpeptides", description="generate peptide list for running mixmhcpred and prime")
fasta_to_subpeptides_parser.add_argument("-input_fasta", required=True,type=str, help='fasta file of viral proteins identified')
fasta_to_subpeptides_parser.add_argument("-output_fasta1", required=True, type=str, help='output fasta file for running PRIME for MHC class I allle')
fasta_to_subpeptides_parser.add_argument("-output_fasta2", required=True, type=str, help='output fasta file for running PRIME for MHC class II allle')
fasta_to_subpeptides_parser.add_argument("-output_table", required=True, type=str, help='output table of all the peptides and location')
fasta_to_subpeptides_parser.add_argument("-human_proteome", required=True, type=str, help='<Required> human fasta file of proteome')
fasta_to_subpeptides_parser.add_argument("-lengths_str1", type=str, default='8,9,10,11', help='comma separated length for class I allele')
fasta_to_subpeptides_parser.add_argument("-lengths_str2", type=str, default='12,13,14,15,16,17,18', help='comma separated length for class II allele')
fasta_to_subpeptides_parser.add_argument("-problematic_n_str", default="", type=str, help='Problematic amino acids at n terminal')
fasta_to_subpeptides_parser.add_argument("-problematic_c_str", default="", type=str, help='Problematic amino acids at c terminal')
fasta_to_subpeptides_parser.set_defaults(func=fasta_to_subpeptides)

if __name__ == "__main__":
    args = parser.parse_args()
    if hasattr(args, "func"):
        # initialize logging now that we have args (so user-specified output_dir is used)
        # fallback to current working dir if output_dir not provided
        outdir_for_logs = getattr(args, "output_dir", None) or os.getcwd()
        logger = setup_logging(outdir_for_logs)

        try:
            args.func(args)
        except Exception as e:
            # log exception and exit non-zero for pipelines
            logger.exception("Unhandled exception during execution")
            sys.exit(2)
    else:
        parser.print_help()  