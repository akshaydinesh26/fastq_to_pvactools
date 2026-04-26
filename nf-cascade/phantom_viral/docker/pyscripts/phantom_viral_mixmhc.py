#!/usr/bin/env python3
# coding: utf8


#define a dataclass for LUT
import pandas as pd
import os
import argparse
import sys
import numpy as np
from Bio import SeqIO
import logging


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
        logfile = os.path.join(output_dir, 'phantom_viral.log')   # <- removed leading dot
        fh = logging.FileHandler(logfile)
        fh.setFormatter(formatter)
        logger.addHandler(fh)
        logger.info(f"Logging to file: {logfile}")
    else:
        logger.info("No output directory provided; logging only to console.")

    return logger

# load LUT table
def importluttable(lutfile: str, filetype: str) -> pd.DataFrame:

    # the type can be virus or human
    if filetype == "human":
        columns = ["gene_id","gene_symbol","gene_biotype","median_transcript_length_kb","strand",
                   "chromosome","start","end","transcript_id","transcript_length","original_gene_id"]
        dtypes = {
            "gene_id": str,
            "gene_symbol": str,
            "gene_biotype": str,
            "median_transcript_length_kb": float,
            "strand": str,
            "chromosome": str,
            "start": int,
            "end": int,
            "transcript_id": str,
           "transcript_length": str,
            "original_gene_id": str
            }
    else:
        columns = ["gene_id","gene_symbol","gene_biotype","median_transcript_length_kb","strand","chr",
                    "start","end","transcript_id","transcript_length_kb","protein_id","virus"]
        dtypes = {
            "gene_id": str,
            "gene_symbol": str,
            "gene_biotype": str,
            "median_transcript_length_kb": float,
            "strand": str,
            "chr": str,
            "start": int,
            "end": int,
            "transcript_id": str,
           "transcript_length_kb": str,
           "protein_id": str,
           "virus": str
            }


    df = pd.read_csv(lutfile, sep="\t", header=None, names=columns,skiprows=1,dtype=dtypes)
    table = df.set_index("gene_id")
    return table

# function parsing star count and record total mapped and unmapped reads
def parseSTARcount(count_matrix: str, stranded: str) -> pd.DataFrame:
    columns = ["gene_id","count_1","count_2","count_3"]
    TotUnmappedReads, TotMappedReads = 0, 0

    dtypes = {
        "gene_id": str,
        "count_1": int,
        "count_2": int,
        "count_3": int
    }

    #Star has 4 colum output strandness - geneid,"no", "yes","reverse"
    count_column = {"no" : "count_1", "yes" : "count_2", "reverse" : "count_3"}
    column = {"no": 1, "yes": 2, "reverse": 3}
    df = pd.read_csv(count_matrix,skiprows=4,header=None,dtype=dtypes,names=columns,sep="\t")
    df = df.set_index("gene_id")
    col = count_column[stranded]
    table = df[[col]].rename(columns={col: "count"})

    # get total mappeed and unmapped reads
    with open(count_matrix, "r") as f:
        for line in f:
            line = line.rstrip().split("\t")
            if line[0].startswith("N_"):
                TotUnmappedReads += int(line[column[stranded]])  # Sum up unmapped/ambiguous/multimapping/noFeature reads
                continue
            
            gid, counts = line[0], int(line[column[stranded]])  # Import data from line and select correct column
            TotMappedReads += counts  # Sum up total of reads counted (for TPM)
    return table, TotMappedReads, TotUnmappedReads
    
def classify_confidence(score, threshold):
    if score == 0:
        return "NE"  # or keep empty
    elif score >= threshold:
        return "HC"
    else:
        return "LC"

#Function to estimate Human RPK for downstream from STAR counts of human alignment
# arguments args.luttable,args.counts,args.stranded = yes,no,reverse 
def QuantifyHumanRPKSUM(count_matrix: str,luttable: str, stranded: str):
    RPKSUM = 0

    #get gene ids from lut table
    gene_table = importluttable(lutfile=luttable, filetype="human")
    gene_list = gene_table.index.to_numpy()

    #get star count
    df, TotMappedReads, TotUnmappedReads= parseSTARcount(count_matrix=count_matrix, stranded=stranded)
    #df = df.reset_index()
    #df["gene_id"] = df["gene_id"].str.split(".").str[0]
    #df = df.set_index("gene_id")
    hitlist = df.index.isin(gene_list)
    human_count = df.loc[hitlist, :]
    human_count_final = human_count.join(gene_table, how="left")
    human_count_final["RPK"] = (human_count_final["count"] / human_count_final["median_transcript_length_kb"])
    RPKSUM = human_count_final["RPK"].sum()
    return RPKSUM


#Function to get viral protein fasta file of genes which are found
# arguments args.luttable,args.counts,args.stranded = yes,no,reverse, args.humancounts
def GetVirusFasta(viral_summary,protein_fasta,luttable,output_dir):
    
    # explode lists into rows
    gene_table = viral_summary.explode("gene_TPM")
    
    # split string into columns
    gene_table[["gene_id", "gene_symbol", "TPM"]] = gene_table["gene_TPM"].str.split(":", expand=True)
    
    # convert TPM back to float
    gene_table["TPM"] = gene_table["TPM"].astype(float)
    gene_table = gene_table[gene_table["TPM"] > 0].copy()
    gene_table = gene_table.drop(columns=["gene_TPM"]).set_index("gene_id")
    
    #get lut table
    luttable = importluttable(lutfile=luttable, filetype="virus")
    luttable = luttable[["protein_id"]]

    #get protein id
    genelist_final = gene_table.join(luttable, how="left").reset_index()
    #the table has gene_id,gene_symbol,virus,TPM,protein_id,confidence,viral_score,viral_count
    genelist_final = genelist_final[genelist_final["confidence"]=="HC"]
    genelist_final = genelist_final[genelist_final["virus"] != "NC_022518.1|Human_endogenous_retrovirus_K113"]
    genelist_final = genelist_final[genelist_final["virus"] != "ref|NC_022518.1|Human_endogenous_retrovirus_K113_complete_genome|gi|548558394"]
    
     # check if the gene list is empty
    if genelist_final.empty:
        outfile = os.path.join(output_dir, "viral_genes_found.fasta")
        # write an empty FASTA
        with open(outfile, "w") as f:
            pass
        logger.info("No high-confidence virus found — wrote empty FASTA and exiting.")
        sys.exit(0)

    # collect all protein IDs we need
    target_proteins = set(genelist_final["protein_id"].dropna().astype(str))

    # make dictionary for metadata lookup
    metadata = genelist_final.set_index("protein_id")[["virus", "confidence", "viral_score", "viral_count"]].to_dict("index")

    selected_records = []
    for record in SeqIO.parse(protein_fasta, "fasta"):
        prot_id = record.id.split()[0]
        if prot_id in target_proteins:
            # fetch metadata
            meta = metadata[prot_id]
            # rewrite header with extra info
            record.description = f"{prot_id} | {meta['virus']} | confidence={meta['confidence']} | score={meta['viral_score']:.4f} | count={meta['viral_count']}"
            selected_records.append(record)
    
    # you can either return them or write to file
    outfile = os.path.join(output_dir, "viral_genes_found.fasta")
    SeqIO.write(selected_records, outfile, "fasta")
    logger.info(f"Wrote {len(selected_records)} protein sequences to {outfile}")

#Function to generate viral expression matrix
# arguments args.luttable,args.counts,args.stranded = yes,no,reverse, args.humancounts
def QuantifyVirusCounts(args):
    RPKSUM = QuantifyHumanRPKSUM(count_matrix=args.humancounts,luttable=args.human_luttable,stranded=args.stranded)
    logger.info(RPKSUM)
    TotReads, confidence_threshold = 0, 0.165352916

    #get gene ids from lut table
    gene_table = importluttable(lutfile=args.luttable, filetype="virus")
    gene_list = gene_table.index.to_numpy()

    # get viral count to table and estimate rpk, tpm, cpm and viral score
    df, TotMappedReads, TotUnmappedReads = parseSTARcount(count_matrix=args.counts, stranded=args.stranded)
    hitlist = df.index.isin(gene_list)
    viral_count = df.loc[hitlist, :]
    #logger.info(viral_count.loc["1489079"])
    viral_count_final = viral_count.join(gene_table, how="left")
    viral_count_final["RPK"] = (viral_count_final["count"] / viral_count_final["median_transcript_length_kb"])
    #logger.info(viral_count_final.loc["1489079"])
    RPKSUM += viral_count_final["RPK"].sum()
    logger.info(RPKSUM)
    # RPKSUM consider all mapped reads in humans as well as virus
    # get total viral reads (sum of mapped and unmapped reads in viral db)
    TotReads += TotMappedReads + TotUnmappedReads
    logger.info(f"Total mapped Reads: {TotMappedReads}")
    logger.info(f"Total unmapped reads: {TotUnmappedReads}")
    logger.info(f"Total reads: {TotReads}")
    logger.info(f"Total RPK: {RPKSUM}")

    # per million values for RPKSUM and TotReads
    RPKSUM_pm, TotReads_pm = RPKSUM / 1000000, TotReads / 1000000
    # Prevent division by zero if no alignment
    TotReads, RPKSUM_pm, TotReads_pm = max(TotReads, 1), max(RPKSUM_pm, 1), max(TotReads_pm, 1)

    viral_count_final["RPKM"] = (viral_count_final["RPK"].astype(float) / TotReads_pm)
    viral_count_final["TPM"] = (viral_count_final["RPK"].astype(float) / RPKSUM_pm)
    #logger.info(viral_count_final.loc["1489079"])
    viral_count_final["CPM"] = ((viral_count_final["count"].astype(float) / TotReads) / viral_count_final["median_transcript_length_kb"].astype(float)) * 1000000

    # Reset index so 'gene_id' becomes a column
    viral_count_final = viral_count_final.reset_index()

    # Create per-gene TPM string
    viral_count_final["gene_TPM_str"] = viral_count_final.apply(
        lambda row: f"{row['gene_id']}:{row['gene_symbol']}:{row['TPM']:.2f}", axis=1
        )
    
    viral_count_final = viral_count_final.sort_values(by='TPM', ascending=False)
    # Group by virus to summarize and collect gene TPMs
    viral_summary = viral_count_final.groupby("virus").agg(
        viral_count=("count", "sum"),
        num_genes=("gene_id", "nunique"),
        total_TPM=("TPM", "sum"),
        total_CPM=("CPM", "sum"),
        gene_TPM=("gene_TPM_str", list)  # list of gene_id:gene_symbol:TPM per virus
        ).reset_index()
    
    # define threshold for HC, LC, NE - based on viral score
    #for stranded HC >/= 0.165352916, LC < 0.165352916, NE ==0
    #for unstranded the threshold is 0.8856791
    #if args.stranded == "no":
    #    confidence_threshold = 0.8856791
    
    viral_summary["viral_score"] = (viral_summary["total_CPM"].astype(float) / viral_summary["num_genes"].astype(int))
    viral_summary["confidence"] = viral_summary["viral_score"].apply(lambda x: classify_confidence(x, confidence_threshold))
    viral_summary_final = viral_summary[["virus","viral_score","confidence","viral_count","gene_TPM"]].copy()
    
    #sort by confidence
    virus_order = {'HC': 0, 'LC': 1, 'NE': 2}
    viral_summary_final['virus_rank'] = viral_summary_final['confidence'].map(virus_order)
    viral_summary_final = viral_summary_final.sort_values(by=['virus_rank','viral_score'], ascending=[True,False]).drop("virus_rank",axis=1)
    #write output to file
    outfile = os.path.join(args.output_dir, "identified_virus.txt")
    viral_summary_final.to_csv(outfile, sep="\t", index=False)

    #get fasta file of proteins expressing
    GetVirusFasta(viral_summary_final,args.protein_fasta,args.luttable,args.output_dir)
    return viral_summary_final

# Instantiate the parser
parser = argparse.ArgumentParser(prog='phantom_viral.py', formatter_class=argparse.RawTextHelpFormatter)
subparsers = parser.add_subparsers(help='sub-command help')

QuantifyVirusCounts_parser = subparsers.add_parser("QuantifyVirusCounts", description="Quantify virus gene expression and calculates confidence of virus infection from STAR")
QuantifyVirusCounts_parser.add_argument("-counts", required=True, type=str, help='<Required> gene count file (STAR counts) for virus')
QuantifyVirusCounts_parser.add_argument("-luttable", required=True, type=str, help='<Required> LUT file (created with GTFToLUT)')
QuantifyVirusCounts_parser.add_argument("-human_luttable", required=True, type=str, help='<Required> LUT file (created with GTFToLUT)')
QuantifyVirusCounts_parser.add_argument("-stranded", required=True, type=str, help='<Required> strandness of protocole for RNAseq [no, yes or reverse]')
QuantifyVirusCounts_parser.add_argument("-humancounts", required=True, type=str, help='<Required> gene count file (STAR counts) for human')
QuantifyVirusCounts_parser.add_argument("-output_dir", required=True, type=str, help='<Required> path of folder where to write output file')
QuantifyVirusCounts_parser.add_argument("-protein_fasta", required=True, type=str, help='<Required> virus protein fasta file')
QuantifyVirusCounts_parser.set_defaults(func=QuantifyVirusCounts)

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
