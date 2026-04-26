#!/usr/bin/env python3
# coding: utf8


#define a dataclass for LUT
import pandas as pd
import os
import argparse
import sys
import numpy as np
from Bio import SeqIO
import ast
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
        logfile = os.path.join(output_dir, 'create_finaltable.log')   # <- removed leading dot
        fh = logging.FileHandler(logfile)
        fh.setFormatter(formatter)
        logger.addHandler(fh)
        logger.info(f"Logging to file: {logfile}")
    else:
        logger.info("No output directory provided; logging only to console.")

    return logger

# this function imports the lookup table
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
    table = df
    return table

# this function loads pvacbind resulrs
def importBindResults(allepitopes: str) -> pd.DataFrame:
    table = pd.read_csv(allepitopes, sep="\t")
    table.columns = table.columns.str.replace(" ", "_")
    return table

# this function loads immunogenicity data
def importImmunePeptides(immunelist: str) -> pd.DataFrame:
    table = pd.read_csv(immunelist,sep="\t")
    return table

# this function add immunogenicity information from iedb epitope to prediction table
def getEvidence(row,immuneinfo):
    peptide = row["Epitope_Seq"]
    HLA_Allele = row["HLA_Allele"]
    hla_restricted = "NA"
    matching_hla_str = "NA"
    matching_hla = []  # store all matching HLAs
    num_of_validations = 0
    if peptide in immuneinfo["Peptide_sequence"].values:
        peptide_data = immuneinfo[immuneinfo["Peptide_sequence"]==peptide]
        num_of_validations = sum(peptide_data["Number_of_papers_reporting_peptide_as_immunogenic"])
        #split_lists = peptide_data["HLA_validated_experimentaly"].str.split(",")
        split_lists = (peptide_data["HLA_validated_experimentaly"].dropna().apply(lambda x: x.split(",") if isinstance(x, str) else []))
        hlalist = [item.strip() for sublist in split_lists for item in sublist]
        #hlalist = [item for sublist in split_lists for item in sublist]
        pep_immuno_status = "High" if "High" in peptide_data["Peptide_immunogenicity"].values else "Medium" if "Medium" in peptide_data["Peptide_immunogenicity"].values else "Low" if "Low" in peptide_data["Peptide_immunogenicity"].values else "unknown"
        pep_id = ";".join(peptide_data["peptide_id"])
        virus_name = "|".join(peptide_data["Virus"])
        gene_name = "|".join(peptide_data["gene_name"])
        matching_hla = [hla for hla in hlalist if hla == HLA_Allele or hla in HLA_Allele.split(",")]
        if matching_hla:
            hla_restricted = "True"
            matching_hla_str = ",".join(matching_hla)
        else:
            hla_restricted = "False"
            matching_hla_str = "NA"
    else:
        num_of_validations = 0
        hla_restricted = "NA"
        pep_id = "NA"
        virus_name ="NA"
        gene_name = "NA"
        pep_immuno_status = "NA"
    return pep_id,virus_name,gene_name,pep_immuno_status,num_of_validations,hla_restricted,matching_hla_str

# this function generate isbinder column for prioritization
def classify_binder(row):
    if (row["Median_Percentile"] <= 1) or (row["EPITOPES: Peptide immunogenicity"] in ["High", "Medium", "Low"]):
        return 1
    else:
        return 2

# check weather the peptide is in human protome by string match    
def peptide_in_fasta(human_sequences,peptide) -> bool:
    return any(peptide in seq for seq in human_sequences)

# this function prioritize epitopes
def sortPredictions(predictiontable,virus_summary_final,luttable,human_proteome):
    
     # explode lists into rows
    df = pd.read_csv(virus_summary_final, sep="\t")
    df["gene_TPM"] = df["gene_TPM"].apply(lambda x: ast.literal_eval(x) if pd.notna(x) else [])
    gene_table = df.explode("gene_TPM")
    
    # split string into columns
    gene_table[["gene_id", "gene_symbol", "TPM"]] = gene_table["gene_TPM"].str.split(":", expand=True)
    
    # convert TPM back to float
    gene_table["TPM"] = gene_table["TPM"].astype(float)
    gene_table = gene_table.drop(columns=["gene_TPM"])

    #get lut table
    luttable = importluttable(lutfile=luttable, filetype="virus")
    luttable = luttable[["gene_id","protein_id"]]

    #get protein id
    genelist_final = gene_table.merge(luttable, how="left",on="gene_id")
    #the table has gene_id,gene_symbol,virus,TPM,protein_id,confidence,viral_score,viral_count
    genelist_final = genelist_final[genelist_final["confidence"]=="HC"]

    genelist_final = genelist_final[["protein_id","gene_symbol","gene_id","TPM","virus"]]
    genelist_final=genelist_final.rename(columns={
        "protein_id":"Mutation",
        "gene_symbol":"gene_symbol",
        "gene_id":"gene_id",
        "TPM":"Expression_TPM",
        "virus": "virus"
    })

    combined_results = predictiontable.merge(
    genelist_final,
    on="Mutation",
    how="left"
    )
    #print(combined_results)
    #create immunogenicity rank for sorting
    immunogenicity = {'High': 0, 'Medium': 1, 'Low': 2, "NA": 4}
    hla_restricted = {'True': 0, 'False': 1, "NA": 3}
    
    human_sequences = [str(record.seq) for record in SeqIO.parse(human_proteome, "fasta")]


    #get frequency of peptide into a new column
    #print(combined_results['EPITOPES: Peptide immunogenicity'])
    combined_results["nb_additonal_alleles"] = combined_results["Additional_HLA_Alleles"].apply(lambda x: len([a for a in str(x).split(",") if a != "NA"]))
    combined_results['pep_immuno_sorting'] = combined_results['EPITOPES: Peptide immunogenicity'].map(immunogenicity)
    combined_results['pep_hla_restricted'] = combined_results['EPITOPES: Patient immunogenic HLA'].map(hla_restricted)
    combined_results["Isbinder"] = combined_results.apply(lambda row: classify_binder(row), axis=1)
    # apply peptide search directly on column
    combined_results["exist_in_human"] = combined_results["Epitope_Seq"].apply(lambda pep: peptide_in_fasta(human_sequences, pep))
    combined_results_virus = combined_results[(combined_results["exist_in_human"]== False)]
    #sort predictions
    combined_results_sorted_all = combined_results_virus.sort_values(
        by=['Isbinder','pep_hla_restricted','pep_immuno_sorting','EPITOPES: Number validations','Expression_TPM','Median_Percentile',"nb_additonal_alleles"], 
        ascending=[True,True,True,False,False,True,False]
    )
    return combined_results_sorted_all

# this function create the final results by combining and sorting prediction with immunogenicity info
def createFinalTable(args):

    #load predictions from pvacind and select colums and apply filters
    predictions = importBindResults(args.allepitopes)
    predictions = predictions[(predictions["difficult_n_terminal_residue"] == False) &
                               (predictions["c_terminal_cysteine"] == False) & 
                               (predictions["c_terminal_proline"] == False) & 
                               (predictions["n_terminal_asparagine"] == False)]
    predictions = predictions[["Mutation","HLA_Allele","Sub-peptide_Position",
                               "Epitope_Seq","Median_IC50_Score","Median_Percentile"]]
    predictions = predictions[(predictions["Median_Percentile"] < 2)]
    predictions_sorted = predictions.sort_values(by="Median_Percentile",ascending=True)

    # retain data of best allele for a peptide (lowest median percentile)
    # summarize othwr valid alleles ad additional alleles ( rank <2)
    final_prediction_best = (predictions_sorted
    .loc[predictions_sorted["Epitope_Seq"].drop_duplicates().index]
    )
    final_prediction_add_alleles = predictions_sorted.groupby("Epitope_Seq", as_index=False).agg({
    "HLA_Allele": lambda x: ",".join(x.iloc[1:]) if len(x) > 1 else "NA"  # skip first
    }).rename(columns={"HLA_Allele": "Additional_HLA_Alleles"})
    final_predictions = final_prediction_best.merge(final_prediction_add_alleles,on="Epitope_Seq",how="left")
    
    
    #load immunogenicity information
    immunedata = importImmunePeptides(args.immunelist)
    #print(immunedata.head())
    # Apply getEvidence row-wise and expand multiple outputs
    final_predictions[["EPITOPES: Peptide ID",
                 "EPITOPES: Virus",
                 "EPITOPES: Gene name",
                 "EPITOPES: Peptide immunogenicity",
                 "EPITOPES: Number validations",
                 "EPITOPES: Patient immunogenic HLA",
                 "EPITOPES: Matching Restricted HLA"]] = \
        final_predictions.apply(lambda row: getEvidence(row, immunedata), axis=1, result_type='expand')
    

    sorted_predictions = sortPredictions(predictiontable=final_predictions,
                                         virus_summary_final=args.virus_summary_final,
                                         luttable=args.luttable,
                                         human_proteome=args.human_fasta)
    sorted_outfile = os.path.join(args.output_dir, "final_prioritized_predictions.tsv")
    predictions_outfile = os.path.join(args.output_dir, "annotated_predictions.tsv")
    sorted_predictions.to_csv(sorted_outfile, sep="\t", index=False)
    final_predictions.to_csv(predictions_outfile, sep="\t", index=False)


# Instantiate the parser
parser = argparse.ArgumentParser(prog='create_finaltable.py', formatter_class=argparse.RawTextHelpFormatter)
subparsers = parser.add_subparsers(help='sub-command help')

createFinalTable_parser = subparsers.add_parser("createFinalTable", description="combine pvacbind predictions with expression and immunogenicity information and give prioritized table")
createFinalTable_parser.add_argument("-allepitopes", required=True,type=str, help='<Required> combined all epitope file from pvacbind')
createFinalTable_parser.add_argument("-luttable", required=True, type=str, help='<Required> LUT file')
createFinalTable_parser.add_argument("-immunelist", required=True, type=str, help='<Required> immunogenicity information table')
createFinalTable_parser.add_argument("-virus_summary_final", required=True, type=str, help='<Required> viral summary final results')
createFinalTable_parser.add_argument("-human_fasta", required=True, type=str, help='<Required> human fasta file of proteome')
createFinalTable_parser.add_argument("-output_dir", required=True, type=str, help='<Required> path of folder where to write output file')
createFinalTable_parser.set_defaults(func=createFinalTable)

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


    
