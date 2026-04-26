#!/usr/bin/env python3
# coding: utf8

import pandas as pd
import os
import argparse
import re
import sys
from statistics import mean
import ast
import logging

# Logging
def setup_logging(output_dir=None, log_level=logging.INFO):
    logger = logging.getLogger()
    logger.setLevel(log_level)

    if logger.hasHandlers():
        logger.handlers.clear()

    formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s',
                                  datefmt='%Y-%m-%d %H:%M:%S')

    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(formatter)
    logger.addHandler(sh)
    return logger

# this function loads lookup table
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

# this function counts the additional alleles in binding prediction result
def count_additional_alleles(row):
    alleles = []
    for col in ["other_significant_alleles", "other_significant_alleles_PRIME"]:
        val = row[col]
        if pd.isna(val) or val == "NA":
            continue
        alleles.extend(str(val).split(","))
    # Remove empty strings and duplicates
    return len({a for a in alleles if a})


#  step1 - this function create a final peptide list by merging peptide table and virus information with gene expression
def create_peptide_final_list(args):

    #get virus and gene expression
    df = pd.read_csv(args.virus_summary_final, sep="\t")
    df["gene_TPM"] = df["gene_TPM"].apply(lambda x: ast.literal_eval(x) if pd.notna(x) else [])
    gene_table = df.explode("gene_TPM")
    gene_table[["gene_id", "gene_symbol", "TPM"]] = gene_table["gene_TPM"].str.split(":", expand=True)
    gene_table["TPM"] = gene_table["TPM"].astype(float)
    gene_table = gene_table.drop(columns=["gene_TPM"])

    # get gene and virus information - LUT table
    virus_gene_info = importluttable(lutfile=args.luttable, filetype="virus")
    peptide_info = pd.read_csv(args.peptide_table, sep="\t")

    peptide_info["Length"] = peptide_info["Length"].astype(int)
    peptide_info["Peptide_Class"] = peptide_info["Length"].apply(lambda x: "HLA_I" if x < 11 else "HLA_II")

    virus_gene_info_subset = virus_gene_info[["protein_id","gene_id","gene_symbol","virus"]]
    merged = peptide_info.merge(virus_gene_info_subset, how="left", left_on="OriginalID", right_on="protein_id")

    # merge and generate final output
    merged = merged.merge(gene_table[["gene_id","TPM"]], how="left", on="gene_id")
    merged = merged.rename(columns={
        "TPM": "Expression",
        "OriginalID": "Identifier",
        "Subpeptide": "pep_sequence",
        "Start": "pep_start",
        "End": "pep_end",
        "Index": "peptide_id",
        "virus": "Organism",
        "gene_symbol": "Gene",
    })

    outcols = ["peptide_id","Identifier","Organism","Gene","Expression",
               "pep_start","pep_end","pep_sequence","Peptide_Class"]
    output_file = os.path.join(args.outdir, "peptide_final_table.tsv")
    merged[outcols].to_csv(output_file, sep="\t", index=False)



# Step 2 - This function combines predictions and petide list
def CombinePrediction(args):
    re_bindingrank = re.compile("%RankBinding_([A-Z0-9_]+)")  # Class I MixMHC
    re_rank = re.compile(r"%Rank_([A-Z0-9_]+)")               # Class I PRIME or Class II MixMHC

    # Read peptide table
    peptides = pd.read_csv(args.table, sep="\t")

    # Read predictions
    preds_tables = {}
    if args.cI_preds:
        preds_tables["CI"] = pd.read_csv(args.cI_preds, sep="\t", comment="#")
    if args.cII_preds:
        preds_tables["CII"] = pd.read_csv(args.cII_preds, sep="\t", comment="#")

    def extract_rank_cols(df, regex):
        return [c for c in df.columns if regex.match(c)]

    out_entries = []

    for _, pepentry in peptides.iterrows():
        pep = pepentry.get("Peptide", pepentry.get("pep_sequence"))
        preds_for_pep = {ctype: df[df["Peptide"] == pep] for ctype, df in preds_tables.items() if pep in df["Peptide"].values}
        if not preds_for_pep:
            continue

        # Best core if available
        best_core = "NA"
        for df in preds_for_pep.values():
            if "Core_best" in df.columns and not df.loc[df["Peptide"] == pep, "Core_best"].empty:
                best_core = df.loc[df["Peptide"] == pep, "Core_best"].iloc[0]
                break

        rank_info = {}


        # MixMHC Class I
        scores_CI = {}
        if "CI" in preds_for_pep:
            df = preds_for_pep["CI"]
            for col in extract_rank_cols(df, re_bindingrank):
                allele = re_bindingrank.match(col).group(1)
                score = float(df.loc[df["Peptide"] == pep, col].iloc[0])
                scores_CI[allele] = score
            if scores_CI:
                best_rank = min(scores_CI.values())
                best_alleles = [a for a, s in scores_CI.items() if s == best_rank]
                additional = [a for a, s in scores_CI.items() if s <= args.BindingScore and a not in best_alleles]
                rank_info["CI"] = {"best_alleles": best_alleles, "additional": additional, "best_rank": best_rank}


        # MixMHC Class II
        scores_CII = {}
        if "CII" in preds_for_pep:
            df = preds_for_pep["CII"]
            for col in extract_rank_cols(df, re_rank):
                allele = re_rank.match(col).group(1)
                score = float(df.loc[df["Peptide"] == pep, col].iloc[0])
                scores_CII[allele] = score
            if scores_CII:
                best_rank = min(scores_CII.values())
                best_alleles = [a for a, s in scores_CII.items() if s == best_rank]
                additional = [a for a, s in scores_CII.items() if s <= args.BindingScore and a not in best_alleles]
                rank_info["CII"] = {"best_alleles": best_alleles, "additional": additional, "best_rank": best_rank}


        # Determine final peptide class
        classes = []
        if "CI" in rank_info:
            classes.append("HLA_I")
        if "CII" in rank_info:
            classes.append("HLA_II")
        final_class = "_".join(classes) if classes else "NA"

        # Skip if no MixMHC passes threshold
        pass_threshold = False
        if "CI" in rank_info and rank_info["CI"]["best_rank"] <= args.ScoreThreshold:
            pass_threshold = True
        if "CII" in rank_info and rank_info["CII"]["best_rank"] <= args.ScoreThreshold:
            pass_threshold = True
        if not pass_threshold:
            continue


        # Prepare output row
        new_entry = pepentry.copy()
        new_entry["best_core"] = best_core
        new_entry["Final_Peptide_Class"] = final_class

        # Combine best and additional alleles
        best_alleles_all = []
        additional_all = []

        if "CI" in rank_info:
            best_alleles_all += rank_info["CI"]["best_alleles"]
            additional_all += rank_info["CI"]["additional"]
        if "CII" in rank_info:
            best_alleles_all += rank_info["CII"]["best_alleles"]
            additional_all += rank_info["CII"]["additional"]

        new_entry["best_alleles"] = ",".join(best_alleles_all) if best_alleles_all else "NA"
        new_entry["other_significant_alleles"] = ",".join(additional_all) if additional_all else "NA"

        # Ranks
        rank_values = []
        if "CI" in rank_info:
            rank_values.append(rank_info["CI"]["best_rank"])
        if "CII" in rank_info:
            rank_values.append(rank_info["CII"]["best_rank"])
        new_entry["rank"] = min(rank_values) if rank_values else "NA"


        # PRIME columns (CI or CII)
        scores_prime = {}
        # CI PRIME
        if "CI" in preds_for_pep:
            df = preds_for_pep["CI"]
            for col in extract_rank_cols(df, re_rank):
                allele = re_rank.match(col).group(1)
                score = float(df.loc[df["Peptide"] == pep, col].iloc[0])
                scores_prime[allele] = score
            prime_class = "HLA_I"
        # CII fallback
        elif "CII" in preds_for_pep:
            df = preds_for_pep["CII"]
            for col in extract_rank_cols(df, re_rank):
                allele = re_rank.match(col).group(1)
                score = float(df.loc[df["Peptide"] == pep, col].iloc[0])
                scores_prime[allele] = score
            prime_class = "HLA_II"
        else:
            prime_class = "NA"

        if scores_prime:
            best_rank = min(scores_prime.values())
            best_alleles = [a for a, s in scores_prime.items() if s == best_rank]
            additional = [a for a, s in scores_prime.items() if s <= args.BindingScore and a not in best_alleles]
            new_entry["Final_Peptide_Class_PRIME"] = prime_class
            new_entry["best_alleles_PRIME"] = ",".join(best_alleles) if best_alleles else "NA"
            new_entry["other_significant_alleles_PRIME"] = ",".join(additional) if additional else "NA"
            new_entry["rank_PRIME"] = best_rank
        else:
            new_entry["Final_Peptide_Class_PRIME"] = "NA"
            new_entry["best_alleles_PRIME"] = "NA"
            new_entry["other_significant_alleles_PRIME"] = "NA"
            new_entry["rank_PRIME"] = "NA"

        out_entries.append(new_entry)

    if not out_entries:
        return

    outdf = pd.DataFrame(out_entries)
    outdf.to_csv(sys.stdout, sep="\t", index=False)


# Step 3 - this function add immunogenicity information AddAnnotation
def AddAnnotation(args):
    epitopes = pd.read_csv(args.vepitopes, sep="\t")
    table = pd.read_csv(args.table, sep="\t")

    hla_I = set(args.cI)
    hla_II = set(args.cII)

    viruses = epitopes.groupby("Peptide sequence")

    results = []
    for _, row in table.iterrows():
        pepseq = row["pep_sequence"]

        virus_name, gene_name, pep_immuno, nb_val, patient_imm_hla, pep_id = "NA","NA","NA","NA","NA", row["peptide_id"]

        if pepseq in viruses.groups:
            subset = viruses.get_group(pepseq)
            #all_tested = set.union(*[set(v.split(",")) for v in subset["HLA validated experimentaly"]])
            all_tested = set.union(*[set(str(v).split(",")) if pd.notna(v) else set() for v in subset["HLA validated experimentaly"]])
            all_immuno = set(subset["Peptide immunogenicity"])
            cI = hla_I.intersection(all_tested)
            cII = hla_II.intersection(all_tested)
            hla = cI.union(cII)

            pep_id = ";".join(subset["peptide_id"].astype(str))
            virus_name = "|".join(subset["Virus"])
            gene_name = "|".join(subset["gene name"])
            if "High" in all_immuno:
                pep_immuno = "High"
            elif "Medium" in all_immuno:
                pep_immuno = "Medium"
            elif "Low" in all_immuno:
                pep_immuno = "Low"
            else:
                pep_immuno = "unknown"
            nb_val = subset["Number of papers reporting peptide as immunogenic"].astype(int).sum()
            patient_imm_hla = ",".join(sorted(hla)) if hla else "NA"

        row["peptide_id"] = pep_id
        row["EPITOPES: Virus"] = virus_name
        row["EPITOPES: Gene name"] = gene_name
        row["EPITOPES: Peptide immunogenicity"] = pep_immuno
        row["EPITOPES: Number validations"] = nb_val
        row["EPITOPES: Patient immunogenic HLA"] = patient_imm_hla
        results.append(row)

    outdf = pd.DataFrame(results)
    outdf = outdf.map(lambda x: pd.NA if (isinstance(x, str) and x.strip() == "") else x)
    outdf.to_csv(sys.stdout, sep="\t", index=False, na_rep="NA")



# Step 4 - this function prioritizes the peptides PrioritizeViralPeptides
def PrioritizeViralPeptides(args):
    # Load table
    df = pd.read_csv(args.table, sep="\t")  # adjust sep if needed
    df = df.sort_values(by="rank")
    # Step 1: compute sorting features exactly like original
    immuno_map = {"High": 1, "Medium": 2, "Low": 3}
    df["pep_immuno_sorting"] = df["EPITOPES: Peptide immunogenicity"].map(immuno_map).fillna(4).astype(int)
    df["pep_nbphla_sorting"] = df["EPITOPES: Patient immunogenic HLA"].apply(lambda x: len(str(x).split(",")) if x != "NA" else 0)
    #df["pep_nbval_sorting"] = df["EPITOPES: Number validations"].replace("NA", -1).astype(int)
    df["pep_nbphla_sorting"] = df["EPITOPES: Patient immunogenic HLA"].apply(
        lambda x: 0 if pd.isna(x) or str(x).strip().upper() == "NA" else len([i for i in str(x).split(",") if i.strip()]))
    df["pep_nbval_sorting"] = df["EPITOPES: Number validations"].replace("NA", -1).fillna(-1).astype(int)
    df["Expression_sorting"] = df["Expression"].replace("NA", -1).astype(float)

    def compute_rank(row):
        vals = []
        for col in ["rank", "rank_PRIME"]:
            try:
                vals.append(float(row[col]))
            except:
                continue
        return mean(vals) if vals else 100
    df["rank_sorting"] = df.apply(compute_rank, axis=1)

    df["IsBinder"] = df.apply(
        lambda row: 1 if row["rank_sorting"] <= args.rank or row["EPITOPES: Peptide immunogenicity"] in ["High", "Medium", "Low"] else 2,
        axis=1
    )

    df["nb_best_alleles_sorting"] = df.apply(
        lambda row: len([al for al in set(str(row["best_alleles"]).split(",") + str(row["best_alleles_PRIME"]).split(",")) if al != "NA"]),
        axis=1
    )

    df["nb_additional_alleles_sorting"] = df.apply(
        lambda row: len([al for al in set(str(row["other_significant_alleles"]).split(",") + str(row["other_significant_alleles_PRIME"]).split(",")) if al != "NA"]),
        axis=1
    )
    # Replace all NaN with string "NA" before converting to dicts
    df = df.fillna("NA")

    # Step 2: convert to list of dicts and sort using Python's sorted() for exact tie-breaking
    entries = df.to_dict(orient="records")

    sorted_entries = sorted(
        entries,
        key=lambda e: (
            e["IsBinder"],
            -e["pep_nbphla_sorting"],
            e["pep_immuno_sorting"],
            -e["pep_nbval_sorting"],
            -e["Expression_sorting"],
            e["rank_sorting"],
            e["nb_best_alleles_sorting"],
            e["nb_additional_alleles_sorting"],
            e["Organism"],
            e["Gene"]
        )
    )

    # Step 3: assign ranks exactly like original
    ci_rank, cii_rank, overall_rank = 1, 1, 1
    print("\t".join(["Rank", "Rank_CI", "Rank_CII"] + list(df.columns)), file=sys.stdout)
    for entry in sorted_entries:
        cir, ciir = "NA", "NA"
        if entry["Final_Peptide_Class"] == "HLA_I_II":
            cir, ciir = ci_rank, cii_rank
            ci_rank += 1
            cii_rank += 1
        elif entry["Final_Peptide_Class"] == "HLA_I":
            cir = ci_rank
            ci_rank += 1
        elif entry["Final_Peptide_Class"] == "HLA_II":
            ciir = cii_rank
            cii_rank += 1

        outline = [str(overall_rank), str(cir), str(ciir)] + [str(entry[col]) for col in df.columns]
        print("\t".join(outline), file=sys.stdout)
        overall_rank += 1

# parse arguments
parser = argparse.ArgumentParser(prog='create_final_table_mixmhc_pd.py', formatter_class=argparse.RawTextHelpFormatter)
subparsers = parser.add_subparsers(help='sub-command help')

create_peptide_final_list_parser = subparsers.add_parser("create_peptide_final_list")
create_peptide_final_list_parser.add_argument("-peptide_table", required=True, type=str, help='Peptide table from get subpeptides')
create_peptide_final_list_parser.add_argument("-luttable", required=True, type=str, help='virus lookup table')
create_peptide_final_list_parser.add_argument("-virus_summary_final", required=True, type=str, help='viral summary - identified virus.txt')
create_peptide_final_list_parser.add_argument("-outdir", required=True, type=str, help='output directory')
create_peptide_final_list_parser.set_defaults(func=create_peptide_final_list)


CombinePrediction_parser = subparsers.add_parser("CombinePrediction")
CombinePrediction_parser.add_argument("-table", required=True, type=str, help='Table output from create_peptide_final_list')
CombinePrediction_parser.add_argument("-cI_preds", required=True, type=str, help='PRIME filtred prediction')
CombinePrediction_parser.add_argument("-cII_preds", required=True, type=str, help='MixMHC2pred filtred prediction')
CombinePrediction_parser.add_argument("-BindingScore", type=float, default=2)
CombinePrediction_parser.add_argument("-ScoreThreshold", type=float, default=100)
CombinePrediction_parser.set_defaults(func=CombinePrediction)

AddAnnotation_parser = subparsers.add_parser("AddAnnotation")
AddAnnotation_parser.add_argument("-table", required=True, type=str, help='Table output from CombinePrediction')
AddAnnotation_parser.add_argument("-vepitopes", required=True, type=str, help='List of Immunogenic Peptides from IEDB epitope')
AddAnnotation_parser.add_argument("-cI", nargs="+", default=[])
AddAnnotation_parser.add_argument("-cII", nargs="+", default=[])
AddAnnotation_parser.set_defaults(func=AddAnnotation)

PrioritizeViralPeptides_parser = subparsers.add_parser("PrioritizeViralPeptides")
PrioritizeViralPeptides_parser.add_argument("-table", required=True, type=str, help='Table output from AddAnnotation')
PrioritizeViralPeptides_parser.add_argument("-rank", type=float, default=1.0)
PrioritizeViralPeptides_parser.set_defaults(func=PrioritizeViralPeptides)

if __name__ == "__main__":
    args = parser.parse_args()
    if hasattr(args,"func"):
        logger = setup_logging(os.getcwd())
        try:
            args.func(args)
        except Exception:
            logger.exception("Unhandled exception during execution")
            sys.exit(2)
    else:
        parser.print_help()
