#!/usr/bin/env python3
import csv
import subprocess
import argparse
import sys
import shutil
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

def get_launch_dir() -> Path:
    """Directory from which the script is launched."""
    return Path.cwd()


def get_project_dir() -> Path:
    """Directory where this script resides."""
    return Path(__file__).resolve().parent


# Permanent Nextflow script location (relative to project dir)
NEXTFLOW_SCRIPT = get_project_dir() / "nf-cascade" / "main.nf"

# check if nextflow worlfow completed for patient
def check_nextflow_log(patient_workdir):
    """Check Nextflow's .nextflow.log to see if run completed."""
    nf_log = Path(patient_workdir) / ".nextflow.log"
    if not nf_log.exists():
        return False
    try:
        with open(nf_log) as f:
            lines = [l.strip() for l in f if l.strip()]
        return bool(lines and "Execution complete" in lines[-1]) # change this tag as per output of nextflow log
    except Exception:
        return False


def get_failed_or_unrun_patients(rows, launch_dir):
    """Detect failed/unrun patients by checking wrapper logs and Nextflow logs."""
    failed_or_unrun = []
    for row in rows:
        patient_id = row["patient"]
        log_file = Path(launch_dir) / "tmp" / patient_id / "nextflow.log"

        if not log_file.exists():
            failed_or_unrun.append(row)
            continue

        with open(log_file) as f:
            lines = [l.strip() for l in f if l.strip()]

        if not lines:
            failed_or_unrun.append(row)
            continue

        if any("Neoantigen analysis finished" in line for line in lines):
            continue
        else:
            failed_or_unrun.append(row)

    return failed_or_unrun

def run_phantom_for_row(row, nextflow_script_path, final_output_dir, resume,launch_dir,custom_workflow,wes_only,pvacseq_only):
    patient_id = row["patient"]

    # temp_base = Path(tempfile.mkdtemp(dir=get_launch_dir()))
    temp_base = Path(launch_dir) / "tmp"
    patient_temp_dir = temp_base / patient_id
    patient_workdir = patient_temp_dir / "work"
    patient_tmp_outdir = patient_temp_dir / "output"

    patient_csv = patient_temp_dir / "patient.csv"
    temp_base.mkdir(parents=True, exist_ok=True)
    patient_temp_dir.mkdir(parents=True, exist_ok=True)

    # Write single-patient CSV
    with open(patient_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=row.keys(),delimiter="\t")
        writer.writeheader()
        writer.writerow(row)

    patient_logfile = final_output_dir / f"{patient_id}.log"

    cmd = [
        "nextflow", "run", str(nextflow_script_path),
        "--input_csv", str(patient_csv),
        "-work-dir", str(patient_workdir),
        "--output_final", str(patient_tmp_outdir),
        "--custom_workflow", str(custom_workflow)
    ]
    if resume:
        cmd.append("-resume")
    if wes_only:
        cmd.append("--wes_only")
    if pvacseq_only:
        cmd.append("--pvacseq_only")

    print(f"Starting patient {patient_id} ...")
    # === LIVE OUTPUT + LOG FILE ===
    # Keep this block if you want to see Nextflow output in the terminal and save it to logs.
    # Replace with simple Popen(..., stdout=log, stderr=log) if you want file-only logging.
    with open(patient_logfile, "a") as log:
        process = subprocess.Popen(cmd, 
                                   stdout=subprocess.PIPE, 
                                   stderr=subprocess.STDOUT, 
                                   text=True,
                                   cwd=patient_temp_dir)
        for line in process.stdout:
            sys.stdout.write(line)      # print live to terminal
            log.write(line)              # save to log file
        retcode = process.wait()
    # === END LIVE OUTPUT BLOCK ===

    if retcode == 0:
        print(f"Patient {patient_id} completed successfully")
        # Move temporary output to final location
        final_patient_dir = final_output_dir / patient_id
        if final_patient_dir.exists():
            shutil.rmtree(final_patient_dir)
        shutil.move(str(patient_tmp_outdir), str(final_patient_dir))
        shutil.move(str(patient_csv), str(final_patient_dir / "patient.csv"))
        #shutil.rmtree(patient_temp_dir, ignore_errors=True)
        return patient_id, True
    else:
        print(f"Patient {patient_id} failed with exit code {retcode}")
        return patient_id, False



def main():
    parser = argparse.ArgumentParser(description="Smart resume wrapper for running Sarek on multiple CSV rows.")
    parser.add_argument("--input_tsv", required=True, help="TSV file with one row per patient/sample")
    parser.add_argument("--output_dir", default="nf_output", help="Final output directory (relative to launch dir)")
    parser.add_argument("--max_parallel", type=int, default=1, help="Max patients to run in parallel")
    parser.add_argument("--stop-on-error", action="store_true", help="Stop all runs if one fails")
    parser.add_argument("--force-all", action="store_true", help="Run all patients regardless of status")
    parser.add_argument("--patients", help="Comma-separated list of patient IDs to run, overrides auto-detect")
    parser.add_argument("-resume", action="store_true", help="Resume previous runs using Nextflow's -resume")
    parser.add_argument("-wes_only", action="store_true", help="Only run wes based in neoantigen prediction")
    parser.add_argument("-pvacseq_only", action="store_true", help="Only run pvacseq in neoantigen prediction")
    parser.add_argument("-custom_workflow",type=str, default="none", help=("custom workflow except the default(which is identified based on if the "
                                                                       "hla alleles are provided or not). The allowed values are "
                                                                       "sarek|rnaseq|rnafusion|virus|sarek,rnaseq|sarek,rnaseq,virus|"
                                                                       "sarek,rnafusion|sarek,rnafusion,virus|rnaseq,rnafusion|"
                                                                       "rnaseq,rnaseq,virus|sarek,rnaseq,rnafusion|sarek,rnaseq,rnafusion,virus|"
                                                                       "sarek,rnaseq,rnafusion,hlatyping,tronflow|"
                                                                       "sarek,rnaseq,rnafusion,hlatyping,tronflow,virus"))
    args = parser.parse_args()

    launch_dir = get_launch_dir()
    final_output_dir = Path(args.output_dir)
    #log_dir = final_output_dir / "logs"

    final_output_dir.mkdir(parents=True, exist_ok=True)
    #log_dir.mkdir(parents=True, exist_ok=True)

    nextflow_script_path = NEXTFLOW_SCRIPT.resolve()

    with open(args.input_tsv) as f:
        reader = csv.DictReader(f,delimiter="\t")
        rows = [row for row in reader]

    if args.patients:
        patient_ids = set(p.strip() for p in args.patients.split(","))
        patients_to_run = [r for r in rows if r["patient"] in patient_ids]
    elif args.force_all:
        patients_to_run = rows
    else:
        patients_to_run = get_failed_or_unrun_patients(rows, launch_dir)

    if not patients_to_run:
        print("All selected patients successfully processed. Nothing to run.")
        sys.exit(0)

    print(f"▶ Running {len(patients_to_run)} patient(s): {[r['patient'] for r in patients_to_run]}")

    results = []
    with ThreadPoolExecutor(max_workers=args.max_parallel) as executor:
        future_to_patient = {
            executor.submit(run_phantom_for_row, row, nextflow_script_path, final_output_dir, args.resume, launch_dir,args.custom_workflow,args.wes_only,args.pvacseq_only): row["patient"]
            for row in patients_to_run
        }

        for future in as_completed(future_to_patient):
            patient_id = future_to_patient[future]
            try:
                pid, success = future.result()
                results.append((pid, success))
                if not success and args.stop_on_error:
                    print(f"Stopping early due to failure in {pid}")
                    executor.shutdown(wait=False, cancel_futures=True)
                    sys.exit(1)
            except Exception as e:
                print(f"Exception in patient {patient_id}: {e}")
                if args.stop_on_error:
                    sys.exit(1)

    print("\n=== Summary ===")
    for pid, success in results:
        print(f"{pid}: {'SUCCESS' if success else 'SUCCESS'}")


if __name__ == "__main__":
    main()
