# PHANTOM  
**Personalized Human ANtigen Target Output Module**

## Overview

**Version:** 01  

PHANTOM is a comprehensive pipeline that integrates:

- Somatic variant calling from tumor-normal WES data  
- RNA-seq expression analysis  
- Gene fusion detection  
- HLA typing  
- Neoantigen prediction  

It combines multiple nf-core workflows along with custom modules to generate prioritized neoantigen candidates for personalized cancer research.

---

## Workflows Included

### nf-core Pipelines

#### 1. Sarek (v3.5.1)
- Somatic mutation prediction from tumor-normal WES
- Alignment: bwa-mem  
- Variant calling: GATK  

#### 2. RNA-seq (v3.19.0)
- Gene expression quantification from tumor RNA-seq  
- Alignment: STAR  
- Quantification: RSEM  

#### 3. RNAfusion (v3.0.2)
- Fusion detection from RNA-seq  
- Alignment: STAR  
- Fusion prediction: Arriba  

#### 4. HLA Typing (optional)
- Predicts MHC Class I alleles from normal WES  

---

### Custom Workflows

#### Tronflow (optional)
- Predicts MHC Class II alleles from normal WES  

#### Phantommini
- Core neoantigen prediction module  
- Integrates outputs from all preprocessing workflows  
- Uses:
  - pvacseq (SNVs/indels)
  - pvacfuse (fusion-derived neoantigens)
  - pvacsplice (splice variants)

---

## Key Features

- Modular Nextflow DSL2-based architecture
- Each workflow runs independently with resume support
- Parallel execution of preprocessing steps
- Docker-based execution (except Tronflow)
- Python wrapper for orchestration and parallelization
- Scalable for multi-patient analysis

---

## Prerequisites

- Linux OS  
- Nextflow ≥ 24.10.6  
- Docker (non-root access configured)  
- Python ≥ 3.12  
- Human genome build: UCSC hg38 / GRCh38  

---

## Reference Data Setup

### Sarek References
Download from Illumina iGenomes:
https://sapac.support.illumina.com/sequencing/sequencing_software/igenome.html

Required:
- Reference genome
- GATK bundle (germline + PON)

Directory structure:
<igenome_base>/Homo_sapiens/GATK/GRCh38/

---

### RNA-seq References

Required:
- Genome FASTA (hg38)
- GTF annotation (GENCODE v47)
- STAR index

STAR Index Generation:

```
STAR \
--runThreadN 8 \
--runMode genomeGenerate \
--genomeDir <path/to/index> \
--genomeFastaFiles genome.fa \
--sjdbGTFfile annotation.gtf \
--sjdbOverhang 100
```

Note:
- Requires ~100GB RAM  
- Minimum 20 threads recommended  

---

### RNAfusion References

```
nextflow run nf-core/rnafusion \
--build_references --all \
--tool arriba \
--genomes_base <PATH> \
--outdir <PATH>
```

---

### Tronflow Setup

GitHub: https://github.com/TRON-Bioinformatics/tronflow  

Requires HLA-HD installation:
https://w3.genome.med.kyoto-u.ac.jp/HLA-HD/

Update path in:
nf-cascade/tronflow-hla-hd/nextflow.config

---

## Phantommini Inputs

- Tumor CRAM (Sarek)
- VCF (Sarek)
- RNA BAM (RNA-seq)
- Fusion results (RNAfusion)
- HLA alleles

Additional Requirements:
- Reference FASTA (hg38)
- GTF (GENCODE v47)

---

## Workflow Structure

Execution Strategy:

1. Parallel preprocessing
   - Sarek  
   - RNA-seq  
   - RNAfusion  
   - HLA typing (optional)  
   - Tronflow (optional)

2. Neoantigen prediction
   - Integrated via Phantommini  

---

## Setup Instructions

1. Extract or clone repository:
   phantom_final/

2. Navigate to:
   nf-cascade/params/

3. Configure parameter files:

Required Edits:

sarek_params.yaml
- Set igenome_base = <absolute_path>

rnaseq_params.yaml
- Set star_index, fasta, gtf

rnafusion_params.yaml
- Set genome_base

neoantigen_params.yaml
- Set dna_fasta, dna_gtf

---

## Execution

Run using Python wrapper:

```
python <path_to>/phantom.py --input_csv input.csv
```

---

## Command Options

| Option | Description |
|------|-------------|
| --input_csv | Input CSV (one row per patient) |
| --output_dir | Output directory (default: nf_output) |
| --max_parallel | Number of patients to run in parallel |
| --stop-on-error | Stop all runs if one fails |
| --force-all | Re-run all patients |
| --patients | Specific patient IDs (comma-separated) |
| -resume | Resume previous runs |

---

## Notes

- HLA typing workflows run only if HLA input is missing  
- Reference caching significantly improves runtime  
- Designed for reproducibility and HPC scalability  
