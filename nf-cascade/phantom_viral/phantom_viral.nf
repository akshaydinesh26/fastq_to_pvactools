#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { star_human }                  from "${projectDir}/modules/human_star.nf"
include { star_virus }                  from "${projectDir}/modules/viral_star.nf"
include { vepitope_predictions }        from "${projectDir}/modules/process_virus.nf"
include { vepitope_predictions_mixmhc } from "${projectDir}/modules/process_virus_mixmhcpred.nf"
include { qc_fastp }                    from "${projectDir}/modules/qc_fastp.nf"
include { qc_fastqc }                   from "${projectDir}/modules/qc_fastqc.nf"

workflow {


  // Check if user asked for help
 if (params.help) {
        println """
        --------------------------------------------------------------
        Nextflow Viral Peptide Workflow
        --------------------------------------------------------------

        Required parameters:
          --sample_sheet      Path to the sample sheet (CSV or TSV)
          --genome            Genome reference to use ('GRCH37' or 'GRCH38'): default GRCH37


        Optional parameters:
          --outdir                              Output directory: default results in current folder
          --cpus                                Number of threads: default 8
          --use_mixmhc                          To run mixmhc as compared to pvacbind which is default for binding prediction
          --fastp_unqualified_percent_limit     fastp unqualified_percent_limit: default 0
          --fastp_qualified_quality_phred       fastp qualified_quality_phred
          --fastp_length_required               fastp length_required: default 0
          --max_memory                          maximum meory usage per process: default '100GB'

        Required parameters(set before run):
          --virus_proteome                      virus proteome created with helper scripts
          --immunogenic_peptide_info            Immunogenic peptide information from IEDB epitopes created with helper script
          --human_proteome                      human proteome
          --virus_reference_lutfile             Virus gene lookup table created with helper script from gtf file
          --human_reference_lutfile             Human gene lookup table prepared with helper script. The data obtained from Biomart
          --virus_reference_gtf                 virus reference gtf file
          --virus_genomeDir                     STAR index for viral fasta, create with sjdbOverhang 99
          --human_genomeDir                     STAR index for human fasta, create with sjdbOverhang 100
          --human_reference_gtf                 human gtf file
          --help for this message

        Example usage:
          nextflow run main.nf --sample_sheet samples.csv --genome GRCH38 --outdir results
        --------------------------------------------------------------
        """
        exit 0
    }

          // checks input file is correct
 if (!params.sample_sheet) {
        log.error " 'sample_sheet' parameter is required. Please provide it via --sample_sheet"
        exit 1
    }

//  Define input reference files
 human_genome_dir = params["${params.genome}"].human_genomeDir
 human_reference_gtf_file = params["${params.genome}"].human_reference_gtf
 virus_genome_dir =  params["${params.genome}"].virus_genomeDir
 virus_referennce_gtf_file = params["${params.genome}"].virus_reference_gtf
 virus_reference_lutfile_loc = params["${params.genome}"].virus_reference_lutfile
 human_reference_lutfile_loc = params["${params.genome}"].human_reference_lutfile
 human_proteome_file = params["${params.genome}"].human_proteome
 immunogenic_peptide_info_file = params["${params.genome}"].immunogenic_peptide_info
 immunogenic_peptide_info_file_mixmhc = params["${params.genome}"].immunogenic_peptide_info_mixmhc

    // Parameter checks - lut tables, reference fasta, star indices, proteome fasta
 if (!human_genome_dir) {
      error "Missing required parameter: --human_genomeDir"
    }

 if (!human_reference_gtf_file) {
      error "Missing required parameter: --human_reference_gtf"
    }

 if (!virus_genome_dir) {
      error "Missing required parameter: --virus_genomeDir"
    }

 if (!virus_referennce_gtf_file) {
      error "Missing required parameter: --virus_reference_gtf"
    }    

if (!virus_reference_lutfile_loc) {
      error "Missing required parameter: --virus_reference_lutfile"
    }

if (!human_reference_lutfile_loc) {
      error "Missing required parameter: --human_reference_lutfile"
    }

if (!params.virus_proteome) {
      error "Missing required parameter: --virus_proteome"
    }

if (!human_proteome_file) {
      error "Missing required parameter: --human_proteome"
    }

if (!immunogenic_peptide_info_file) {
      error "Missing required parameter: --immunogenic_peptide_info"
    }

if (!immunogenic_peptide_info_file_mixmhc) {
      error "Missing required parameter: --immunogenic_peptide_info_aa"
    }



// initial channel from tsv
Channel
    .fromPath(params.sample_sheet)
    .splitCsv(header: true, sep: '\t')
    .map { row ->
        tuple(
            row.patient_id,
            file(row.fastq_1),
            file(row.fastq_2),
            row.strandness,
            row.hla_alleles
        )
    }
    .set { input_data }

// prepare input for humans star alignment
//input_data
//    .map {patient_id, fastq_1, fastq_2, strandness, hla_alleles -> tuple(patient_id,[fastq_1,fastq_2])}
//    .set { human_star_input }

// prepare input for qc_check and filtering
input_data
    .map {patient_id, fastq_1, fastq_2, strandness, hla_alleles -> tuple(patient_id,[fastq_1,fastq_2])}
    .set { qc_input }

qc_fastqc(qc_input)
fastp_ch = qc_fastp(qc_input,params.fastp_length_required,params.fastp_qualified_quality_phred,params.fastp_unqualified_percent_limit)

// prepare input for humans star alignment
fastp_ch.hq_reads
    .map {patient_id, hq_fastqs -> tuple(patient_id,hq_fastqs)}
    .set { human_star_input }

//prepare input for vepitope_predictions
input_data
    .map {patient_id, fastq_1, fastq_2, strandness, hla_alleles -> tuple(patient_id,strandness,hla_alleles)}
    .set { vepitope_val_info }

// pass values to process - star alignment
star_alignment_ch = star_human(human_star_input,human_genome_dir,human_reference_gtf_file)

// prepare input file for star virus alignment 
star_alignment_ch.unmapped_reads
    .map {patient_id, unmapped_fastqs -> tuple(patient_id,unmapped_fastqs)}
    .set { viral_star_input }

// get human count to channel 
star_alignment_ch.human_count
    .map {patient_id, human_star_count -> tuple(patient_id,human_star_count)}
    .set { star_human_count }

// pass input file for star virus alignment
virus_alignment_ch = star_virus(viral_star_input,virus_genome_dir,virus_referennce_gtf_file)

virus_alignment_ch.virus_count
    .map {patient_id, virus_star_count -> tuple(patient_id,virus_star_count)}
    .join (star_human_count)
    .join (vepitope_val_info)
    .set { vepitope_predictions_input }

if (params.use_mixmhc){

    vepitope_predictions_final = vepitope_predictions_mixmhc(vepitope_predictions_input,
    virus_reference_lutfile_loc,human_reference_lutfile_loc,
    params.virus_proteome,human_proteome_file,
    immunogenic_peptide_info_file_mixmhc,
    params.genome
    )

} else {

    vepitope_predictions_final = vepitope_predictions(vepitope_predictions_input,
    virus_reference_lutfile_loc,human_reference_lutfile_loc,
    params.virus_proteome,human_proteome_file,
    immunogenic_peptide_info_file,
    params.genome
    )
}


}

workflow.onComplete {
    println """
    ===========================
        Pipeline Completed
    ---------------------------
    Completed at : ${workflow.complete}
    Duration     : ${workflow.duration}
    Success      : ${workflow.success}
    Work Dir     : ${workflow.workDir}
    Output Dir   : ${params.outdir}
    Exit Status  : ${workflow.exitStatus}
    ===========================
    """.stripIndent()
}

