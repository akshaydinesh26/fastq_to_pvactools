#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// run as nextflow run <path to phantom.nf> --sample_sheet <input_sample_sheet>
//--cpus <value for cpu> --max_memory <memory as '200 GB' default> --outdir <output folder> 


include { annotate_pvacseq }                   from "${projectDir}/workflows/annotate_pvacseq.nf"
include { annotate }                           from "${projectDir}/workflows/annotate.nf"
include { neoantigen_fusion }                  from "${projectDir}/workflows/neoantigen_fusion.nf"
include { neoantigen_prediction_pvacseq }      from "${projectDir}/workflows/neoantigen_prediction_pvacseq.nf"
include { neoantigen_prediction }              from "${projectDir}/workflows/neoantigen_prediction.nf"
include { create_empty_files }                 from "${projectDir}/modules/placeholders.nf"

workflow {

      // checks input file is correct
 if (!params.sample_sheet) {
        log.error " 'sample_sheet' parameter is required. Please provide it via --sample_sheet"
        exit 1
    }

     // Parameter checks - dna, gtf, outdir, max_memory
 if (!params.dna_fasta) {
      error "Missing required parameter: --dna_fasta"
    }
    
if (!params.dna_gtf) {
      error "Missing required parameter: --dna_gtf"
    }
    // Optional parameters: warn if unset
   
if (!params.outdir) {
      log.warn "Output directory (--outdir) not set. Using default working directory."
    }

if (!params.max_memory) {
      log.warn "Maximum memory (--max_memory) not set. Processes may use default memory limits."
    }

// Print input provided and directories used    
 println "=== Pipeline Starting ==="
    println "Launch Directory : ${workflow.launchDir}"
    println "Working Directory: ${workflow.workDir}"
    println "Output Directory : ${params.outdir}"
    println "-------------------------"
    println "Sample Sheet     : ${params.sample_sheet}"
    println "Reference FASTA  : ${params.dna_fasta}"
    println "GTF Annotation   : ${params.dna_gtf}"
    println "CPUs Requested   : ${params.cpus}"
    println "Max Memory       : ${params.max_memory}"
    println "========================="

// make input references to channel
dna_fasta_ch   = Channel.value(file(params.dna_fasta))
dna_gtf_ch     = Channel.value(file(params.dna_gtf))

//empty_files = create_empty_files()

// initial channel from tsv
    Channel
    .fromPath(params.sample_sheet)
    .splitCsv(header: true, sep: '\t')
    .map { row ->
        tuple(
            row.PATIENT,
            row.MUTECT_VCF ? file(row.MUTECT_VCF) : file('empty.vcf'),
            row.DNA_CRAM ? file(row.DNA_CRAM) : file('empty.cram'),
            row.RNA_BAM ? file(row.RNA_BAM) : file('empty.bam'),
            row.RSEM_GENE ? file(row.RSEM_GENE) : file('empty.rsem_gene'),
            row.RSEM_TX ? file(row.RSEM_TX) : file('empty.rsem_tx'),
            row.HLA_ALLELES ? row.HLA_ALLELES : val('NA'),
            row.RNA_FUSION ? file(row.RNA_FUSION) : file('empty_fusion')
       )
    }
    .set { annotate_input }

    main:
    
    if (params.wes_only){
 
        annotate_run = annotate_pvacseq(annotate_input,dna_fasta_ch,dna_gtf_ch)
        neoantigen_prediction_complete = neoantigen_prediction_pvacseq(
            pvacseq_input_combined = annotate_run.pvacseq_input_combined_ch,
            )
    } else {
    
    annotate_run = annotate(annotate_input,dna_fasta_ch,dna_gtf_ch)
    neoantigen_prediction_complete = neoantigen_prediction(
        pvacseq_input_combined = annotate_run.pvacseq_input_combined_ch,
        pvacsplice_input_combined = annotate_run.pvacsplice_input_combined_ch,
        dna_fasta_ch = dna_fasta_ch,
        dna_gtf_ch = dna_gtf_ch
    )

    // create input for fusion
    annotate_input
    .map { patient_id, col2, col3, col4, col5, col6, col7, rna_fusion -> tuple(patient_id,rna_fusion) }
    .set { fusion_input }

    // Prepare input channel for hla for fusion
    annotate_input
        .map { patient_id, col2, col3, col4, col5, col6, hla_alleles, col8 -> tuple(patient_id, hla_alleles) }
        .set { hla_allele_input }
    }
  
   if (!params.pvacseq_only && !params.wes_only ){

   neoantigen_fusion_run = neoantigen_fusion(hla_allele_input,fusion_input)

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
