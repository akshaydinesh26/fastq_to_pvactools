#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { run_pvacseq }        from "${projectDir}/modules/pvacseq.nf"
include { run_pvacsplice }     from "${projectDir}/modules/pvacsplice.nf"


workflow neoantigen_prediction {
    take:
      pvacseq_input_combined
      pvacsplice_input_combined
      dna_fasta_ch
      dna_gtf_ch

    main:
      pvacseq_run    = run_pvacseq(pvacseq_input_combined)
      pvacsplice_run = run_pvacsplice(pvacsplice_input_combined,dna_fasta_ch,dna_gtf_ch)
      
    emit:
      pvacseq_out    = pvacseq_run.pvacseq_dir ?: Channel.empty()
      pvacsplice_out = pvacsplice_run.pvacsplice_dir ?: Channel.empty()
}
