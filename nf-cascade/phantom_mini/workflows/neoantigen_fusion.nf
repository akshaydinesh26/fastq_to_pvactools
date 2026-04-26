#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { run_pvacfuse }       from "${projectDir}/modules/pvacfuse.nf"

workflow neoantigen_fusion {
    take:
      hla_allele_input
      fusion_input

    main:
      pvacfuse_run = run_pvacfuse(hla_allele_input,fusion_input)

    emit:
      pvacseq_out    = pvacfuse_run.pvacfuse_dir
}
