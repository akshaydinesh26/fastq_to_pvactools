#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { run_pvacseq }        from "${projectDir}/modules/pvacseq.nf"


workflow neoantigen_prediction_pvacseq {
    take:
      pvacseq_input_combined

    main:
      pvacseq_run    = run_pvacseq(pvacseq_input_combined)

    emit:
      pvacseq_out    = pvacseq_run.pvacseq_dir ?: Channel.empty()
}
