process qc_fastqc {

    tag "$patient_id"
    publishDir "${params.outdir}/${patient_id}/quality_control", mode: 'copy'

    input:
    tuple val(patient_id), path(fastqs)

    output:
    tuple val(patient_id), path("fastqc")
    script:
    """
    mkdir -p fastqc

    fastqc -t "${task.cpus}" -o fastqc "${fastqs[0]}" "${fastqs[1]}"
    """
}
