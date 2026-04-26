process run_pvacfuse {

    tag "$patient_id"
    publishDir "${params.outdir}/${patient_id}", mode: 'copy'

    input:
    tuple val(patient_id), val(hla)
    tuple val(patient_id), path(rna_fusion)

    output:
    tuple val(patient_id), path("pvacfuse"), emit: pvacfuse_dir
    path "pvacfuse.log", emit: log_file


    // Error handling
    errorStrategy 'ignore'   // workflow continues even if this process fails
    maxRetries 1             // optional: do not retry failed tasks

    script:
    """
    # create output folder
    mkdir -p pvacfuse

    # pvacfuse run
    pvacfuse run \\
        ${rna_fusion} \\
        ${patient_id} \\
        ${hla} \\
        all \\
        pvacfuse \\
        --n-threads ${params.cpus} \\
        -e1 8,9,10,11 \\
        -e2 12,13,14,15,16,17,18 \\
        --iedb-install-directory /opt/iedb/ \\
        > pvacfuse.log 2>&1
    """
}
