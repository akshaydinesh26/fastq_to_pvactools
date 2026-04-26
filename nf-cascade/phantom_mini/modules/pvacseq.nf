process run_pvacseq {

    tag "$patient_id"
    publishDir "${params.outdir}/${patient_id}", mode: 'copy'

    input:
    // tuple val(patient_id), val(hla)
    // tuple val(patient_id), path(expr_transcript_vcf)
    tuple val(patient_id), val(hla), path(expr_transcript_vcf)

    output:
    tuple val(patient_id), path("pvacseq"), emit: pvacseq_dir
    path "pvacseq.log", emit: log_file

    // Error handling
    errorStrategy 'ignore'   // workflow continues even if this process fails
    maxRetries 1             // optional: do not retry failed tasks

    script:
    """
    # get tumor and normal name
    tumor_sample=\$( grep "^#CHROM" ${expr_transcript_vcf} | awk '{print \$NF}' )
    normal_sample=\$( grep "^#CHROM" ${expr_transcript_vcf} | awk '{print \$(NF-1)}' )
    
    # create output folder
    mkdir -p pvacseq

    # run pvacseq
    pvacseq run \\
        ${expr_transcript_vcf} \\
        \${tumor_sample} \\
        "${hla}" \\
        all \\
        pvacseq \\
        --n-threads ${params.cpus} \\
        --normal-sample-name \${normal_sample} \\
        -e1 8,9,10,11 \\
        -e2 12,13,14,15,16,17,18 \\
        --pass-only \\
        --iedb-install-directory /opt/iedb/ \\
        > pvacseq.log 2>&1
    """
}
