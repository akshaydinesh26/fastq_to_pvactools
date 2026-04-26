process generate_junction {
    tag "${patient_id}"
    publishDir "${params.outdir}/${patient_id}/final_annotation", mode: 'copy'

    input:
    tuple val(patient_id), path(rna_bam), path(expr_transcript_vcf)
    path(rna_fasta)
    path(rna_gtf)

    output:
    tuple val(patient_id), path("${patient_id}_splice_annotation.tsv"), emit: junction_list
    path "generate_junction.log", emit: log_file

    // Error handling
    errorStrategy 'ignore'   // workflow continues even if this process fails
    maxRetries 1             // optional: do not retry failed tasks

    script:
    """
    {
    # index bam file
    samtools index ${rna_bam}
    
    # generate junction list
    regtools cis-splice-effects identify \\
        -o ${patient_id}_splice_annotation.tsv \\
        -s intron-motif \\
        -C \\
        ${expr_transcript_vcf} \\
        ${rna_bam} \\
        ${rna_fasta} \\
        ${rna_gtf}
    } > generate_junction.log 2>&1
    
    """
}
