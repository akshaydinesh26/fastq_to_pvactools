process decompose_vcf {

    tag "$patient_id"
    publishDir "${params.outdir}/${patient_id}/final_annotation", mode: 'copy'
    
    input:
    tuple val(patient_id), path(vep_vcf)

    output:
    tuple val(patient_id), path("${patient_id}_decomposed.vcf"), emit: decompose_vcf
    path "decompose_vcf.log", emit: log_file

    // Error handling
    errorStrategy 'ignore'   // workflow continues even if this process fails
    maxRetries 1             // optional: do not retry failed tasks

    script:
    """
    /usr/local/bin/vt decompose -s "$vep_vcf" -o "${patient_id}_decomposed.vcf" > decompose_vcf.log 2>&1
    """
}
