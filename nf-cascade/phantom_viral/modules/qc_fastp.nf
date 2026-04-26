process qc_fastp {

    tag "$patient_id"
    publishDir "${params.outdir}/${patient_id}/quality_control/fastp", mode: 'copy', pattern: "${patient_id}_fastp.html"

    input:
    tuple val(patient_id), path(fastqs)
    val(length_required)
    val(qualified_quality_phred)
    val(unqualified_percent_limit)

    output:
    tuple val(patient_id), path({["${patient_id}_R1_filtered.fastq.gz", "${patient_id}_R2_filtered.fastq.gz"]}), emit: hq_reads
    tuple val(patient_id), path("${patient_id}_fastp.html"), path("${patient_id}_fastp.json")
    script:
    """

    fastp \
     -i "${fastqs[0]}" \
     -I "${fastqs[1]}" \
     -o "${patient_id}_R1_filtered.fastq.gz" \
     -O "${patient_id}_R2_filtered.fastq.gz" \
     --html "${patient_id}_fastp.html" \
     --json "${patient_id}_fastp.json" \
     --thread "${task.cpus}" \
     --detect_adapter_for_pe \
     --disable_quality_filtering \
     --disable_adapter_trimming \
     --disable_length_filtering \
     --disable_trim_poly_g \
     --length_required "${length_required}" \
     --qualified_quality_phred "${qualified_quality_phred}" \
     --unqualified_percent_limit "${unqualified_percent_limit}" \
    """
}
