process star_virus {

    tag "$patient_id"
    publishDir "${params.outdir}/${patient_id}/viral_alignment", mode: 'copy'
    
    input:
    tuple val(patient_id), path(fastqs)
    path(index)
    path(gtf)

    output:
    tuple val(patient_id), path("${patient_id}_virus_ReadsPerGene.out.tab"), emit: virus_count
    tuple val(patient_id), path("${patient_id}_virus_Aligned.out.bam"), emit: alignment_bam
    script:
    """
    # Check if both unmapped FASTQs are empty
    if [ ! -s "${fastqs[0]}" ] && [ ! -s "${fastqs[1]}" ]; then
        echo "No unmapped reads found for $patient_id. Skipping downstream steps."
        touch "${patient_id}_virus_ReadsPerGene.out.tab"
        touch "${patient_id}__virus_Aligned.out.bam"
        exit 0
    fi

    STAR --runMode alignReads \
     --readFilesCommand zcat \
     --genomeDir $index \
     --sjdbGTFfile $gtf \
     --twopassMode Basic \
     --outSAMtype BAM Unsorted \
     --outReadsUnmapped Fastx \
     --quantMode GeneCounts \
     --sjdbOverhang 99 \
     --readFilesIn "${fastqs[0]}" "${fastqs[1]}" \
     --outFileNamePrefix "${patient_id}_virus_" \
     --runThreadN ${task.cpus} \
    """
}
