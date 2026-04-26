process star_human {

    tag "$patient_id"
    publishDir "${params.outdir}/${patient_id}/human_alignment", mode: 'copy'
    
    input:
    tuple val(patient_id), path(fastqs)
    path(index)
    path(gtf)

    output:
    tuple val(patient_id), path("${patient_id}_human_ReadsPerGene.out.tab"), emit: human_count
    //tuple val(patient_id), path(["${patient_id}_human_Unmapped_R1.fastq.gz", "${patient_id}_human_Unmapped_R2.fastq.gz"]), emit: unmapped_reads
    tuple val(patient_id), path({ ["${patient_id}_human_Unmapped_R1.fastq.gz", "${patient_id}_human_Unmapped_R2.fastq.gz"] }), emit: unmapped_reads
    tuple val(patient_id), path("${patient_id}_human_Aligned.out.bam"), emit: alignment_bam
    script:
    """
    STAR --runMode alignReads \
     --readFilesCommand zcat \
     --genomeDir $index \
     --sjdbGTFfile $gtf \
     --twopassMode Basic \
     --outSAMtype BAM Unsorted \
     --outReadsUnmapped Fastx \
     --quantMode GeneCounts \
     --sjdbOverhang 100 \
     --readFilesIn "${fastqs[0]}" "${fastqs[1]}" \
     --outFileNamePrefix "${patient_id}_human_" \
     --runThreadN ${task.cpus} \

    if [ ! -s "${patient_id}_human_Unmapped.out.mate1" ] || [ ! -s "${patient_id}_human_Unmapped.out.mate2" ]; then
        echo "Error: There is no unmapped reads for $patient_id!" >&2
    fi
    
    gzip -c "${patient_id}_human_Unmapped.out.mate1" > "${patient_id}_human_Unmapped_R1.fastq.gz"
    gzip -c "${patient_id}_human_Unmapped.out.mate2" > "${patient_id}_human_Unmapped_R2.fastq.gz"
    """
}
