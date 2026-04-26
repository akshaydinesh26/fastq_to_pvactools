process generate_site_list {
    tag "${patient_id}"
    publishDir "${params.outdir}/${patient_id}/final_annotation", mode: 'copy'
    
    input:
    tuple val(patient_id), path(decompose_vcf), path(dna_cram), path(rna_bam)
    path(dna_fasta)

    output:
    tuple val(patient_id), path("site_list_dir_dna"), emit: dna_site_list
    tuple val(patient_id), path("site_list_dir_rna"), emit: rna_site_list
    path("generate_site_list.log"), emit: log_file

    script:
    """
    {
    
    # create output folders
    mkdir -p site_list_dir_dna
    mkdir -p site_list_dir_rna

    # index bam and cram files
    samtools index ${dna_cram}
    
    # get tumor sample ID
    tumor_sample=\$( grep "^#CHROM" ${decompose_vcf} | awk '{print \$NF}' )

    # generate site list
    /usr/bin/bam_readcount_helper.py ${decompose_vcf} \${tumor_sample} ${dna_fasta} ${dna_cram} NOPREFIX site_list_dir_dna
    
    #for rnaseq
    if [ -s "${rna_bam}" ]; then
    samtools index ${rna_bam}
    /usr/bin/bam_readcount_helper.py ${decompose_vcf} \${tumor_sample} ${dna_fasta} ${rna_bam} NOPREFIX site_list_dir_rna
    fi


    } > generate_site_list.log 2>&1
    
    
    """
}
