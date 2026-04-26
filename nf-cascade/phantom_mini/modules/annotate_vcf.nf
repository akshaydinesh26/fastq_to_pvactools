process vcf_annotate_all {
    tag "$patient_id"
    publishDir "${params.outdir}/${patient_id}/final_annotation", mode: 'copy'

    input:
    tuple val(patient_id), path(decomposed_vcf), path(rsem_gene), path(rsem_tx), path(site_list_dir_dna), path(site_list_dir_rna)


    output:
    tuple val(patient_id), path("${patient_id}_vcf_expr_transcript.vcf"), emit: vcf_annotated, optional:true
    path "${patient_id}_vcf_annotated_DNA_SNV.vcf"
    tuple val(patient_id), path("${patient_id}_vcf_annotated_DNA_SNV_INDEL.vcf"), emit: vcf_for_wes_only
    path "${patient_id}_vcf_annotated_DNA_RNA_SNV.vcf", optional: true
    path "${patient_id}_vcf_annotated_DNA_RNA_SNV_INDEL.vcf", optional: true
    path "${patient_id}_vcf_expr_gene.vcf", optional: true
    path "annotate_vcf.log", emit: log_file

    // Error handling
    errorStrategy 'ignore'   // workflow continues even if this process fails
    maxRetries 1             // optional: do not retry failed tasks

    script:
    """
    {
    tumor_sample=\$( grep "^#CHROM" ${decomposed_vcf} | awk '{print \$NF}' )

    SNP_SITES_SNV_DNA="${site_list_dir_dna}/\${tumor_sample}_bam_readcount_snv.tsv"
    SNP_SITES_INDEL_DNA="${site_list_dir_dna}/\${tumor_sample}_bam_readcount_indel.tsv"
    SNP_SITES_SNV_RNA="${site_list_dir_rna}/\${tumor_sample}_bam_readcount_snv.tsv"
    SNP_SITES_INDEL_RNA="${site_list_dir_rna}/\${tumor_sample}_bam_readcount_indel.tsv"


    # Step 1: DNA SNV
    vcf-readcount-annotator \\
        ${decomposed_vcf} \${SNP_SITES_SNV_DNA} DNA \\
        -s \${tumor_sample} -t snv \\
        -o ${patient_id}_vcf_annotated_DNA_SNV.vcf

    # Step 2: DNA INDEL
    vcf-readcount-annotator \\
        ${patient_id}_vcf_annotated_DNA_SNV.vcf \${SNP_SITES_INDEL_DNA} DNA \\
        -s \${tumor_sample} -t indel \\
        -o ${patient_id}_vcf_annotated_DNA_SNV_INDEL.vcf

    if [ -s "${rsem_gene}" ] && [ -s "${rsem_tx}" ]; then

    # Step 3: RNA SNV
    vcf-readcount-annotator \\
        ${patient_id}_vcf_annotated_DNA_SNV_INDEL.vcf \${SNP_SITES_SNV_RNA} RNA \\
        -s \${tumor_sample} -t snv \\
        -o ${patient_id}_vcf_annotated_DNA_RNA_SNV.vcf

    # Step 4: RNA INDEL
    vcf-readcount-annotator \\
        ${patient_id}_vcf_annotated_DNA_RNA_SNV.vcf \${SNP_SITES_INDEL_RNA} RNA \\
        -s \${tumor_sample} -t indel \\
        -o ${patient_id}_vcf_annotated_DNA_RNA_SNV_INDEL.vcf

    # Step 5: Gene expression
    vcf-expression-annotator \\
        ${patient_id}_vcf_annotated_DNA_RNA_SNV_INDEL.vcf ${rsem_gene} custom gene \\
        -i gene_id -e TPM -s \${tumor_sample} \\
        -o ${patient_id}_vcf_expr_gene.vcf \\
        --ignore-ensembl-id-version

    # Step 6: Transcript expression
    vcf-expression-annotator \\
        ${patient_id}_vcf_expr_gene.vcf ${rsem_tx} custom transcript \\
        -i transcript_id -e TPM -s \${tumor_sample} \\
        -o ${patient_id}_vcf_expr_transcript.vcf \\
        --ignore-ensembl-id-version
    
    fi
    } > annotate_vcf.log 2>&1

    """
}
