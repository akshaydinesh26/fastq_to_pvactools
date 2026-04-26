process run_pvacsplice {

    tag "$patient_id"
    publishDir "${params.outdir}/${patient_id}", mode: 'copy'

    input:
    tuple val(patient_id), val(hla), path(expr_transcript_vcf), path(junctions)
    path(dna_fasta)
    path(dna_gtf)

    output:
    path "pvacsplice", emit: pvacsplice_dir
    path "pvacsplice.log", emit: log_file

    // Error handling
    errorStrategy 'ignore'   // workflow continues even if this process fails
    maxRetries 1             // optional: do not retry failed tasks

    script:
    """
    # get tumor and normal sample name
    tumor_sample=\$( grep "^#CHROM" ${expr_transcript_vcf} | awk '{print \$NF}' )
    normal_sample=\$( grep "^#CHROM" ${expr_transcript_vcf} | awk '{print \$(NF-1)}' )
    
    # bgzip and index the vcf
    bgzip ${expr_transcript_vcf} -o "${expr_transcript_vcf}.gz"
    tabix -p vcf "${expr_transcript_vcf}.gz"

    # create ouput directory
    mkdir -p pvacsplice
    
    # run pvacsplice
    pvacsplice run \\
        ${junctions} \\
        \${tumor_sample} \\
        "${hla}" \\
        all \\
        pvacsplice \\
        "${expr_transcript_vcf}.gz" \\
        ${dna_fasta} \\
        ${dna_gtf} \\
        --n-threads ${params.cpus} \\
        --normal-sample-name \${normal_sample} \\
        -e1 8,9,10,11 \\
        -e2 12,13,14,15,16,17,18 \\
        --pass-only \\
        --iedb-install-directory /opt/iedb/ \\
        > pvacsplice.log 2>&1 || true

     # Check if pvacsplice generated output, or exit silently
        if grep -q "No valid" pvacsplice.log; then \\
         echo "No valid mutations found for ${patient_id}, skipping."
         exit 0
        fi

    # Also ensure output directory exists even if empty
    touch pvacsplice/.placeholder

    """
}
