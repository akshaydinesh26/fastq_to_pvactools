process vep_annotate {
    
    tag "$patient_id"
    publishDir "${params.outdir}/${patient_id}/final_annotation", mode: 'copy'

    input:
    tuple val(patient_id), path(mutect_vcf)

    output:
    tuple val(patient_id), path("${patient_id}_vep.vcf"), emit: vep_vcf
    path "vep.log", emit: log_file

    // Error handling
    errorStrategy 'ignore'   // workflow continues even if this process fails
    maxRetries 1             // optional: do not retry failed tasks
    
    script:

    """
    vep -i "$mutect_vcf" -o "${patient_id}_vep.vcf" \
    --format vcf --vcf --symbol --protein --canonical --terms SO --tsl --biotype --hgvs \
    --fasta /data/homo_sapiens/113_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa \
    --gene_version --transcript_version \
    --offline --dir_cache /data \
    --plugin Frameshift --plugin Wildtype \
    --dir_plugins /data/Plugins > vep.log 2>&1
    """
}

