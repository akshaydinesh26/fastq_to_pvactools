process create_empty_files {
    
   output:
    tuple path('empty.vcf'), path('empty.cram'),path('empty.bam'),path('empty.rsem_gene'),path('empty.rsem_tx'),path('empty.fusion'), emit: empty_files_ch
    """
    touch empty.vcf empty.cram empty.bam empty.rsem_gene empty.rsem_tx empty.fusion
    """
}
