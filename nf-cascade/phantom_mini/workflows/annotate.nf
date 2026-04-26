#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { vep_annotate }       from "${projectDir}/modules/vep_annotate.nf"
include { decompose_vcf }      from "${projectDir}/modules/decompose_vcf.nf"
include { generate_site_list } from "${projectDir}/modules/generate_site_list.nf"
include { vcf_annotate_all }   from "${projectDir}/modules/annotate_vcf.nf"
include { generate_junction }  from "${projectDir}/modules/generate_junction.nf"

workflow annotate {
    take:
      annotate_input
      dna_fasta_ch
      dna_gtf_ch

    main:
    annotate_input
        .map { patient_id, mutect_vcf, col3, col4, col5, col6, col7, col8 -> tuple(patient_id, mutect_vcf) }
        .set { vep_input }  

    // run vep annotation - create empty channel if failed
    vep_results = vep_annotate(vep_input)

    vep_vcf_ch = vep_results.vep_vcf ?: Channel.empty()

    // run vep decomposition
    vep_decomposed = decompose_vcf(vep_vcf_ch) 
    decompose_vcf_ch = vep_decomposed.decompose_vcf ?: Channel.empty()

    //  Prepare input channel for generate_site_list
    annotate_input
        .map { patient_id, col2, dna_cram, rna_bam, col5, col6, col7, col9 ->
            tuple(patient_id, dna_cram, rna_bam)
        }
        .set { generate_site_list_input_bam }
    
    generate_site_list_input_combined = decompose_vcf_ch
        .join(generate_site_list_input_bam)

    // run bam readcount for site list generation
    generated_site_list = generate_site_list(generate_site_list_input_combined,file(params.dna_fasta))
    dna_site_list_ch = generated_site_list.dna_site_list ?: Channel.empty()
    rna_site_list_ch = generated_site_list.rna_site_list ?: Channel.empty()

    // Prepare input channel for vcf annotation
    annotate_input
        .map { patient_id, col2, col3, col4, rsem_gene, rsem_tx, col7, col8 -> tuple(patient_id, rsem_gene, rsem_tx) }
        .set { expression_input }

    // annotate vcf with read count and expression
    
    vcf_annotate_input = decompose_vcf_ch
        .join(expression_input)
        .map {patient_id,decomposed_vcf,rsem_gene,rsem_tx -> 
            tuple(patient_id,decomposed_vcf,rsem_gene,rsem_tx)}
        .join(dna_site_list_ch)
        .map {patient_id,decomposed_vcf,rsem_gene,rsem_tx,site_list_dir_dna -> 
            tuple(patient_id,decomposed_vcf,rsem_gene,rsem_tx,site_list_dir_dna)}
        .join(rna_site_list_ch)
        .map {patient_id,decomposed_vcf,rsem_gene,rsem_tx,site_list_dir_dna,site_list_dir_rna ->
            tuple(patient_id,decomposed_vcf,rsem_gene,rsem_tx,site_list_dir_dna,site_list_dir_rna)}

    final_annotated_vcf = vcf_annotate_all(
        vcf_annotate_input
    )
    vcf_annotated_ch = final_annotated_vcf.vcf_annotated ?: Channel.empty()

    // Prepare input channel for junction annotation
    annotate_input
        .map { patient_id, col2, col3, rna_bam, col5, col6, col7, col9 ->
            tuple(patient_id,rna_bam)
        }
        .set { rna_bam_input } ?: Channel.empty()


    generate_junction_input_combined = rna_bam_input
        .join(vcf_annotated_ch)

    // Prepare input channel for hla
    annotate_input
        .map { patient_id, col2, col3, col4, col5, col6, hla_alleles, col8 -> tuple(patient_id, hla_alleles) }
        .set { hla_allele_input }


    final_annotated_vcf.vcf_annotated
        .map { patient_id, vcf_expr_transcript -> tuple(patient_id,vcf_expr_transcript) }
        .set { vcf_expr_transcript_input }

    pvacseq_input_combined = hla_allele_input 
        .join(vcf_expr_transcript_input) 
    pvacseq_input_combined_ch = pvacseq_input_combined ?: Channel.empty()
    pvacsplice_input_combined_ch = Channel.empty()

    if (!params.pvacseq_only) {

        // generate junction list
        splice_junctions = generate_junction(
            generate_junction_input_combined,
            params.dna_fasta,
            params.dna_gtf)
        junction_list_ch = splice_junctions.junction_list ?: Channel.empty()

        // pvacsplice input (conditional)
        pvacsplice_input_combined = hla_allele_input
                .join(vcf_expr_transcript_input)
                .join(junction_list_ch)

        pvacsplice_input_combined_ch = pvacsplice_input_combined  
        
    }

    
   
    
    emit:
      pvacsplice_input_combined_ch
      pvacseq_input_combined_ch
}
