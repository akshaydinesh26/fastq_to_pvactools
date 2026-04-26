process vepitope_predictions {

    tag "$patient_id"
    publishDir "${params.outdir}/${patient_id}/pvacbind", mode: 'copy'
    
    input:
    tuple val(patient_id), path(virus_star_count), path(human_star_count), val(strandness), val(hla_alleles)
    path(virus_reference_lutfile)
    path(reference_lutfile)
    path(virus_proteome)
    path(human_proteome)
    path(immunogenic_peptide_info)
    val(genome)


    output:
    tuple val(patient_id), path("vepitope_prediction")
    script:
    """
    # create output folder
    mkdir -p "vepitope_prediction"

     # check viral count exist
    if [ ! -s "$virus_star_count" ]; then
    echo "Exiting: No virus STAR counts for $patient_id. Skipping vepitope prediction."
    exit 0
    fi
    
    
    mkdir -p "vepitope_prediction/${patient_id}_pvacbind"
    mkdir -p "vepitope_prediction/prioritized_epitopes"

    if [ $genome == "GRCH38" ]; then
       PHANTOM_VIRAL_SCRIPT_LOC="/pyscripts/phantom_viral.py" 
    else
       PHANTOM_VIRAL_SCRIPT_LOC="/pyscripts/phantom_viral_mixmhc.py"
    fi 

    # identify virus and generate fasta file of identified proteins
    python3 "\$PHANTOM_VIRAL_SCRIPT_LOC" QuantifyVirusCounts \
    -counts "$virus_star_count" \
    -luttable "$virus_reference_lutfile" \
    -human_luttable "$reference_lutfile" \
    -stranded "$strandness" \
    -humancounts "$human_star_count" \
    -output_dir vepitope_prediction \
    -protein_fasta "$virus_proteome"

    # check if any proteins found expressing
    if [ -e "vepitope_prediction/viral_genes_found.fasta" ]; then
        pvacbind run --iedb-install-directory /opt/iedb -t ${task.cpus} \
        -e1 8,9,10,11,12 -e2 12,13,14,15,16,17,18 \
        "vepitope_prediction/viral_genes_found.fasta" "${patient_id}_virus" "$hla_alleles" \
        all "vepitope_prediction/${patient_id}_pvacbind"
    else
        echo "Exiting: no virus identified fasta file found"
        exit 0
    fi

    # if only hla I take all epitope from MHC_class_I folder
    if [[ "${hla_alleles}" =~ (DP|DR|DQ) ]]; then
	epitope_prediction="vepitope_prediction/${patient_id}_pvacbind/MHC_Class_I/${patient_id}_virus.all_epitopes.tsv"
    else
        epitope_prediction="vepitope_prediction/${patient_id}_pvacbind/combined/${patient_id}_virus.all_epitopes.tsv"
    fi   

    # generate final epitope prioritized list
    if [ -e "\$epitope_prediction" ]; then
    python3 /pyscripts/create_finaltable.py createFinalTable \
        -allepitopes "\$epitope_prediction" \
        -luttable ${virus_reference_lutfile} \
        -immunelist "$immunogenic_peptide_info" \
        -virus_summary_final "vepitope_prediction/identified_virus.txt" \
        -human_fasta "$human_proteome" \
        -output_dir "vepitope_prediction/prioritized_epitopes"
    mv "vepitope_prediction/prioritized_epitopes/final_prioritized_predictions.tsv" "vepitope_prediction/prioritized_epitopes/${patient_id}_final_prioritized_predictions.tsv"
    mv "vepitope_prediction/prioritized_epitopes/annotated_predictions.tsv" "vepitope_prediction/prioritized_epitopes/${patient_id}_annotated_predictions.tsv"
    mv "vepitope_prediction/identified_virus.txt" "vepitope_prediction/${patient_id}_identified_virus.txt"
    mv "vepitope_prediction/viral_genes_found.fasta" "vepitope_prediction/${patient_id}_viral_genes_found.fasta"
    else
        echo "Exiting: no pvacbind result found"
        exit 0
    fi    
    """
}
