process vepitope_predictions_mixmhc {
    errorStrategy 'ignore'

    tag "$patient_id"
    publishDir "${params.outdir}/${patient_id}/mixmhc", mode: 'copy'
    
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

    #function to pair MHC class II alleles for BA predictor
    pair_all() {
        local alpha="\$1"
        local beta="\$2"
        local result=""
        IFS=',' read -ra ALPHA_ARR <<< "\$alpha"
        IFS=',' read -ra BETA_ARR <<< "\$beta"
        for a in "\${ALPHA_ARR[@]}"; do
            for b in "\${BETA_ARR[@]}"; do
                result+="\${a}__\${b},"
            done
        done
        echo "\${result%,}"
    }


    # create output folder
    #mkdir -p "vepitope_prediction"

    # check viral count exist
    if [ ! -s "$virus_star_count" ]; then
    echo "Exiting: No virus STAR counts for $patient_id. Skipping vepitope prediction."
    exit 0
    fi
    
    
    mkdir -p "vepitope_prediction/${patient_id}_mixmhc2pred"
    mkdir -p "vepitope_prediction/${patient_id}_prime"
    mkdir -p "vepitope_prediction/${patient_id}_subpeptides"
    mkdir -p "vepitope_prediction/${patient_id}_final_results"
    mkdir -p "vepitope_prediction/prioritized_epitopes"

    if [ $genome == "GRCH38" ]; then
       PHANTOM_VIRAL_SCRIPT_LOC="/opt/pyscripts/phantom_viral.py"
    else
       PHANTOM_VIRAL_SCRIPT_LOC="/opt/pyscripts/phantom_viral_mixmhc.py"
    fi

    python3 "\$PHANTOM_VIRAL_SCRIPT_LOC" QuantifyVirusCounts \
    -counts "$virus_star_count" \
    -luttable "$virus_reference_lutfile" \
    -human_luttable "$reference_lutfile" \
    -stranded "$strandness" \
    -humancounts "$human_star_count" \
    -output_dir vepitope_prediction \
    -protein_fasta "$virus_proteome"

    # if any proteins found generate subpeptides - which not found in human proteome
    # default - -lengths_str1 8,9,10,11 -lengths_str2 12,13,14,15,16,17,18
    # default no -problematic_n_str -problematic_c_str add if necessary
    if [ -e "vepitope_prediction/viral_genes_found.fasta" ]; then
        python3 /opt/pyscripts/get_subpeptides.py fasta_to_subpeptides \
        -input_fasta "vepitope_prediction/viral_genes_found.fasta" \
        -output_fasta1 "vepitope_prediction/${patient_id}_subpeptides/MHC_class_I_peptides.fasta" \
        -output_fasta2 "vepitope_prediction/${patient_id}_subpeptides/MHC_class_II_peptides.fasta" \
        -output_table "vepitope_prediction/${patient_id}_subpeptides/Generated_peptide_list.tsv" \
        -human_proteome "${human_proteome}"
    else
        echo "Exiting: no virus identified fasta file found"
        exit 0
    fi

    # prepare alleles
    
    #variable to track progress
    run_prime=0
    run_mixmhc2pred=0

    MHC_class_I=""
    MHC_class_II_final=""

    #check if alleles are empty
    if [[ -n "$hla_alleles" && "$hla_alleles" =~ HLA-[ABC] ]]; then

        # class I
        MHC_class_I=\$(echo "${hla_alleles}" | tr ',' '\n' \
        |sed 's/^[[:space:]]*//' \
        | grep -E '^(HLA-A|HLA-B|HLA-C)' \
        | sed 's/^HLA-//' \
        | tr -d '*:' \
        | paste -sd ",")

        # run PRIME
        /opt/PRIME-2.0/PRIME -i "vepitope_prediction/${patient_id}_subpeptides/MHC_class_I_peptides.fasta" \
        -o "vepitope_prediction/${patient_id}_prime/prime_results.txt" -a "\${MHC_class_I}" -mix "/opt/MixMHCpred-3.0/MixMHCpred"

        run_prime=1
    else
        echo "Exiting: MHC class I allele defined"
    fi

    if [[ -n "$hla_alleles" && "$hla_alleles" =~ D(R|Q|P) ]]; then
        
        #class II
        processed_class_II=\$(echo "${hla_alleles}" | tr ',' '\n' \
        | sed 's/^[[:space:]]*//' \
        | grep -E '^(DR|DQ|DP)' \
        | sed 's/[*:]/_/g')

        DPA=\$(echo "\$processed_class_II" | grep '^DPA' | paste -sd "," -)
        DPB=\$(echo "\$processed_class_II" | grep '^DPB' | paste -sd "," -)
        DQA=\$(echo "\$processed_class_II" | grep '^DQA' | paste -sd "," -)
        DQB=\$(echo "\$processed_class_II" | grep '^DQB' | paste -sd "," -)
        DR=\$(echo "\$processed_class_II" | grep '^DR' | paste -sd "," -)


        DP_pairs=\$(pair_all "\$DPA" "\$DPB")
        DQ_pairs=\$(pair_all "\$DQA" "\$DQB")
        #echo "\$DQ_pairs"

        MHC_class_II_final="\$DP_pairs,\$DQ_pairs"
        [ -n "\$DR" ] && MHC_class_II_final="\$MHC_class_II_final,\$DR"
        MHC_class_II_final=\$(echo "\$MHC_class_II_final" | sed 's/^,*//;s/,*\$//' | tr ',' ' ')

        # run MixMHC2pred
        MixMHC2pred --input "vepitope_prediction/${patient_id}_subpeptides/MHC_class_II_peptides.fasta" \
        --output "vepitope_prediction/${patient_id}_mixmhc2pred/mixMHC2pred_results.txt" --alleles \$MHC_class_II_final --no_context

        run_mixmhc2pred=1
    else
        echo "Exiting: MHC class II allele defined"
    fi

    prime_pred=""
    mixmhc2pred_pred=""

    if [ -s "vepitope_prediction/${patient_id}_prime/prime_results.txt" ]; then
        prime_pred="vepitope_prediction/${patient_id}_prime/prime_results.txt"
    fi

    if [ -s "vepitope_prediction/${patient_id}_mixmhc2pred/mixMHC2pred_results.txt" ]; then
        mixmhc2pred_pred="vepitope_prediction/${patient_id}_mixmhc2pred/mixMHC2pred_results.txt"
    fi    


    if [[ "\$run_prime" -eq 1 || "\$run_mixmhc2pred" -eq 1 ]]; then
        echo "2 \${MHC_class_II_final}"
        echo "1 \${MHC_class_I}"
    
        # create peptide list table
        python3 /opt/pyscripts/create_final_table_mixmhc.py create_peptide_final_list \
            -peptide_table "vepitope_prediction/${patient_id}_subpeptides/Generated_peptide_list.tsv" \
            -luttable "${virus_reference_lutfile}" \
            -virus_summary_final "vepitope_prediction/identified_virus.txt" \
            -outdir "vepitope_prediction/${patient_id}_final_results"

        # filter the results to remove non binders
        awk -F'\t' '{ if( \$0 ~ /^#/ ){ print \$0 } else if(\$1 ~/^Peptide/){print \$0} else if(\$4 <= 2.0){ print \$0 } }' \$prime_pred > "vepitope_prediction/${patient_id}_prime/prime_results_filtered.txt"
        awk -F'\t' '{ if( \$0 ~ /^#/ ){ print \$0 } else if(\$1 ~/^Peptide/){print \$0} else if(\$4 <= 2.0){ print \$0 } }' \$mixmhc2pred_pred > "vepitope_prediction/${patient_id}_mixmhc2pred/mixMHC2pred_results_filtered.txt"

        # combine predictions
        python3 /opt/pyscripts/create_final_table_mixmhc.py CombinePrediction \
            -table "vepitope_prediction/${patient_id}_final_results/peptide_final_table.tsv" \
            -cI_preds "vepitope_prediction/${patient_id}_prime/prime_results_filtered.txt" \
            -cII_preds "vepitope_prediction/${patient_id}_mixmhc2pred/mixMHC2pred_results_filtered.txt" \
            -BindingScore 1 -ScoreThreshold 2 > "vepitope_prediction/${patient_id}_final_results/combined_peptide_prediction.tsv"

        # add immunogenic annotation
        python3 /opt/pyscripts/create_final_table_mixmhc.py AddAnnotation \
            -table "vepitope_prediction/${patient_id}_final_results/combined_peptide_prediction.tsv" \
            -vepitopes "${immunogenic_peptide_info}" > "vepitope_prediction/${patient_id}_final_results/annotated_results.tsv"

        # prioritization
        python3 /opt/pyscripts/create_final_table_mixmhc.py PrioritizeViralPeptides \
            -table "vepitope_prediction/${patient_id}_final_results/annotated_results.tsv" \
            -rank 1 > "vepitope_prediction/${patient_id}_final_results/final_annotated_sorted.results.tsv"
    else
        echo "Exiting: MHC class II allele defined"
        exit 0
    fi
    """
}
