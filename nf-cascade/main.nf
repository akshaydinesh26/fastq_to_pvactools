// include all workflows and functions - for nfcore, custom workflows and functions
include { NEXTFLOW_RUN as NFCORE_SAREK          } from "$projectDir/modules/local/nextflow/run/main"
include { NEXTFLOW_RUN as NFCORE_RNASEQ      } from "$projectDir/modules/local/nextflow/run/main"
include { NEXTFLOW_RUN as NFCORE_RNAFUSION        } from "$projectDir/modules/local/nextflow/run/main"
include { NEXTFLOW_RUN as NFCORE_HLATYPING        } from "$projectDir/modules/local/nextflow/run/main"
include { NEXTFLOW_TRONFLOW_RUN as TRONFLOW        } from "$projectDir/modules/local/nextflow/run/main"
include { NEXTFLOW_NEOANTIGEN_RUN as NEOANTIGEN        } from "$projectDir/modules/local/nextflow/run/main"
include { PHANTOM_VIRAL_RUN as VIRAL_PEP        } from "$projectDir/modules/local/nextflow/run/main"
include { readWithDefault                      } from "$projectDir/functions/local/utils"
include { resolveFileFromDir as getSamplesheet } from "$projectDir/functions/local/utils"
include { readCsv                              } from "$projectDir/functions/local/utils"
include { writeTemplates                       } from "$projectDir/functions/local/utils"
include { writeTemplatesTsv                       } from "$projectDir/functions/local/utils"
include { getHLAHDFileFromDir                  } from "$projectDir/functions/local/utils"
include { getOptiTypeFileFromDir               } from "$projectDir/functions/local/utils"
include { parseHLAFromFiles                    } from "$projectDir/functions/local/utils"
include { getvcfFileFromDir                    } from "$projectDir/functions/local/utils"
include { gettumorcramFileFromDir              } from "$projectDir/functions/local/utils"
include { getrnabamFileFromDir                 } from "$projectDir/functions/local/utils"
include { getgenecountFileFromDir              } from "$projectDir/functions/local/utils"
include { gettrcountFileFromDir                } from "$projectDir/functions/local/utils"
include { getarribaFileFromDir                 } from "$projectDir/functions/local/utils"





workflow {

    if( !params.input_csv ) {
    exit 1, "ERROR: You must provide --input_csv <file>"
    }

    def samples = readCsv(params.input_csv)
    
    // Generate per-patient template CSVs and return maps of csv files for each patients mapped with patient id
    def sarekInputs       = writeTemplates(samples, params.sarek_template,     "${params.pipeline_input}", "sarek")
    def rnaseqInputs      = writeTemplates(samples, params.rnaseq_template,    "${params.pipeline_input}", "rnaseq")
    def rnafusionInputs   = writeTemplates(samples, params.rnafusion_template, "${params.pipeline_input}", "rnafusion")
    def hlatypingInputs   = writeTemplates(samples, params.hlatyping_template, "${params.pipeline_input}", "hlatyping")
    def tronflowInputs    = writeTemplates(samples, params.tronflow_template, "${params.pipeline_input}", "tronflow")
        

    // loop through each patient and updates worflow input file paramters
    // then runs all the worflows per patients
    // output the results to per patient folder defined by updated params.outdir
    samples*.patient.unique().each { patient ->
    
    log.info "=== Running all workflows for patient: ${patient} ==="

    // create output directory
    // get the single sample for this patient
    def patient_sample     = samples.find { it.patient == patient }
    def patient_outdir     = "${params.output_final}"
    def sarek_outdir       = "${patient_outdir}/sarek"
    def rnaseq_outdir      = "${patient_outdir}/rnaseq"
    def rnafusion_outdir   = "${patient_outdir}/rnafusion"
    def hlatyping_outdir   = "${patient_outdir}/hlatyping"
    def tronflow_outdir    = "${patient_outdir}/tronflow"
    def neoantigen_outdir  = "${patient_outdir}/phantom_mini"
    def phantom_viral_outdir       = "${patient_outdir}/phantom_viral"

    // define worflow
    def workflow_to_process = "${params.workflows}"

    // Check if 'hla' is defined and non-empty
    if (patient_sample?.hla?.trim()) {
        log.info "HLA data is present for patient: ${patient}  ~ skipping HLA-related workflows"
        workflow_to_process = 'sarek,rnaseq,rnafusion,neoantigen,virus'
    } else {
        log.warn "No HLA data for patient: ${patient} — doing HLA-related workflows"
        workflow_to_process = 'sarek,rnaseq,rnafusion,hlatyping,tronflow,neoantigen,virus'
    }

    if (params.custom_workflow != "none") {
        log.warn "Custom workflow found: ${patient} — overriding and executimng custom workflow"
        workflow_to_process = params.custom_workflow
    }

    // Update params for this patient
    params.sarek.input             = sarekInputs[patient]
    params.rnaseq.input            = rnaseqInputs[patient]
    params.rnafusion.input         = rnafusionInputs[patient]
    params.hlatyping.input         = hlatypingInputs[patient]
    params.tronflow.input_fastqs   = tronflowInputs[patient]

    // Run workflows in sequence for this patient
    // Validate possible pipeline chains
    def valid_chains = [
        'sarek',
        'rnaseq',
        'rnafusion',
        'sarek,rnaseq',
        'sarek,rnafusion',
        'sarek,rnaseq,rnafusion',
        'rnaseq,rnafusion',
        'virus',
        'hlatyping,tronflow',
        'virus,hlatyping,tronflow',
        'sarek,rnaseq,rnafusion,virus',
        'sarek,rnaseq,rnafusion,hlatyping,tronflow',
        'sarek,rnaseq,rnafusion,hlatyping,tronflow,virus',
        'sarek,rnaseq,rnafusion,neoantigen',
        'sarek,rnaseq,rnafusion,neoantigen,virus',
        'sarek,rnaseq,rnafusion,hlatyping,tronflow,neoantigen',
        'sarek,rnaseq,rnafusion,hlatyping,tronflow,neoantigen,virus'
    ]
    assert workflow_to_process in valid_chains
    def wf_chain = workflow_to_process.tokenize(',')

    // Initialise undefined channels
    def hlatyping_ch                   = null
    def tronflow_ch                    = null
    def rnaseq_ch                      = null
    def rnafusion_ch                   = null
    def sarek_ch                       = null
    def neoantigen_ch                  = null
// def neoantigen_ch                  = Channel.value([])

    // Run pipelines
   if ( 'sarek' in wf_chain ){
        // SAREK
        NFCORE_SAREK (
            'nf-core/sarek',
            "${params.general.wf_opts?: ''} ${params.sarek.wf_opts?: ''}",  // workflow opts
            readWithDefault( params.sarek.params_file, Channel.value([]) ), // params file
            readWithDefault( params.sarek.input, Channel.value([]) ),       // samplesheet
            readWithDefault( params.sarek.add_config, Channel.value([]) ),  // custom config
            workflow.workDir.resolve('nf-core/sarek').toUriString(),
            sarek_outdir
        )
        sarek_ch = NFCORE_SAREK.out.output  ?: Channel.value([])
    }

   if ('rnaseq' in wf_chain ){
        // RNASEQ
       NFCORE_RNASEQ (
           'nf-core/rnaseq',
           "${ params.general.wf_opts?: ''} ${params.rnaseq.wf_opts?: ''}",     // workflow opts
           readWithDefault( params.rnaseq.params_file, Channel.value([]) ),     // params file
           readWithDefault( params.rnaseq.input, Channel.value([]) ), // samplesheet
           readWithDefault( params.rnaseq.add_config, Channel.value([]) ),      // custom config
           workflow.workDir.resolve('nf-core/rnaseq').toUriString(),
           rnaseq_outdir
       )
       rnaseq_ch = NFCORE_RNASEQ.out.output  ?: Channel.value([])
   }
   
   if ('rnafusion' in wf_chain ){
        // RNAFUSION
       NFCORE_RNAFUSION (
           'nf-core/rnafusion',
           "${ params.general.wf_opts?: ''} ${params.rnafusion.wf_opts?: ''}",     // workflow opts
           readWithDefault( params.rnafusion.params_file, Channel.value([]) ),     // params file
           readWithDefault( params.rnafusion.input, Channel.value([]) ), // samplesheet
           readWithDefault( params.rnafusion.add_config, Channel.value([]) ),      // custom config
           workflow.workDir.resolve('nf-core/rnafusion').toUriString(),
           rnafusion_outdir
       )
       rnafusion_ch = NFCORE_RNAFUSION.out.output  ?: Channel.value([])
   }

   if ('hlatyping' in wf_chain ){
        // OPTITYPE
       NFCORE_HLATYPING (
           'nf-core/hlatyping',
           "${ params.general.wf_opts?: ''} ${params.hlatyping.wf_opts?: ''}",     // workflow opts
           readWithDefault( params.hlatyping.params_file, Channel.value([]) ),     // params file
           readWithDefault( params.hlatyping.input, Channel.value([]) ), // samplesheet
           readWithDefault( params.hlatyping.add_config, Channel.value([]) ),      // custom config
           workflow.workDir.resolve('nf-core/hlatyping').toUriString(),
           hlatyping_outdir
       )
       hlatyping_ch = NFCORE_HLATYPING.out.output ?: Channel.value([])
   }

   if ('tronflow' in wf_chain ){
        // TRONFLOW
       TRONFLOW (
           "${projectDir}/tronflow-hla-hd/main.nf",
           "${ params.general.wf_opts?: ''} ${params.tronflow.wf_opts?: ''}",     // workflow opts
           readWithDefault( params.tronflow.params_file, Channel.value([]) ),     // params file
           readWithDefault( params.tronflow.input_fastqs, Channel.value([]) ), // samplesheet
           readWithDefault( params.tronflow.add_config, Channel.value([]) ),      // custom config
           workflow.workDir.resolve('tronflow').toUriString(),
           tronflow_outdir
       )
       tronflow_ch  = TRONFLOW.out.output  ?: Channel.value([])
   }

   
   
//    create default variables
   def hlaMap = null
   def neoantigen_input = []
//    def neoantigenInputs = []
   
   if (['hlatyping', 'tronflow', 'neoantigen'].every { it in wf_chain }) {
    
    // Combine the outputs via channel (DSL2 way) takeout output files for tronflow and hlatyping and get HLAstring
    // combined_ch = hlatyping_ch.combine(
    //     tronflow_ch,
    //     sarek_ch,
    //     rnaseq_ch,
    //     rnafusion_ch).filter {
    //         hlatyping, tronflow, sarek, rnaseq,rnafusion ->
    //         // check all are non empty list or non-null
    //         [hlatyping, tronflow, sarek, rnaseq, rnafusion].every {it && !it.isEmpty() }
    //     }


    hlatyping_ch.view()
    tronflow_ch.view()

    combined_ch = hlatyping_ch
        .combine(tronflow_ch)
        .combine(sarek_ch)
        .combine(rnaseq_ch)
        .combine(rnafusion_ch)
    
    combined_ch.map {tronDir, hlaDir, sarekDir, rnaseqDir, rnafusionDir ->


    def hlaHDFile    = getHLAHDFileFromDir(hlaDir,"${patient}")
    def optiTypeFile = getOptiTypeFileFromDir(tronDir,"${patient}")
    def vcf          = getvcfFileFromDir(sarekDir,"${patient}")
    def tumorCram    = gettumorcramFileFromDir(sarekDir,"${patient}")
    def rnaBam       = getrnabamFileFromDir(rnaseqDir,"${patient}")
    def geneCount    = getgenecountFileFromDir(rnaseqDir,"${patient}")
    def trCount      = gettrcountFileFromDir(rnaseqDir,"${patient}")
    def fusionFile   = getarribaFileFromDir(rnafusionDir,"${patient}")
    
    println "HLA-HD file: ${hlaHDFile?.absolutePath ?: 'Not found'}"
    println "OptiType file: ${optiTypeFile?.absolutePath ?: 'Not found'}"
    println "vcf file: ${vcf?.absolutePath ?: 'Not found'}"
    println "tumorCram file: ${tumorCram?.absolutePath ?: 'Not found'}"
    println "rnaBam file: ${rnaBam?.absolutePath ?: 'Not found'}"
    println "geneCount file: ${geneCount?.absolutePath ?: 'Not found'}"
    println "trCount file: ${trCount?.absolutePath ?: 'Not found'}"
    println "fusionFile file: ${fusionFile?.absolutePath ?: 'Not found'}"

    
    
    // get hla
    hlaMap = parseHLAFromFiles("${patient}", hlaHDFile, optiTypeFile)
    def hlaString = hlaMap[patient] ?: ""
     neoantigen_input << [
        patient: patient,
        hlaString : hlaString,
        vcf: vcf.toString(),
        tumorCram: tumorCram.toString(),
        rnaBam: rnaBam.toString(),
        geneCount: geneCount.toString(),
        trCount: trCount.toString(),
        fusionFile: fusionFile.toString()
        ]

    def neoantigenInputs = writeTemplatesTsv(neoantigen_input, params.neoantigen_template, "${params.pipeline_input}", "neoantigen")
    return neoantigenInputs[patient]
    

        }.set{neoantigen_tsv}

    
    } else if ('neoantigen' in wf_chain && !('hlatyping' in wf_chain) && !('tronflow' in wf_chain)) {
        hlaString = patient_sample.hla

        combined_ch = sarek_ch
        .combine(rnaseq_ch)
        .combine(rnafusion_ch)

        // combined_ch = sarek_ch.combine(
        // rnaseq_ch,
        // rnafusion_ch).filter {
        //     sarek, rnaseq,rnafusion ->
        //     // check all are non empty list or non-null
        //     [sarek, rnaseq, rnafusion].every {it && !it.isEmpty() }
        // }
        
        combined_ch.map {sarekDir, rnaseqDir, rnafusionDir ->

        def vcf          = getvcfFileFromDir(sarekDir,"${patient}")
        def tumorCram    = gettumorcramFileFromDir(sarekDir,"${patient}")
        def rnaBam       = getrnabamFileFromDir(rnaseqDir,"${patient}")
        def geneCount    = getgenecountFileFromDir(rnaseqDir,"${patient}")
        def trCount      = gettrcountFileFromDir(rnaseqDir,"${patient}")
        def fusionFile   = getarribaFileFromDir(rnafusionDir,"${patient}")
        
        
        println "vcf file: ${vcf?.absolutePath ?: 'Not found'}"
        println "tumorCram file: ${tumorCram?.absolutePath ?: 'Not found'}"
        println "rnaBam file: ${rnaBam?.absolutePath ?: 'Not found'}"
        println "geneCount file: ${geneCount?.absolutePath ?: 'Not found'}"
        println "trCount file: ${trCount?.absolutePath ?: 'Not found'}"
        println "fusionFile file: ${fusionFile?.absolutePath ?: 'Not found'}"
        
        
        
        neoantigen_input << [
            patient: patient,
            hlaString : hlaString,
            vcf: vcf.toString(),
            tumorCram: tumorCram.toString(),
            rnaBam: rnaBam.toString(),
            geneCount: geneCount.toString(),
            trCount: trCount.toString(),
            fusionFile: fusionFile.toString()
        ]
        def neoantigenInputs = writeTemplatesTsv(neoantigen_input, params.neoantigen_template, "${params.pipeline_input}", "neoantigen")
        return neoantigenInputs[patient]
        // if ('neoantigen' in wf_chain ){
            // NEOANTIGEN
        
//    }

        // neoantigen_ch
        // .collect()
        // .subscribe { items ->
        //     if (items && !items.isEmpty()) {
        //         println "Neoantigen analysis finished (${items.size()} items)"
        //     } else {
        //         println "Warning: Neoantigen channel is empty"
        //     }
        // }

        }.set{neoantigen_tsv}

    } else {

        neoantigen_tsv = Channel.empty()
    }

    def neoantigen_final_opts = params.neoantigen.wf_opts ?: ''
    
    if (params.wes_only){
        neoantigen_final_opts +=  " --wes_only"
    }
    if (params.pvacseq_only){
        neoantigen_final_opts +=  " --pvacseq_only"
    }

    // if (file(neoantigen_tsv).exists()) {
       NEOANTIGEN (
            "${projectDir}/phantom_mini/phantom.nf",
            "${ params.general.wf_opts?: ''} ${neoantigen_final_opts?: ''}",     // workflow opts
            readWithDefault( params.neoantigen.params_file, Channel.value([]) ),     // params file
            neoantigen_tsv,// readWithDefault( neoantigen_tsv, Channel.value([]) ), // samplesheet
            readWithDefault( params.neoantigen.add_config, Channel.value([]) ),      // custom config
            workflow.workDir.resolve('neoantigen').toUriString(),
            neoantigen_outdir
        )

    

        if ('virus' in wf_chain ){

            Channel
             .of([patient_sample])
             .set { ch_patient_sample }

            if (patient_sample?.hla?.trim()) {
                ch_patient_sample.map {patient_sample_new_data ->

                def phantomviralInput = writeTemplatesTsv(patient_sample_new_data, params.phantom_viral_template, "${params.pipeline_input}", "phantom_viral")
                return phantomviralInput[patient]

                }.set{phantom_viral_input_ch}

            }else{
                
                ch_patient_sample
                    .combine(hlatyping_ch)
                    .combine(tronflow_ch)
                    .set {combined_input_viral_ch}

                combined_input_viral_ch.map{patient_sample_new_data, hlaDir, tronDir ->
                    def hlaHDFile    = getHLAHDFileFromDir(hlaDir,"${patient}")
                    def optiTypeFile = getOptiTypeFileFromDir(tronDir,"${patient}")
                    hlaMap = parseHLAFromFiles("${patient}", hlaHDFile, optiTypeFile)
                    def hlaString = hlaMap[patient] ?: ""
                    patient_sample.hla = hlaString

                    def phantomviralInput = writeTemplatesTsv(patient_sample_new_data, params.phantom_viral_template, "${params.pipeline_input}", "phantom_viral")
                    return phantomviralInput[patient]


                }.set{phantom_viral_input_ch}
            }

            


        // viral
       VIRAL_PEP (
            "${projectDir}/phantom_viral/phantom_viral.nf",
            "${ params.general.wf_opts?: ''} ${params.phantom_viral.wf_opts?: ''}",     // workflow opts
            readWithDefault( params.phantom_viral.params_file, Channel.value([]) ),     // params file
            phantom_viral_input_ch,// readWithDefault( neoantigen_tsv, Channel.value([]) ), // samplesheet
            readWithDefault( params.phantom_viral.add_config, Channel.value([]) ),      // custom config
            workflow.workDir.resolve('phantom_viral').toUriString(),
            phantom_viral_outdir
        )
        }

    // }
    }    


}




