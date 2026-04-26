/**
 * Returns a channel with the path if it's defined, otherwise returns a default channel.
 * 
 * @param path             The path to include into the channel
 * @param default_channel  A channel to use as the default if no path is defined.
 * @return                 A channel with a path, or the default channel
 */
def readWithDefault( String path, Object default_channel ) {
    path ? Channel.fromPath( path, checkIfExists: true ) : default_channel
}

/**
 * Returns a channel with the file defined by the path resolved against the directory.
 * 
 * @param path  The path of the file relative to the directory in dir
 * @param dir   A channel with a directory.
 * @return      A channel with a path relative to the dir path
 */
def resolveFileFromDir ( String path, Object dir ){
    dir.map{ results -> file( results.resolve( path ) ) }
}

/**
 * Returns a List of maps 
 * 
 * For Each patient in the input csv file there is a map
 */
def readCsv(path) {
    def lines = file(path).readLines()
    def header = lines[0].split('\t').collect { it.trim() }
    lines.drop(1).collect { line ->
        def cols = line.split('\t').collect { it.trim() }
        [header, cols].transpose().collectEntries()
    }
}

/**
 * Returns a List of input csv files for the pipelines and stores the files per patient in a map
 * 
 * For Each patient a map of input csv files to pipelines
 */
def writeTemplates(samples, templateFile, outputDir, workflowName) {
    def templateText = file(templateFile).text
    def outdir = file(outputDir)
    outdir.mkdirs()

    def patientFiles = [:]  // patient -> file path
    samples.each { row ->
        def rendered = templateText
        row.each { key, value ->
            rendered = rendered.replace("\${${key}}", value)
        }
        def outfile = file("${outdir}/${row.patient}_${workflowName}.csv")
        outfile.text = rendered
        patientFiles[row.patient] = outfile.toString()
    }
    return patientFiles
}

def writeTemplatesTsv(samples, templateFile, outputDir, workflowName) {
    def templateText = file(templateFile).text
    def outdir = file(outputDir)
    outdir.mkdirs()

    def patientFiles = [:]  // patient -> file path
    samples.each { row ->
        def rendered = templateText
        row.each { key, value ->
            rendered = rendered.replace("\${${key}}", value)
        }
        def outfile = file("${outdir}/${row.patient}_${workflowName}.tsv")
        outfile.text = rendered
        patientFiles[row.patient] = outfile.toString()
    }
    return patientFiles
}
/**
 * Returns a  file object* 
 * For Each patient a map of input csv files to pipelines
 */
// For HLA-HD file
def getHLAHDFileFromDir(emittedDir, patientId) {
    def baseDir = emittedDir instanceof File ? emittedDir : new File(emittedDir.toString())
    // If the file is inside the baseDir/patientId subfolder
    def filePath = new File(baseDir, "${patientId}/${patientId}_final.result.txt")
    if (filePath.exists()) {
        return filePath
    } else {
        println "Warning: HLA-HD file not found at ${filePath}"
        return null
    }
}

// For OptiType file
def getOptiTypeFileFromDir(emittedDir, patientId) {
    def baseDir = emittedDir instanceof File ? emittedDir : new File(emittedDir.toString())
    // If the file is inside baseDir/patientId subfolder
    def filePath = new File(baseDir, "optitype/${patientId}/${patientId}_result.tsv")
    if (filePath.exists()) {
        return filePath
    } else {
        println "Warning: OptiType file not found at ${filePath}"
        return null
    }
}

// For filtered vcf file
def getvcfFileFromDir(emittedDir, patientId) {
    def baseDir = emittedDir instanceof File ? emittedDir : new File(emittedDir.toString())
    // If the file is inside baseDir/patientId subfolder
    def filePath = new File(baseDir, "variant_calling/mutect2/${patientId}/${patientId}.mutect2.filtered.vcf.gz")
    if (filePath.exists()) {
        return filePath
    } else {
        println "Warning: filtered vcf file not found at ${filePath}"
        return null
    }
}

// For tumor cram file
def gettumorcramFileFromDir(emittedDir, patientId) {
    def baseDir = emittedDir instanceof File ? emittedDir : new File(emittedDir.toString())
    // If the file is inside baseDir/patientId subfolder
    def filePath = new File(baseDir, "preprocessing/recalibrated/${patientId}T/${patientId}T.recal.cram")
    if (filePath.exists()) {
        return filePath
    } else {
        println "Warning: tumor cram file not found at ${filePath}"
        return null
    }
}


// For rna bam file
def getrnabamFileFromDir(emittedDir, patientId) {
    def baseDir = emittedDir instanceof File ? emittedDir : new File(emittedDir.toString())
    // If the file is inside baseDir/patientId subfolder
    def filePath = new File(baseDir, "star_rsem/${patientId}.markdup.sorted.bam")
    if (filePath.exists()) {
        return filePath
    } else {
        println "Warning: rna bam file not found at ${filePath}"
        return null
    }
}

// For rna gene results file
def getgenecountFileFromDir(emittedDir, patientId) {
    def baseDir = emittedDir instanceof File ? emittedDir : new File(emittedDir.toString())
    // If the file is inside baseDir/patientId subfolder
    def filePath = new File(baseDir, "star_rsem/${patientId}.genes.results")
    if (filePath.exists()) {
        return filePath
    } else {
        println "Warning: gene count results file not found at ${filePath}"
        return null
    }
}

// For rna transcript results file
def gettrcountFileFromDir(emittedDir, patientId) {
    def baseDir = emittedDir instanceof File ? emittedDir : new File(emittedDir.toString())
    // If the file is inside baseDir/patientId subfolder
    def filePath = new File(baseDir, "star_rsem/${patientId}.isoforms.results")
    if (filePath.exists()) {
        return filePath
    } else {
        println "Warning: transcript count results file not found at ${filePath}"
        return null
    }
}


// For rna fusion results file
def getarribaFileFromDir(emittedDir, patientId) {
    def baseDir = emittedDir instanceof File ? emittedDir : new File(emittedDir.toString())
    // If the file is inside baseDir/patientId subfolder
    def filePath = new File(baseDir, "arriba/${patientId}.arriba.fusions.tsv")
    if (filePath.exists()) {
        return filePath
    } else {
        println "Warning: arriba results file not found at ${filePath}"
        return null
    }
}



/**
 * Returns a object patient:"HLAstring"* 
 * For Each patient takes in TRONFLOW AND NFCORE HLATYPING output and returns hlastring
 */
// For HLA-HD file
def parseHLAFromFiles(patientId, hlahdFile, optitypeFile) {
    def result = [:]

    if (!hlahdFile.exists()) {
        println "Warning: HLA-HD file not found: ${hlahdFile}"
        return result
    }
    if (!optitypeFile.exists()) {
        println "Warning: OptiType file not found: ${optitypeFile}"
        return result
    }


    def optiAlleles = []
    optitypeFile.eachLine { line ->
        if (!line.startsWith("Allele") && line.trim()) {
            def cols = line.split("\t")
            cols.each { col ->
                col = col.trim()
                if (col) {
                    // Extract just the gene and first two numeric fields, ignore the rest
                    def matcher = col =~ /^(A|B|C)\*(\d{2,3}):(\d{2,3})(?::\d+)?$/
                    if (matcher.matches()) {
                        // extract HLA till second numeric value
                       def simplifiedAllele1 = "HLA-${matcher[0][1]}*${matcher[0][2]}:${matcher[0][3]}"
                       // Only add if no optiAlleles end with this simplified allele (or its last field)
                       optiAlleles << simplifiedAllele1
                    }
                }
            }
        } 
    }


    def hlahdAlleles = []
    hlahdFile.eachLine { line ->
        if (!line.startsWith("#") && line.trim()) {
            def cols = line.split("\t")
            cols.each { col ->
                col = col.trim()
                if (col) {
                    // Extract just the gene and first two numeric fields, ignore the rest
                    //def matcher = col =~ /^(DRB1|DQA1|DQB1|DPA1|DPB1)\*(\d{2,3}):(\d{2,3})(?::\d+)?$/
                    def matcher = col =~ /^HLA-(DRB1|DQA1|DQB1|DPA1|DPB1)\*(\d{2,3}):(\d{2,3})(?::\d+)?$/
                    if (matcher.matches()) {
                       def simplifiedAllele2 = "${matcher[0][1]}*${matcher[0][2]}:${matcher[0][3]}"
                       // Only add if no optiAlleles end with this simplified allele (or its last field)
                       hlahdAlleles << simplifiedAllele2
                    }
                }
            }
        }
    }


    // Combine ABC from OptiType + others from HLA-HD
    def combined = []
    combined.addAll(optiAlleles)
    combined.addAll(hlahdAlleles)

    // Remove duplicates by unique from groovy
    def uniqueCombined = combined.unique()

    // Join as comma-separated string
    def finalHLAString = uniqueCombined.join(",")

    // Return as map with patient id
    result[patientId] = finalHLAString
    return result
}



