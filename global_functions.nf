include { softwareVersionsToYAML } from "../subworkflows/nf-core/utils_nfcore_pipeline/"

// colors and other
def colorCodes() {
    def colorcodes = [:]
    // font
    colorcodes['bold']     = "\033[1m"
    colorcodes['reset']    = "\033[0m"
    colorcodes['line']     = "\033[0m"

    // colors
    colorcodes['yellow']     = "\033[0;33m"
    colorcodes['black']      = "\033[0;30m"
    colorcodes['red']        = "\033[0;31m"
    colorcodes['green']      = "\033[0;32m"
    colorcodes['green']      = "\033[0;32m"
    return(colorcodes)
}

// notify
def notify_slack(text, hook) {
    def myFile = file('./notify.json')
    myFile << '{"text": "'
    myFile << text.replace("\n",'\\n')
    myFile << '"}'
    println "curl -X POST -H 'Content-type: application/json' -d @./notify.json ${hook}".execute().text
    myFile.delete()

}

def makeSoftwareVersionYamlFile(ch_versions) {
        def ch_collated_versions = softwareVersionsToYAML(ch_versions)
        .collectFile(
            name: 'nf_core_'  + 'pipeline_software_' +  'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
         )

    return ch_collated_versions
}


// reverse complement DNA sequence
def revCompDNA(seq) {
	def newseq = seq.tr("ATGCatgc","TACGtacg").reverse()
	return(newseq)

}

def trim_NF_date(nfdate){
    newdate = nfdate.substring(0, nfdate.indexOf(".")).replaceAll("T", " ")
    return(newdate)
}

def final_message(title="") {
	def ostart = "${workflow.start}"
	def ostop = "${workflow.complete}"
	def start = trim_NF_date(ostart)
        def stop = trim_NF_date(ostop)
        def error = ""
        if (workflow.errorReport) {
            error = "\n```${workflow.errorReport}```\n"
        }

	def message =  "-"*51 + "\n"
	message = message + "*Pipeline ${title} completed!*".center(51) + "\n"
        message = message + "-"*51 + "\n"

	message = message + "- Launched by `$workflow.userName`" + "\n"
	message = message + "- Started at $start" + "\n"
    	message = message + "- Finished at $stop" + "\n"
    	message = message + "- Time elapsed: $workflow.duration" + "\n"
    	message = message + "- Execution status: ${ workflow.success ? 'OK' : 'failed' }" + "\n"
    	message = message + "```$workflow.commandLine```"+ "\n"
    	message = message + error + "-"*51 + "\n"
return (message)

}

// read a fasta file
def readFasta(fastapath) {
	def fastaMap = [:]
	def currentId = null
	def currentSeq = ""

	file(fastapath).readLines()
    .each {     
		if (it.startsWith(">")) {
			if (currentId) {
				fastaMap[currentId] = currentSeq
			}
			currentId = it.substring(1).trim()
			currentSeq = ""
		} else {
        	currentSeq = currentSeq + it.trim()
    	}
	}
	fastaMap[currentId] = currentSeq
		
	return(fastaMap)
}

// make named pipe
def unzipNamedPipe(filename) {
    def cmd = filename.toString()
    if (cmd[-3..-1] == ".gz") {
    	cmd = "<(zcat ${filename})"
    }
    return cmd
}

def zcatOrCat(filename) {
    def fname = filename.toString()
    def cmd = "cat ${filename}"
    if (fname[-3..-1] == ".gz") {
    	cmd = "zcat ${filename}"
    }
    return cmd
}


// evaluate input string for making SE or PE read channel with meta info
// suitable for nf-core modules
def fromStringToNFCoreSeqs(input_string, parseid = false) {
	def myseqs = channel.of()
	if (input_string.contains("{") && input_string.contains(",") && input_string.contains("}")) {
    	  myseqs = channel.fromFilePairs( input_string, checkIfExists: true ) 
    	  .map {[ [id: it[0], single_end:false ],  it[1] ] }
	} else {
          mypars = channel.fromFilePairs( input_string, size: 1, checkIfExists: true)
	  if (parseid) {
             myseqs = mypars.map {[ [id: it[1][0].name, single_end: true], it[1] ] }
	  } else {
     	     myseqs = mypars.map {[ [id: it[0], single_end:true],  it[1] ] }
	  }
    }
    return (myseqs)
}

// evaluate input string for making a file path or an empty map
def fromParToValueFileChannel(input_string) {
	myfile = []
	if (input_string) {
		myfile = channel.fromPath( input_string,  checkIfExists:true).first()
	}
    return (myfile)
}

// from nf-core standard channel to value channel 
def fromNFcoreToValueChannel(input_channel) {
     val_channel = input_channel.map{it[1]}.first()
     return (val_channel)
}

// evaluate input string for making a file path or an empty map
def addPrefixToFiles(input_nf_ch, prefix) {
	output_nf_ch = []
	output_nf_ch = input_nf_ch.map{ 
		meta, files ->  
    	[ [id: "${meta.id}${prefix}", single_end: meta.single_end], files ]  
	}
    return (output_nf_ch)
}

// remove the metamap from a channel 
def metaToCanonical(input_nf_ch, flat=false) {
    input_nf_ch.map { meta, files ->
        if (flat) {
            // produce [id, file1, file2, ...]
            [ meta.id ] + files.flatten()
        } else {
            // produce [id, [file1, file2, ...]]
            [ meta.id, files.flatten() ]
        }
    }
}

def flattenToTuple(input_ch) {
	tuples = input_ch.map{ 
		[ it[0], it[1..-1]  ]
	}
	return(tuples)
}


// remove the metamap from a channel 
def canonicalToMeta(input_ch) {
	nfcore_ch = input_ch.map { id, files ->
		if( files.size() == 2 ) {
			[ [ id: id, single_end : false], files ]

		} else {
			[ [ id: id, single_end : true], files ]
		}
	}
   return (nfcore_ch)
}

def subsetReads (reads, subset_num) {

    reads_types = reads.branch { meta, files ->
        pe: meta.single_end == false
        se: meta.single_end == true
    }

    flat_reads_pe = metaToCanonical(reads_types.pe, true)
    flat_reads_se = metaToCanonical(reads_types.se, true)

    splitted_pe = flat_reads_pe.splitFastq( by: subset_num, decompress:true, file:true, limit:subset_num, pe: true )
    splitted_se = flat_reads_se.splitFastq( by: subset_num, decompress:true, file:true, limit:subset_num )
	splitted = splitted_pe.mix(splitted_se)
	tuples = flattenToTuple(splitted)
	
	return (canonicalToMeta(tuples))


}



