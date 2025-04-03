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

def separateSEandPE(fastq) {

    def result = fastq.branch {
    	pe: it[1].size() > 1  // Checks if there are two reads (paired-end)
    	se: it[1].size() == 1  // If only one read, it's single-end
    }
    return (result)
}

def checkSEorPE(fastq){
   num = fastq.map{it[1]}.first().size()
   res = num.map{   
	   it -> ( it == 2 ? "PE" : (it == 1 ? "SE" : "ERR" ))
   }
   return(res)
}

/*
* unzip command: how to use it
* script:
*  def unzip      = unzipCmd(file)
*  def file_name  = unzip[0]
*  def cmd_unzip  = unzip[1]
*  def cmd_clean  = unzip[2]
*  """
*  ${cmd_unzip}
*  do something with ${file_name}
*  ${cmd_clean}
*  """
*/

def unzipCmd(filename) {
    def ext = filename.getExtension()
    def fname = filename
    def cmd = ""
    def clean = ""
    if (ext == "gz") {
    	fname = filename.baseName
    	cmd = "zcat ${filename} > ${fname}"
        clean = "rm ${fname}"
    }
    return [fname, cmd, clean]
}
