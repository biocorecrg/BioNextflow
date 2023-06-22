// notify

def notify_slack(text, hook) {
    def myFile = file('./notify.json')
    myFile << '{"text": "'
    myFile << text.replace("\n",'\\n')
    myFile << '"}'
    println "curl -X POST -H 'Content-type: application/json' -d @./notify.json ${hook}".execute().text
    myFile.delete()

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
        if ("${ workflow.success}") {
            error = "\n```${workflow.errorReport}```\n"
        }

	def message =  "-"*51 + "\n"
	message = message + "*Pipeline ${title} completed!*".center(51) + "\n"
    message = message + "- Started at $start" + "\n"
    message = message + "- Finished at $stop" + "\n"
    message = message + "- Time elapsed: $workflow.duration" + "\n"
    message = message + "- Execution status: ${ workflow.success ? 'OK' : 'failed' }" + "\n"
    message = message + "```$workflow.commandLine```"+ "\n"
    message = message + error + "-"*51 + "\n"
/*
---------------------------------------------------------------------------------
                        *Pipeline ${title} completed!*
```$workflow.commandLine```
- Launched by `$workflow.userName`
- Started at $start
- Finished at $stop
- Time elapsed: $workflow.duration 
- Execution status: ${ workflow.success ? 'OK' : 'failed' } ${error}
---------------------------------------------------------------------------------
*/

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
        pe: it[1][[1]]
        se: !it[1][[1]]
    }
    return (result)
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
