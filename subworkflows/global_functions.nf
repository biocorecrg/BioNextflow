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
        pe: it[1][1]
        se: !it[1][1]
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
