// make named pipe 
def unzipNamedPipe(filename) { 
    def cmd = filename.toString()
    if (cmd[-3..-1] == ".gz") {
    	cmd = "<(zcat ${filename})"
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



// unzip command 
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

