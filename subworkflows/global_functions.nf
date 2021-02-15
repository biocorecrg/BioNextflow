// make named pipe 
def unzipNamedPipe(filename) { 
    def cmd = filename.toString()
    if (cmd[-3..-1] == ".gz") {
    	cmd = "<(zcat ${filename})"
    }
    return cmd
}

// unzip command 
def unzipCmd(filename, unzippedname, copy="") { 
    def cmd = "ln -s ${filename} ${unzippedname}"
	if (copy!="") {
		cmd = "cp ${filename} ${unzippedname}"
	}
    def ext = filename.getExtension()
    if (ext == "gz") {
    	cmd = "zcat ${filename} > ${unzippedname}"
    }
    return cmd
}

