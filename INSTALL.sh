#!/usr/bin/env sh

# Check pipeline dependency
bionext_ver="0.7.2"

# Check presence of a lib
if [ -d "lib" ]; then {
	echo "BioNextflow installed"
} else {
  	echo "BioNextflow not installed. Installing it..."
	wget https://github.com/biocorecrg/BioNextflow/archive/${bionext_ver}.tar.gz
	tar -zvxf ${bionext_ver}.tar.gz
	mv "BioNextflow-${bionext_ver}/lib" .
	rm ${bionext_ver}.tar.gz
	rm -fr BioNextflow-${bionext_ver}
}
fi



