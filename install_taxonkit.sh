#!/usr/bin/bash -l
# this setups folder for taxonkit runs
mkdir -p taxdump
pushd taxdump
curl -C - -O ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -zxvf taxdump.tar.gz

mkdir -p $HOME/.taxonkit
rsync -a names.dmp nodes.dmp delnodes.dmp merged.dmp $HOME/.taxonkit

rm *dmp readme.txt gc.prt
