#!/usr/bin/bash -l

mkdir db
pushd db
URL=ftp://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Fungi
fungi.ITS.fna.gz
for db in fungi.ITS.fna.gz fungi.ITS.gbff.gz fungi.18SrRNA.fna.gz fungi.18SrRNA.gbff.gz fungi.28SrRNA.fna.gz fungi.28SrRNA.gbff.gz
do
    curl -O $URL/$db
done

popd