#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import gzip
import re
import requests
from subprocess import Popen, PIPE
from Bio import SeqIO
import pytaxonkit

fungi_url = 'https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Fungi'
fungi_files = ('fungi.ITS.fna.gz', 'fungi.ITS.gbff.gz',
               'fungi.18SrRNA.fna.gz', 'fungi.18SrRNA.gbff.gz',
               'fungi.28SrRNA.fna.gz', 'fungi.28SrRNA.gbff.gz')

def download_file(url, filename):
    response = requests.get(url,timeout=30)
    with open(filename, 'wb') as f:
        f.write(response.content)


def softwrap(string, every=80):
    lines = []
    for i in range(0, len(string), every):
        lines.append(string[i:i+every])
    return '\n'.join(lines)


def main(indb="db"):
    if not os.path.exists(indb):
        os.mkdir(indb)
    for fname in fungi_files:
        localfile = os.path.join(indb,fname)
        if not os.path.exists(localfile):
            download_file(os.path.join(fungi_url,fname), localfile)
    for filename in os.listdir(indb):
        if filename.endswith(".gbff.gz"):
            seqs = {}
            acc2org = {}
            ofilename = filename.replace(".gbff.gz", ".fasta")
            ifilename = os.path.join(indb, filename)
            with gzip.open(ifilename, "rt") as handle, \
                open(ofilename, "w") as ofile:
                for seq_record in SeqIO.parse(handle, "genbank"):
                    accession = seq_record.id                    
                    seqs[accession] = seq_record.seq
                    # print(seq_record.annotations['organism'])
                    acc2org[accession] = seq_record.annotations['organism']
                # print('species names are',
                # "\n".join(list(set(acc2org.values()))))
                ptx = pytaxonkit.name2taxid(list(set(acc2org.values())))
                org2id = {}
                id2org = {}
                ids = set()
                for idx, r in ptx.iterrows():
                    name = r['Name']
                    taxid = r['TaxID']
                    org2id[name] = {'id': str(taxid), 'taxstring': ''}
                    id2org[str(taxid)] = {'name': name, 'taxstring': ''}
                    ids.add(str(taxid))
                # amptk formatting string
                cmd = ["taxonkit", 'reformat', '-I', '1', '-f',
                       "'k:{K},p:{p},c:{c},o:{o},f:{f},g:{g},s:{s}'"]

                p = Popen(cmd, stdout=PIPE, stdin=PIPE, stderr=PIPE)
                idsstr = "\n".join(ids) + "\n"              
                resdata = p.communicate(input=idsstr.encode())[0]
                for line in resdata.decode().splitlines():
                    (taxonid, taxonstr) = line.split("\t")
                    taxonstr = taxonstr.replace('\'', '')
                    taxonstr = re.sub(r',$', '', taxonstr)
                    print(taxonstr)
                    # drop empty taxon levels
                    taxonstr = re.sub("([kpcofgs]):,",'',taxonstr)
                    if taxonid in id2org:
                        id2org[taxonid]['taxstring'] = taxonstr
                        org2id[id2org[taxonid]['name']]['taxstring'] = taxonstr
                for acc in seqs:
                    org = acc2org[acc]
                    taxstring = org2id[org]['taxstring']
                    seq = seqs[acc]
                    ofile.write(f">{acc};tax={taxstring}\n")
                    ofile.write(softwrap(str(seq))+"\n")


if __name__ == '__main__':
    main()
