#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import gzip
import re
import subprocess
import logging
from subprocess import Popen, PIPE, STDOUT
from Bio import SeqIO
import pytaxonkit

def softwrap(string, every=80):
    lines = []
    for i in range(0, len(string), every):
        lines.append(string[i:i+every])
    return '\n'.join(lines)

def main(indb="db"):
    for filename in os.listdir(indb):
        if filename.endswith(".gbff.gz"):
            seqs = {}
            acc2org = {}
            ofilename = filename.replace(".gbff.gz",".fasta")
            with gzip.open(os.path.join(indb,filename), "rt") as handle, open(ofilename,"w") as ofile:
                for seq_record in SeqIO.parse(handle, "genbank"):
                    accession = seq_record.id                    
                    seqs[accession] = seq_record.seq
                    # print(seq_record.annotations['organism'])                    
                    acc2org[accession] = seq_record.annotations['organism']
                #logging.info(f'species names are {"\n".join(list(set(acc2org.values())))}')
                ptx = pytaxonkit.name2taxid(list(set(acc2org.values())))
                #logging.info(ptx)
                orgs2ids = {}
                ids2orgs = {}
                ids = set()
                for idx,r in ptx.iterrows():
                    name = r['Name']
                    taxid = r['TaxID']
                    orgs2ids[name] = { 'id': str(taxid), 'taxstring': '' }
                    ids2orgs[str(taxid)] = { 'name': name, 'taxstring': '' }
                    ids.add(str(taxid))
                # logging.info(f'the ids are {"\n".join(ids)}')
                # amptk formatting string
                cmd = ["taxonkit",'reformat','-I','1', '-f', "'k:{K},p:{p},c:{c},o:{o},f:{f},g:{g},s:{s}'"]

                p = Popen(cmd, stdout=PIPE, stdin=PIPE, stderr=PIPE)
                idsstr = "\n".join(ids) + "\n"              
                resdata = p.communicate(input=idsstr.encode())[0]
                for line in resdata.decode().splitlines():
                    (taxonid,taxonstr) = line.split("\t")
                    taxonstr = taxonstr.replace('\'','')
                    taxonstr = re.sub(r',$','',taxonstr)
                    print(taxonstr)
                    # drop empty taxon levels
                    for t in taxonstr.split(","):
                        (k,v) = t.split(":")
                        if len(v) == 0:
                            taxonstr = taxonstr.replace(f',{k}:,','')
                    if taxonid in ids2orgs:
                        ids2orgs[taxonid]['taxstring'] = taxonstr
                        orgs2ids[ids2orgs[taxonid]['name']]['taxstring'] = taxonstr
                for acc in seqs:
                    org = acc2org[acc]
#                    taxid = orgs2ids[org]['id']
                    taxstring = orgs2ids[org]['taxstring']
                    seq = seqs[acc]
                    ofile.write(f">{acc};tax={taxstring}\n")
                    ofile.write(softwrap(str(seq))+"\n")

if __name__ == '__main__':
    main()
