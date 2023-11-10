#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Create a 18S DB for amptk from EukRibo"""

import os
import gzip
import re
import pandas as pd
from subprocess import Popen, PIPE
from Bio import SeqIO
import pytaxonkit
import requests

eukribotsv = "https://zenodo.org/records/6896896/files/46346_EukRibo-02_2022-07-22.tsv.gz?download=1"
eukribofas = "https://zenodo.org/records/6896896/files/46346_EukRibo-02_full_seqs_2022-07-22.fas.gz?download=1"
def download_file(url, filename):
    response = requests.get(url,timeout=30)
    with open(filename, 'wb') as f:
        f.write(response.content)



def softwrap(string, every=80):
    lines = []
    for i in range(0, len(string), every):
        lines.append(string[i:i+every])
    return '\n'.join(lines)


def main():
    taxotable = os.path.basename(eukribotsv).replace("?download=1","")
    taxofas = os.path.basename(eukribofas).replace("?download=1","")
    if not os.path.exists(taxotable):
        download_file(eukribotsv, taxotable)
    if not os.path.exists(taxofas):
        download_file(eukribofas, taxofas)
    with gzip.open(taxofas, "rt", encoding='latin-1') as handle, \
        open("eukribo.fasta", "w") as ofile:
        seqs = {}
        org2id = {}
        id2org = {}
        ids = set()
        i = 0
        for seq_record in SeqIO.parse(handle, "fasta"):
            accession = seq_record.id
            seqs[accession] = seq_record
            desc = seq_record.description
            taxonlst = desc.split('|')
            lookup = taxonlst.pop()
            # this searches for IDs one at time which is slower but necessary when missing? 
            while lookup:
                if lookup.startswith('isolate') or lookup.startswith('strain') or lookup.startswith("clone"):
                    lookup = taxonlst.pop()
                else:
                    lookup = lookup.replace("+", " ")
                    lookup = re.sub(r'^[kpcofgs]:','',lookup)
# this is slower as we are reading in one at a time
                    ptx = pytaxonkit.name2taxid([lookup])
                    keep = False
                    for idx, r in ptx.iterrows():
                        name = r['Name']
                        taxid = r['TaxID']
                        if pd.isna(taxid):
                            print(f"cannot find taxid for {name} {accession} {desc} ... will go up a level")
                            lookup = taxonlst.pop()
                        else:
                            org2id[name] = {'id': str(taxid), 'taxstring': ''}
                            id2org[str(taxid)] = {'name': name, 'taxstring': ''}
                            print(f"storing {name} {taxid} based on {lookup} -> {desc}")
                            ids.add(str(taxid))
                            keep = True
                if keep:  # break out of loop if we have found an ID
                    break
            if i > 0 and i % 1000 == 0:
                print(f"processed {i} sequences")
            i += 1
        print(sorted(ids))
        cmd = ["taxonkit", 'reformat', '-I', '1', '-F', '-f',
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
#            if not m:
#                print(taxonstr)
#            for t in taxonstr.split(","):
#                if len(t) == 0:
#                    continue
#                # remove special case of o:Penicillaria (in: tube anenomes)
#                t = re.sub(r'\(in:.+\)','',t)
#                if t.count(':') != 1:
#                    print(t)
#                (k, v) = t.split(":")
#                if len(v) == 0:
#                    taxonstr = taxonstr.replace(f'{k}:,', '')
            if taxonid in id2org:
                id2org[taxonid]['taxstring'] = taxonstr
                org2id[id2org[taxonid]['name']]['taxstring'] = taxonstr
        for acc in seqs:
            org = acc2org[acc]
            taxstring = org2id[org]['taxstring']
            seq = seqs[acc]
            ofile.write(f">{acc};tax={taxstring}\n")
            ofile.write(softwrap(str(seq.seq))+"\n")


if __name__ == '__main__':
    main()

# import csv
#        tblin = csv.reader(tblfh, delimiter="\t")
#        header = next(tblin)
#        print(header)
#        names = set()
#        for row in tblin:
#            print(row)
#            acc = row[0]
#            supergroup = row[1]
#            UniEuk_taxonomy_string = row[4]
#            print(acc,UniEuk_taxonomy_string)