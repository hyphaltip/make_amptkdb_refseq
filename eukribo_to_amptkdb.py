#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Create a 18S DB for amptk from EukRibo"""

import os
import gzip
import re
import csv
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
    seqs = {}
    acc2org = {}
    with gzip.open(taxofas, "rt", encoding='latin-1') as handle, \
        open("eukribo.fasta", "w") as ofile:
        for seq_record in SeqIO.parse(handle, "fasta"):
                    accession = seq_record.id
                    seqs[accession] = seq_record
                    desc = seq_record.description
                    taxonlst = desc.split('|')
                    lookup = taxonlst[-1]
                    while taxonlst:
                        if lookup.startswith('isolate') or lookup.startswith('strain'):
                            lookup = taxonlst.pop()
                        else:
                            lookup = lookup.replace("+", " ")
                            break
                    acc2org[accession] = lookup

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
            for t in taxonstr.split(","):
                if len(t) == 0:
                    continue
                (k, v) = t.split(":")
                if len(v) == 0:
                    taxonstr = taxonstr.replace(f'{k}:,', '')
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