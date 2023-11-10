#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Create a 18S DB for amptk from EukRibo"""

import argparse
import os
import time
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

def seq_to_ids_single(seq_dict,limit=0,verbose=False):
    acc2id = {}
    id2acc = {}
    i = 0
    for acc in seq_dict:
        seq_record = seq_dict[acc]
        desc = seq_record.description
        acc = seq_record.id
        taxonlst = desc.split('|')
        lookup = taxonlst.pop()
        # this searches for IDs one at time which is slower but necessary when missing?
 #       tm_pylookup_start = time.time()
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
                        print(f"cannot find taxid for {name} {seq_record.id} {desc} ... will go up a level")
                        lookup = taxonlst.pop()
                    else:
                        acc2id[acc] = [str(taxid), '']
                        id2acc[str(taxid)] = acc
                        print(f"storing acc: {acc} to {taxid} and {taxid} -> {acc} based on {lookup} -> {desc}")
                        keep = True
            if keep:  # break out of loop if we have found an ID
                break
#        tm_pylookup_end = time.time()
        if i > 0 and i % 1000 == 0:
            print(f"processed {i} sequences")
        i += 1
        if i > limit and limit > 0:
            break
    return (acc2id,id2acc)


def seq_to_ids_combined(seq_dict,limit=0,verbose=False):
    acc2id = {}
    id2acc = {}
    i = 0
    for acc in seq_dict:
        seq_record = seq_dict[acc]
        desc = seq_record.description
        acc = seq_record.id
        taxonlst = desc.split('|')
        lookup = taxonlst.pop()
        # this searches for IDs one at time which is slower but necessary when missing?
 #       tm_pylookup_start = time.time()
 
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
                        print(f"cannot find taxid for {name} {seq_record.id} {desc} ... will go up a level")
                        lookup = taxonlst.pop()
                    else:
                        acc2id[acc] = [str(taxid), '']
                        id2acc[str(taxid)] = acc
                        print(f"storing acc: {acc} to {taxid} and {taxid} -> {acc} based on {lookup} -> {desc}")
                        keep = True
            if keep:  # break out of loop if we have found an ID
                break
#        tm_pylookup_end = time.time()
        if i > 0 and i % 1000 == 0:
            print(f"processed {i} sequences")
        i += 1
        if i > limit and limit > 0:
            break
    return (acc2id,id2acc)

def main(limit=0,verbose=False):
    taxotable = os.path.basename(eukribotsv).replace("?download=1","")
    taxofas = os.path.basename(eukribofas).replace("?download=1","")
    if not os.path.exists(taxotable):
        download_file(eukribotsv, taxotable)
    if not os.path.exists(taxofas):
        download_file(eukribofas, taxofas)
    with gzip.open(taxofas, "rt", encoding='latin-1') as handle, \
        open("eukribo.fasta", "w") as ofile:
        seqs = {}
        tm_parse_start = time.time()
        for seq_record in SeqIO.parse(handle, "fasta"):
            accession = seq_record.id
            seqs[accession] = seq_record
        tm_parse_end = time.time()
        if verbose:
            print(f"done reading in sequences that took {tm_parse_end - tm_parse_start} seconds")

        tm_lookup_start = time.time()
        (acc2id,id2acc) = seq_to_ids_single(seqs,limit,verbose)
        tm_lookup_end = time.time()
        idnum = len(acc2id)
        if verbose:
            print(f"Looking up {idnum} IDs took {tm_lookup_end - tm_lookup_start} seconds")
        ids = list(id2acc.keys())
        cmd = ["taxonkit", 'reformat', '-I', '1', '-F', '-f',
               "'k:{K},p:{p},c:{c},o:{o},f:{f},g:{g},s:{s}'"]

        p = Popen(cmd, stdout=PIPE, stdin=PIPE, stderr=PIPE)
        idsstr = "\n".join(ids) + "\n"
        resdata = p.communicate(input=idsstr.encode())[0]
        for line in resdata.decode().splitlines():
            (taxonid, taxonstr) = line.split("\t")
            taxonstr = taxonstr.replace('\'', '')
            taxonstr = re.sub(r',$', '', taxonstr)
            if verbose:
                print(f'taxonid is {taxonid} and strmatch is {taxonstr}')
            # drop empty taxon levels
            taxonstr = re.sub("([kpcofgs]):,",'',taxonstr)
            if taxonid in id2acc:
                acc = id2acc[taxonid]
                acc2id[acc][1] = taxonstr  # store the taxstring associated with acccession and taxonid
            else:
                print(f'Taxonid {taxonid} is not found already in acc2id')
                break
        for acc in seqs:
            if acc in acc2id:
                (taxid,taxstring) = acc2id[acc]
            else:
                taxstring = ''
            seq = seqs[acc]
            ofile.write(f">{acc};tax={taxstring}\n")
            ofile.write(softwrap(str(seq.seq))+"\n")

parser = argparse.ArgumentParser(description='Process EukRibo 18S db')

# Add arguments
parser.add_argument('--limit', type=int, default=0, required=False, help='specify a test run of this many lookups')
parser.add_argument('-v', '--verbose', action='store_true')  
# Parse the arguments
args = parser.parse_args()

if __name__ == '__main__':
    main(args.limit,args.verbose)