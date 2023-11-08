#!/usr/bin/env python3

import sys
import os
import gzip
import subprocess
from subprocess import Popen, PIPE, STDOUT
from Bio import SeqIO
import pytaxonkit

def main(indb="db"):
    for filename in os.listdir(indb):
        if filename.endswith(".gbff.gz"):
            ids_to_lookup = {}
            seqs = {}
            with gzip.open(os.path.join(indb,filename), "rt") as handle:
                for seq_record in SeqIO.parse(handle, "genbank"):
                    accession = seq_record.id                    
                    seqs[accession] = seq_record.seq
                    # print(seq_record.annotations['organism'])
                    ids_to_lookup[seq_record.annotations['organism']] = accession
                    break
                print('ids are',list(ids_to_lookup.keys()))
                ids = pytaxonkit.name2taxid(list(ids_to_lookup.keys()))
                print(ids)
                cmd = "taxonkit reformat -I 1 -f 'k:{k};p:{p};c:{c};o:{o};f:{f};g:{g};s:{s}'"

                p = Popen(cmd, stdout=PIPE, stdin=PIPE, stderr=PIPE)
                stdout_data = p.communicate(input="\n".join(ids))[0]
                print(stdout_data)

if __name__ == '__main__':
    main()
