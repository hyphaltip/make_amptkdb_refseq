#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Create a 18S DB for amptk from EukRibo"""

import argparse
import gzip
import os
import re
import time
from subprocess import PIPE, Popen

import pandas as pd
import pytaxonkit
import requests
from Bio import SeqIO

eukribotsv = "https://zenodo.org/records/6896896/files/46346_EukRibo-02_2022-07-22.tsv.gz?download=1"
eukribofas = "https://zenodo.org/records/6896896/files/46346_EukRibo-02_full_seqs_2022-07-22.fas.gz?download=1"


def download_file(url, filename):
    response = requests.get(url, timeout=30)
    with open(filename, "wb") as f:
        f.write(response.content)


def softwrap(string, every=80):
    lines = []
    for i in range(0, len(string), every):
        lines.append(string[i : i + every])
    return "\n".join(lines)


def seq_to_ids_single(seq_dict, limit=0, verbose=False):
    acc2id = {}
    id2acc = {}
    i = 0
    for acc in seq_dict:
        seq_record = seq_dict[acc]
        desc = seq_record.description
        acc = seq_record.id
        taxonlst = desc.split("|")
        lookup = taxonlst.pop()
        # this searches for IDs one at time which is slower but necessary when missing?
        #       tm_pylookup_start = time.time()
        while lookup:
            if (
                lookup.startswith("isolate")
                or lookup.startswith("strain")
                or lookup.startswith("clone")
            ):
                lookup = taxonlst.pop()
            else:
                lookup = lookup.replace("+", " ")
                lookup = re.sub(r"^[kpcofgs]:", "", lookup)
                # this is slower as we are reading in one at a time
                ptx = pytaxonkit.name2taxid([lookup])
                keep = False
                for idx, r in ptx.iterrows():
                    name = r["Name"]
                    taxid = r["TaxID"]
                    if pd.isna(taxid):
                        if verbose:
                            print(
                                f"cannot find taxid for {name} {seq_record.id} {desc} ... will go up a level"
                            )
                        lookup = taxonlst.pop()
                    else:
                        str_taxid = str(taxid)
                        acc2id[acc] = [str_taxid, ""]
                        if str_taxid not in id2acc:
                            id2acc[str_taxid] = set()
                        id2acc[str(taxid)].add(acc)
                        if verbose:
                            print(
                                f"storing acc: {acc} to {taxid} and {taxid} -> {acc} based on {lookup} -> {desc}"
                            )
                        keep = True
            if keep:  # break out of loop if we have found an ID
                break
        #        tm_pylookup_end = time.time()
        if i > 0 and i % 1000 == 0:
            if verbose:
                print(f"processed {i} sequences")
        i += 1
        if i > limit and limit > 0:
            break
    return (acc2id, id2acc)


def lookup_set(lookup_needed, verbose=False):
    query_order = set()
    name2taxid = {}
    for acc in lookup_needed:
        if len(lookup_needed[acc]['taxonid']) == 0:
            # need to run this through
            # take last name
            if len(lookup_needed[acc]['taxonq']):
                lookup = lookup_needed[acc]['taxonq'].pop()
            else:
                continue
            # if this is a strain, isolate, clone, etc. then take the next name
            if (
                lookup.startswith("isolate")
                or lookup.startswith("strain")
                or lookup.startswith("clone")
            ):
                lookup = lookup_needed[acc]['taxonq'].pop()
            # replace + with space
            lookup = lookup.replace("+", " ")
            # replace taxon level prefix that exists in some of these
            lookup = re.sub(r"^[kpcofgs]:", "", lookup)
            lookup_needed[acc]['currentsearch'] = lookup
        else:
            lookup = lookup_needed[acc]['currentsearch']
            name2taxid[lookup] = lookup_needed[acc]['taxonid']
    # if we know the right answer, eg that this name has aleady been assigned a taxid, let's reuse it
    for acc in lookup_needed:
        lookup = lookup_needed[acc]['currentsearch']
        if lookup in name2taxid:
            lookup_needed[acc]['taxonid'] = name2taxid[lookup]
        else:
            query_order.add(lookup)

    updated_number = 0
    if len(query_order) > 0:
        ptx = pytaxonkit.name2taxid(query_order)
        for idx, r in ptx.iterrows():
            name = r["Name"]
            taxid = r["TaxID"]
            if pd.isna(taxid):
                print(f"cannot find taxid for {name} ... will go up a level")
            else:
                if name not in name2taxid:
                    name2taxid[name] = str(taxid)
        # now we will try to assign these new taxids to the accs that need them
        for acc in lookup_needed:
            if len(lookup_needed[acc]['taxonid']) == 0:  # we don't have a taxid yet
                if lookup_needed[acc]['currentsearch'] in name2taxid:
                    # assign value to column 1 for the lookup name2taxid from the search term used
                    lookup_needed[acc]['taxonid'] = name2taxid[lookup_needed[acc]['currentsearch']]
                    updated_number += 1
                else:
                    print(f"no hit found for {lookup_needed[acc]['currentsearch']} still need to search in {lookup_needed[acc]['taxonq']}")
    # this might need another iteration but that will just mean calling this function again
    return (lookup_needed, len(query_order), updated_number)


def seq_to_ids_multi(seq_dict, max_iter=7, limit=0, verbose=False):
    i = 0
    if verbose:
        print(f'there are {len(seq_dict)} sequences')
    lookups_to_do = {}
    for acc in seq_dict:
        if i > limit and limit > 0:
            break
        seq_record = seq_dict[acc]
        (idacc,desc) = seq_record.description.split(" ",1)
        acc = seq_record.id
        taxonlst = desc.split("|")
        # col 0 is taxonlist which will be updated as we test
        # col 1 will be the taxid when this works
        # col 2 will be the current taxon string we are searching with (as we pop off col 0)
        lookups_to_do[acc] = { 'taxonq': taxonlst, 'taxonid': '', 'currentsearch': ''}
        i += 1

    updated_number = -1
    iters = 0
    while updated_number != 0:
        if iters > max_iter:
            break
        tm_pylookup_start = time.time()
        (lookups_to_do, querycount, updated_number) = lookup_set(lookups_to_do, verbose)
        tm_pylookup_end = time.time()
        if verbose:
            print(f"ran lookup, asked for {querycount} out of ({len(lookups_to_do)}) and did {updated_number} updates in {tm_pylookup_end - tm_pylookup_start} seconds")
        iters += 1

    lastlookacc2id = {}
    lastlookid2acc = {}
    for acc in lookups_to_do:
        id = lookups_to_do[acc]['taxonid']
        lastlookacc2id[acc] = [id, ""]
        if id not in lastlookid2acc:
            lastlookid2acc[id] = set()
        lastlookid2acc[id].add(acc)
    return (lastlookacc2id, lastlookid2acc)


def main(method="simple", limit=0, verbose=False):
    taxotable = os.path.basename(eukribotsv).replace("?download=1", "")
    taxofas = os.path.basename(eukribofas).replace("?download=1", "")
    if not os.path.exists(taxotable):
        download_file(eukribotsv, taxotable)
    if not os.path.exists(taxofas):
        download_file(eukribofas, taxofas)
    with gzip.open(taxofas, "rt", encoding="latin-1") as handle, \
        open("eukribo.fasta", "w") as ofile:
        seqs = {}
        tm_parse_start = time.time()
        for seq_record in SeqIO.parse(handle, "fasta"):
            accession = seq_record.id
            seqs[accession] = seq_record
        tm_parse_end = time.time()
        if verbose:
            print(f"done reading in sequences that took {tm_parse_end - tm_parse_start} seconds")

        (tm_lookup_start, tm_lookup_end) = (0, 0)
        acc2id = {}
        id2acc = {}
        if method == "simple":
            tm_lookup_start = time.time()
            (acc2id, id2acc) = seq_to_ids_single(seqs, limit, verbose)
            tm_lookup_end = time.time()
        elif method == "multi":
            tm_lookup_start = time.time()
            (acc2id, id2acc) = seq_to_ids_multi(seqs, 7, limit, verbose)
            tm_lookup_end = time.time()
        else:
            print(f"unknown method {method} specified")
            return

        idnum = len(acc2id)
        if verbose:
            print(
                f"Looking up {idnum} IDs with {method} took {tm_lookup_end - tm_lookup_start} seconds"
            )
        ids = list(id2acc.keys())
        cmd = [
            "taxonkit",
            "reformat",
            "-I",
            "1",
            "-F",
            "-f",
            "'k:{K},p:{p},c:{c},o:{o},f:{f},g:{g},s:{s}'",
        ]

        p = Popen(cmd, stdout=PIPE, stdin=PIPE, stderr=PIPE)
        idsstr = "\n".join(ids) + "\n"
        resdata = p.communicate(input=idsstr.encode())[0]
        for line in resdata.decode().splitlines():
            (taxonid, taxonstr) = line.split("\t")
            taxonstr = taxonstr.replace("'", "")
            taxonstr = re.sub(r",$", "", taxonstr)
            if verbose:
                print(f"taxonid is {taxonid} and strmatch is {taxonstr}")
            # drop empty taxon levels
            taxonstr = re.sub("([kpcofgs]):,", "", taxonstr)
            if taxonid in id2acc:
                for acc in id2acc[taxonid]:
                    acc2id[acc][1] = taxonstr  # store the taxstring associated with acccession and taxonid
            else:
                break
        for acc in seqs:
            taxstring = ""
            if acc in acc2id:
                (taxid, taxstring) = acc2id[acc]
            seq = seqs[acc]
            ofile.write(f">{acc};tax={taxstring}\n")
            ofile.write(softwrap(str(seq.seq)) + "\n")


parser = argparse.ArgumentParser(description="Process EukRibo 18S db")

# Add arguments
parser.add_argument(
    "--limit",
    type=int,
    default=0,
    required=False,
    help="specify a test run of this many lookups",
)
parser.add_argument("-v", "--verbose", action="store_true")
parser.add_argument(
    "-m",
    "--method",
    default="simple",
    choices=["simple", "multi"],
    help="single or multi method for lookup speed",
)
# Parse the arguments
args = parser.parse_args()
if __name__ == "__main__":
    main(args.method, args.limit, args.verbose)
