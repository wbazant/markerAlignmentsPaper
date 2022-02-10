import sys
import re
import os
import argparse
import pysam
import re
import json
import statistics
import logging
import subprocess
import sqlite3

from marker_alignments.store import SqliteStore
from ete3 import NCBITaxa

from pathlib import Path

def read_marker_to_taxon(path):
    result = {}
    with open(path, 'r') as f:
        for line in f:
            (marker, taxon) = line.rstrip().split("\t")
            result[marker] = taxon
    return result

# removed $ from the end compared to marker_alignments version
#other_MGCollapse_EPSP3-12_382_886_0:0:0_0:0:0_0
eukprot_refdb_regex_taxon = "^[a-z]+-(.*_[a-z-]+(?:_.*)?)-[0-9]+[ab]t2759-[A-Z]\d|^[a-z]+-(.*)...Collapse_[^_]*|^()other_MGCollapse_[^_]*"
pattern_taxon = re.compile(eukprot_refdb_regex_taxon)

def next_g(search):
    return next(g for g in search.groups() if g is not None)


def read_scrambled_name_to_taxid_from_marker_to_taxon(path):
    result = {}
    with open(path, 'r') as f:
        for line in f:
            if line.find("Collapse") > -1:
                continue
            (marker, taxon_str) = line.rstrip().split("\t")
            taxon = int(taxon_str)
            search = pattern_taxon.search(marker)
            scrambled_name = next_g(search)
            if scrambled_name in result and result[scrambled_name] != taxon:
                raise ValueError(line, taxon, result[scrambled_name])
            result[scrambled_name] = taxon

    return result

trad_ranks = {"superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"}

def base_taxid(ncbi, taxids):
    tree = ncbi.get_topology(taxids)
    r = tree.get_tree_root()
    return str(r.name)

def get_match_type(ncbi, taxids):
    tree = ncbi.get_topology(taxids)

    r = tree.get_tree_root()
    if hasattr(r, 'rank') and r.rank in trad_ranks:
        result = r.rank
    else:
        lineage = ncbi.get_lineage(r.name)
        ranks = ncbi.get_rank(lineage)
        ranks_increasing = [ranks[x] for x in reversed(lineage) if ranks[x] in trad_ranks]

        if not ranks_increasing:
            raise ValueError(lineage, ranks)
        result = ranks_increasing[0]
    return result

def ct(n, rank):
    if n == 1 and rank == "genus":
        return "One result, correct genus"
    if n > 1 and rank == "genus":
        return "Many results, correct genus"
    if n == 1 and rank != "genus":
        return "One result, incorrect genus"
    if n > 1 and rank != "genus":
        return "Many results, incorrect genus"

header = [
"No results",
"One result, correct genus" ,
"Many results, correct genus",
"One result, incorrect genus",
"Many results, incorrect genus",
"Signal detected",
"Detected signal is one species"
]
header_shortcuts = {
        "No results": "NR",
        "One result, correct genus": "OC" ,
        "Many results, correct genus": "MC",
        "One result, incorrect genus": "OI",
        "Many results, incorrect genus": "MI",
}

def add_stats(counts):
    s = sum(counts.values())
    counts["Signal detected"] = round(1.0 * (s - counts["No results"] ) / s, 3)
    counts["Detected signal is one species"] = round(1.0 * (counts["One result, correct genus"] + counts["One result, incorrect genus"]) / (s - counts["No results"]), 3)

def sc(all_results, input_file, taxid):
    if taxid not in all_results[input_file]:
        return "XX"
    return header_shortcuts[all_results[input_file][taxid]]

# it could start with a ?
# if yes, it could be a list of comma-separated results
# it could have taxid at the beginning
def get_taxid_from_string(ncbi, x):
    if x[0] == '?':
        ids = [get_taxid_from_string(ncbi, _x) for _x in x[1:].split(",")]
        return base_taxid(ncbi, ids)

    if len(x.split("|")) > 1:
        return x.split("|")[0]

    name2taxid = ncbi.get_name_translator([x])
    vs = [x for xx in name2taxid.values() for x in xx]
    return vs[0]

def do_one_line(scrambled_name_to_taxid, ncbi, xs):
    source_taxon = xs.pop(0)
    source_taxid = scrambled_name_to_taxid[source_taxon]
    n = int(xs.pop(0))
    if n:
        vs = [get_taxid_from_string(ncbi, x) for x in xs]
        vs.append(str(source_taxid))
        rank = get_match_type(ncbi, vs)
        return source_taxid, ct(n, rank)
    else:
        return source_taxid, "No results"

def do_one(scrambled_name_to_taxid, ncbi, input_file):
    counts = {h:0 for h in header}
    cs_for_species = {}
    with open(input_file, 'r') as f:
        for l in f:
            xs = l.rstrip().split("\t")
            if len(xs) < 2 or not xs[0]:
                continue
            try:
                source_taxid, c = do_one_line(scrambled_name_to_taxid, ncbi, xs)
            except KeyError as e:
                raise ValueError(input_file, xs, l)
            counts[c]+=1
            cs_for_species[source_taxid] = c

    add_stats(counts)
    return counts, cs_for_species

def do(refdb_ncbi, refdb_marker_to_taxon_path, input_files):
    scrambled_name_to_taxid = read_scrambled_name_to_taxid_from_marker_to_taxon(refdb_marker_to_taxon_path)

    ncbi = NCBITaxa(refdb_ncbi)
    lines = []
    lines.append("# Summary: ")
    lines.append("# | Name | " + " | ".join(header) + " |")
    lines.append("# | -- | " + " | ".join(["--" for h in header]) + " |")

    all_results = {}
    input_names = []
    for input_file in input_files:
        if ":" in input_file:
            (input_name, input_path) = input_file.split(":")
        else:
            input_name = input_file
            input_path = input_file
        input_names.append(input_name)
        counts, cs_for_species = do_one(scrambled_name_to_taxid, ncbi, input_path)
        lines.append("# | " + input_name + " | " + " | ".join([str(counts[h]) for h in header]) + " |")
        all_results[input_name] = cs_for_species

    all_taxids = set([x for xx in all_results.values() for x in xx ])
    taxid_to_name = ncbi.get_taxid_translator(all_taxids)

    lines.append("# Values legend: ")
    for h in header:
        if h not in header_shortcuts:
            continue
        lines.append("# " + header_shortcuts[h] + ": " + h)

    lines.append("#")
    lines.append("species\t" + "\t".join(input_names))
    for taxid in sorted(all_taxids):
        lines.append(str(taxid) + "|" + taxid_to_name[taxid] + "\t" + "\t".join([sc(all_results, input_name, taxid) for input_name in input_names]))
    for line in lines:
        print(line)



def main(argv=sys.argv[1:]):
    parser = argparse.ArgumentParser(
      description="simulate_and_align",
      formatter_class = argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--refdb-marker-to-taxon-path", type=str, action="store", dest="refdb_marker_to_taxon_path", help = "Lookup file, two columns - marker name, taxon name", required = True)
    parser.add_argument("--refdb-ncbi", type=str, action="store", dest="refdb_ncbi", help = "argument for ete.NCBITaxa", required = True)
    parser.add_argument("--input", type=str, action="append", dest="input_files", help = "results summary inputs", required = True)



    options=parser.parse_args(argv)
    do(**vars(options))


if __name__ == '__main__':
    main()
