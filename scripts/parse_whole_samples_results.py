import sys
import re
import os
sys.path.append(os.path.join(os.path.realpath(os.path.dirname(__file__)), "../"))
from lib.ncbi2 import NCBITaxa2
import argparse
import pysam
import re
import json
import statistics
import logging
import subprocess
import sqlite3
import pandas

from marker_alignments.store import SqliteStore

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


def ct(n, rank, all_good):
    result = []
    result.append("One result" if n == 1 else "Many results")
    result.append("same genus" if rank == "genus" else "not same genus")
    result.append("taxonomically closest" if all_good else "not taxonomically closest")
    return ", ".join(result)

header_data = [   ('No results', 'NR'),
    ('One result, same genus, all taxonomically closest', 'PGC'),
    ('Many results, same genus, all taxonomically closest', 'pGC'),
    ('One result, not same genus, all taxonomically closest', 'PgC'),
    ('Many results, not same genus, all taxonomically closest', 'pgC'),
    ('One result, same genus, not taxonomically closest', 'PGc'),
    ('Many results, same genus, not taxonomically closest', 'pGc'),
    ('One result, not same genus, not taxonomically closest', 'Pgc'),
    ('Many results, not same genus, not taxonomically closest', 'pgc'),
    ('Signal detected', 'SD = (total - NR) / total'),
    ('Detected signal is one species', 'P.. / SD'),
    ('Detected signal is same genus', '.G. / SD'),
    ('Detected signal is one same genus species', 'PG. / SD'),
    ('Detected signal is taxonomically closest', '..C / SD'),
    ('Detected signal is one taxonomically closest species', 'P.C / SD'),
    (   'Detected signal is one same genus taxonomically closest species',
        'PGC / SD')]

header_shortcuts = dict(header_data)
header = [x for x,y in header_data]

def pat(xs, counts):
    s = 0
    for l, c in counts.items():
        if l == 'No results' or c == 0 or len(l.split(", ")) < len(xs):
            continue
        m = True
        for n in range(0, len(xs)):
            m = m and (xs[n] == "." or xs[n] == l.split(", ")[n])
        if m:
            s += c
    return s


def add_stats(counts):
    s = sum(counts.values())
    S = pat([".", ".", "."], counts)

    counts_with_stats = counts.copy()
    counts_with_stats["Signal detected"] = round(1.0 * pat([".", ".", "."], counts) / s, 3)
    counts_with_stats["Detected signal is one species"] = round(1.0 * pat(["One result", ".", "."], counts) / pat([".", ".", "."], counts), 3)
    counts_with_stats["Detected signal is same genus"] = round(1.0 * pat([".", "same genus", "."], counts) / pat([".", ".", "."], counts), 3)
    counts_with_stats["Detected signal is one same genus species"] = round(1.0 * pat(["One result", "same genus", "."], counts) / pat([".", ".", "."], counts), 3)
    counts_with_stats["Detected signal is all taxonomically closest"] = round(1.0 * pat([".", ".", "all taxonomically closest"], counts) / pat([".", ".", "."], counts), 3)
    counts_with_stats["Detected signal is one taxonomically closest species"] = round(1.0 * pat(["One result", ".", "all taxonomically closest"], counts) / pat([".", ".", "."], counts), 3)
    counts_with_stats["Detected signal is one same genus taxonomically closest species"] = round(1.0 * pat(["One result", "same genus", "all taxonomically closest"], counts) / pat([".", ".", "."], counts), 3)

    return counts_with_stats

def sc(all_results, input_file, taxid):
    if taxid not in all_results[input_file]:
        return "XXX"
    return header_shortcuts[all_results[input_file][taxid]]


def do_one_line(scrambled_name_to_taxid, ncbi, good_matches, xs):
    source_taxon = xs.pop(0)
    source_taxid = scrambled_name_to_taxid[source_taxon]
    n = int(xs.pop(0))
    if n:
        vs = [ncbi.get_taxid_from_string(x) for x in xs]
        all_good = all([(source_taxid, int(v)) in good_matches for v in vs]) 
        vs.append(str(source_taxid))
        rank = ncbi.get_match_type(vs)
        return source_taxid, ct(n, rank, all_good)
    else:
        return source_taxid, "No results"

def do_one(scrambled_name_to_taxid, ncbi, good_matches, input_file):
    counts = {h:0 for h in header}
    cs_for_species = {}
    with open(input_file, 'r') as f:
        for l in f:
            xs = l.rstrip().split("\t")
            if len(xs) < 2 or not xs[0]:
                continue
            try:
                source_taxid, c = do_one_line(scrambled_name_to_taxid, ncbi, good_matches, xs)
            except KeyError as e:
                raise ValueError(input_file, xs, l)
            counts[c]+=1
            cs_for_species[source_taxid] = c

    counts_with_stats = add_stats(counts)
    return counts_with_stats, cs_for_species

def read_good_matches(good_matches_path):
    df = pandas.read_csv(good_matches_path, sep = "\t")
    return set(zip(df["Holdout taxid"].astype(int), df["Reference taxid"].astype(int)))

def read_all(refdb_ncbi, refdb_marker_to_taxon_path, good_matches_path, input_files):
    scrambled_name_to_taxid = read_scrambled_name_to_taxid_from_marker_to_taxon(refdb_marker_to_taxon_path)
    good_matches = read_good_matches(good_matches_path)

    ncbi = NCBITaxa2(refdb_ncbi)
    lines = []
    lines.append("# Summary: ")
    lines.append("# | Name | " + " | ".join(header) + " |")
    lines.append("# | -- | " + " | ".join(["--" for h in header]) + " |")

    all_results = {}
    input_names = []
    all_counts = {}
    for input_file in input_files:
        if ":" in input_file:
            (input_name, input_path) = input_file.split(":")
        else:
            input_name = input_file
            input_path = input_file
        input_names.append(input_name)
        counts, cs_for_species = do_one(scrambled_name_to_taxid, ncbi, good_matches, input_path)
        all_counts[input_name] = counts
        all_results[input_name] = cs_for_species

    all_taxids = set([x for xx in all_results.values() for x in xx ])
    taxid_to_name = ncbi.get_taxid_translator(all_taxids)

#    return (input_names, all_results, all_taxids, taxid_to_name)

    columns_detailed = ["species"] + input_names
    data_detailed = [[str(taxid) + "|" + taxid_to_name[taxid]] + [sc(all_results, input_name, taxid) for input_name in input_names] for taxid in sorted(all_taxids)]
    df_detailed = pandas.DataFrame(columns = columns_detailed, data = data_detailed)

    columns_summary = ["Name"] + ["{}: {}".format(header_shortcuts[h], h) for h in header]
    data_summary = [[input_name] + [ all_counts[input_name][h] for h in header ] for input_name in input_names]
    df_summary = pandas.DataFrame(columns = columns_summary, data = data_summary)
    return df_detailed, df_summary


#https://stackoverflow.com/a/40535454
def fix_column_width(writer, sheet_name, df):
    worksheet = writer.sheets[sheet_name]  # pull worksheet object
    lengths = [len(str(df[col].name)) for idx, col in enumerate(df)]
    worksheet.set_column(0, len(lengths), max(lengths))

def do(refdb_ncbi, refdb_marker_to_taxon_path, good_matches_path, input_files, output_tsv, output_xlsx):
    df_detailed, df_summary = read_all(refdb_ncbi, refdb_marker_to_taxon_path, good_matches_path, input_files)

    df_detailed.to_csv(output_tsv, sep = "\t", index = False)

    if output_xlsx:
        with pandas.ExcelWriter(output_xlsx, engine="xlsxwriter") as writer:
            sheet_name_summary = "Summary"
            df_summary.to_excel(writer, sheet_name=sheet_name_summary, index = False)
            worksheet_summary = writer.sheets[sheet_name_summary]
            cell_format_summary = writer.book.add_format({'bold': True})
            worksheet_summary.set_column(0, 0, 30,cell_format_summary )
            fix_column_width(writer, sheet_name=sheet_name_summary, df = df_summary)

            sheet_name_detailed = "Results by species"
            df_detailed.to_excel(writer, sheet_name=sheet_name_detailed, index = False)
            worksheet_detailed = writer.sheets[sheet_name_detailed]
            cell_format_detailed = writer.book.add_format({'bold': True})
            worksheet_detailed.set_column(0, 0, 30,cell_format_detailed )
            fix_column_width(writer, sheet_name=sheet_name_detailed, df = df_detailed)




def main(argv=sys.argv[1:]):
    parser = argparse.ArgumentParser(
      description="parse_whole_samples_results",
      formatter_class = argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--refdb-marker-to-taxon-path", type=str, action="store", dest="refdb_marker_to_taxon_path", help = "Lookup file, two columns - marker name, taxon name", required = True)
    parser.add_argument("--refdb-ncbi", type=str, action="store", dest="refdb_ncbi", help = "argument for ete.NCBITaxa", required = True)
    parser.add_argument("--good-matches-path", type=str, action="store", dest="good_matches_path", help = "List of matches we consider good", required = True)
    parser.add_argument("--input", type=str, action="append", dest="input_files", help = "results summary inputs", required = True)
    parser.add_argument("--output-tsv", type=str, action="store", dest="output_tsv", help = "result tsv", required = True)
    parser.add_argument("--output-xlsx", type=str, action="store", dest="output_xlsx", help = "result xlsx", default = None)



    options=parser.parse_args(argv)
    do(**vars(options))


if __name__ == '__main__':
    main()
