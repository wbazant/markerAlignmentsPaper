import sys
import os
import argparse
import pysam
import re
import pandas


import os
sys.path.append(os.path.join(os.path.realpath(os.path.dirname(__file__)), "../"))
from lib.ncbi2 import NCBITaxa2
from lib.utils import read_scrambled_name_to_taxid_from_marker_to_taxon

def ct(n, rank, all_good, some_good):
    result = []
    result.append("One result" if n == 1 else "Many results")
    result.append("same genus" if rank == "genus" else "not same genus")
    result.append("all taxonomically closest" if all_good else "part taxonomically closest" if some_good else "not taxonomically closest")
    return ", ".join(result)

header_data = [
  ('No results', 'NR'),
  ('Multiple taxa including correct', 'MC'),
  ('Single correct taxon', 'SC'),
  ('Multiple incorrect taxa', 'MI'),
  ('Single incorrect taxon', 'SI'),
  ('Outcomes containing correct taxon', '(SC+MC)/total'),
  ('Ratio of single correct taxon to all reported results', 'SC/(total-NR)'),
]

header_shortcuts = dict(header_data)
header = [x for x,y in header_data]

def frac_to_score(x):
    return int(round(x*5, 0)) or 1

def scorecard(df_summary):
    result = pandas.DataFrame()
    result["Method"] = df_summary["Name"]
    result["Detects signal"] = df_summary["(SC+MC)/total: Outcomes containing correct taxon"].map(frac_to_score)
    result["Reports only correct result"] = df_summary["SC/(total-NR): Ratio of single correct taxon to all reported results"].map(frac_to_score)
    return result
        

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

    counts_with_stats = counts.copy()
    counts_with_stats["Outcomes containing correct taxon"] = round(1.0 * (counts["Single correct taxon"] + counts["Multiple taxa including correct"] ) / s, 3)
    counts_with_stats["Ratio of single correct taxon to all reported results"] = round(1.0 * counts["Single correct taxon"] / (s - counts["No results"] ), 3) if s > counts["No results"] else 0

    return counts_with_stats

def sc(all_results, input_file, taxid):
    if taxid not in all_results[input_file]:
        return "XXX"
    return header_shortcuts[all_results[input_file][taxid]]


def do_one_line(scrambled_name_to_taxid, ncbi, xs):
    source_taxon = xs.pop(0)
    source_taxid = scrambled_name_to_taxid[source_taxon]
    n = int(xs.pop(0))
    if n:
        vs = [int(ncbi.get_taxid_from_string(x)) for x in xs]
        has_correct = source_taxid in vs
        has_many = len(vs) > 1
        if has_correct and has_many:
            return source_taxid, "Multiple taxa including correct"
        elif has_correct and not has_many:
            return source_taxid, "Single correct taxon"
        elif not has_correct and has_many:
            return source_taxid, "Multiple incorrect taxa"
        elif not has_correct and not has_many:
            return source_taxid, "Single incorrect taxon"
        else:
            raise ValueError(source_taxid)
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

    counts_with_stats = add_stats(counts)
    return counts_with_stats, cs_for_species

def read_all(refdb_ncbi, refdb_marker_to_taxon_path,  input_files, **kwargs):
    scrambled_name_to_taxid = read_scrambled_name_to_taxid_from_marker_to_taxon(refdb_marker_to_taxon_path)
    ncbi = NCBITaxa2(refdb_ncbi)

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
        counts, cs_for_species = do_one(scrambled_name_to_taxid, ncbi, input_path)
        all_counts[input_name] = counts
        all_results[input_name] = cs_for_species

    all_taxids = set([x for xx in all_results.values() for x in xx ])
    taxid_to_name = ncbi.get_taxid_translator(all_taxids)

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

def do(refdb_ncbi, refdb_marker_to_taxon_path, input_files, output_tsv, output_xlsx):
    df_detailed, df_summary = read_all(refdb_ncbi, refdb_marker_to_taxon_path, input_files)

    df_detailed.to_csv(output_tsv, sep = "\t", index = False)

    if output_xlsx:
        df_scorecard = scorecard(df_summary)
        with pandas.ExcelWriter(output_xlsx, engine="xlsxwriter") as writer:
            sheet_name_scorecard = "Scorecard"
            df_scorecard.to_excel(writer, sheet_name=sheet_name_scorecard, index = False)
            worksheet_scorecard = writer.sheets[sheet_name_scorecard]
            cell_format_scorecard = writer.book.add_format({'bold': True})
            worksheet_scorecard.set_column(0, 0, 30,cell_format_scorecard )
            fix_column_width(writer, sheet_name=sheet_name_scorecard, df = df_scorecard)

            sheet_name_summary = "Counts and stats"
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


def parse(argv):
    parser = argparse.ArgumentParser(
      description="parse_whole_samples_results",
      formatter_class = argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--refdb-marker-to-taxon-path", type=str, action="store", dest="refdb_marker_to_taxon_path", help = "Lookup file, two columns - marker name, taxon name", required = True)
    parser.add_argument("--refdb-ncbi", type=str, action="store", dest="refdb_ncbi", help = "argument for ete.NCBITaxa", required = True)
    parser.add_argument("--input", type=str, action="append", dest="input_files", help = "results summary inputs", required = True)
    parser.add_argument("--output-tsv", type=str, action="store", dest="output_tsv", help = "result tsv", required = True)
    parser.add_argument("--output-xlsx", type=str, action="store", dest="output_xlsx", help = "result xlsx", default = None)
    return parser.parse_args(argv)


def main(argv=sys.argv[1:]):
    options=parse(argv)
    do(**vars(options))


if __name__ == '__main__':
    main()
