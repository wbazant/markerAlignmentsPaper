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

header_data = [
  ('No results', 'NR'),
  ('has A, has B, has others', 'ABO'),
  ('has A, has B, no others', 'ABo'),
  ('has A, misses B, has others', 'AbO'),
  ('has A, misses B, no others', 'Abo'),
  ('misses A, has B, has others', 'aBO'),
  ('misses A, has B, no others', 'aBo'),
  ('misses A, has B, no others', 'aBO'),
]

header_shortcuts = dict(header_data)
header = [x for x,y in header_data]

def frac_to_score(x):
    return int(round(x*5, 0)) or 1

def scorecard(df_summary):
    result = pandas.DataFrame()
    result["Method"] = df_summary["Name"]
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

    return counts_with_stats

def sc(all_results, input_file, taxid_a, taxid_b):
    if (taxid_a, taxid_b) not in all_results[input_file]:
        return "XXX"
    return header_shortcuts[all_results[input_file][(taxid_a, taxid_b)]]

def do_one_line(ncbi, xs):
    sample = xs.pop(0)

    m = re.search('taxonA(.*)taxonB(.*)', sample)
    if not m:
        raiseValueError(sample)
    taxon_a = int(m.group(1))
    taxon_b = int(m.group(2))

    n = int(xs.pop(0))
    if n:
        vs = [int(ncbi.get_taxid_from_string(x)) for x in xs]
        result = []
        result.append("has A" if taxon_a in vs else "misses A")
        result.append("has B" if taxon_b in vs else "misses B")
        result.append("has others" if set([taxon_a, taxon_b]) != set(vs) else "no others")
        return taxon_a, taxon_b, ", ".join(result)
    else:
        return taxon_a, taxon_b, "No results"

def do_one(ncbi, input_file):
    counts = {h:0 for h in header}
    cs_for_species = {}
    with open(input_file, 'r') as f:
        for l in f:
            xs = l.rstrip().split("\t")
            if len(xs) < 2 or not xs[0]:
                continue
            try:
                taxon_a, taxon_b, c = do_one_line(ncbi, xs)
            except KeyError as e:
                raise ValueError(input_file, xs, l)
            counts[c]+=1
            cs_for_species[(taxon_a, taxon_b)] = c

    counts_with_stats = add_stats(counts)
    return counts_with_stats, cs_for_species

def read_all(refdb_ncbi, refdb_marker_to_taxon_path,  input_files, **kwargs):
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
        counts, cs_for_species = do_one(ncbi, input_path)
        all_counts[input_name] = counts
        all_results[input_name] = cs_for_species

    all_taxids = set([x for xxx in all_results.values() for xx in xxx for x in xx ])
    all_pairs = set([x for xx in all_results.values() for x in xx ])
    taxid_to_name = ncbi.get_taxid_translator(all_taxids)

    columns_detailed = ["taxon_a", "taxon_b"] + input_names
    data_detailed = [[str(taxid_a) + "|" + taxid_to_name[taxid_a], str(taxid_b) + "|" + taxid_to_name[taxid_b] ] + [sc(all_results, input_name, taxid_a, taxid_b) for input_name in input_names] for taxid_a, taxid_b in sorted(all_pairs)]
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
