import pandas
import argparse
import sys
import matplotlib.pyplot as plt

import re
import os
sys.path.append(os.path.join(os.path.realpath(os.path.dirname(__file__)), "../"))
from lib.ncbi2 import NCBITaxa2
from lib.utils import read_scrambled_name_to_taxid_from_marker_to_taxon

header_data = [
  ('has A, has B, has others, good LCA, good genus', 'ABOLG'),
  ('has A, has B, has others, no good LCA, good genus', 'ABOlG'),
  ('has A, has B, has others, good LCA, no good genus', 'ABOLg'),
  ('has A, has B, has others, no good LCA, no good genus', 'ABOlg'),
  ('has A, has B, has others, good LCA, n/a genus', 'ABOL'),
  ('has A, has B, has others, no good LCA, n/a genus', 'ABOl'),
  ('has A, has B, no others', 'ABo'),
  ('has A, misses B, has others', 'AbO'),
  ('has A, misses B, no others', 'Abo'),
  ('misses A, has B, has others', 'aBO'),
  ('misses A, has B, no others', 'aBo'),
  ('misses A, misses B, has others', 'abO'),
  ('No results', 'NR'),
]

header_shortcuts = dict(header_data)
header = [x for x,y in header_data]

def add_stats(counts):
    s = sum(counts.values())

    counts_with_stats = counts.copy()

    return counts_with_stats

def sc(ncbi, all_results, input_file, taxid_a, taxid_b, lca_taxid, lca_name, lca_rank, lca_genus):
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

    pair_lca_taxid, pair_lca_name, pair_lca_rank, pair_lca_genus = ncbi.lca_info([taxon_a, taxon_b])
    n = int(xs.pop(0))
    if n:
        vs = [int(ncbi.get_taxid_from_string(x)) for x in xs]
        result = []
        result.append("has A" if taxon_a in vs else "misses A")
        result.append("has B" if taxon_b in vs else "misses B")
        if all([v in set([taxon_a, taxon_b]) for v in vs]):
            result.append("no others")
        else:
            result.append("has others")
            if taxon_a in vs and taxon_b in vs:
                result_lca_taxid, result_lca_name, result_lca_genus, result_lca_genus = ncbi.lca_info(sorted(vs))
                result.append("good LCA" if result_lca_taxid == pair_lca_taxid else "no good LCA")
                result.append("good genus" if pair_lca_genus and pair_lca_genus == result_lca_genus else "no good genus" if pair_lca_genus else "n/a genus")
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

def read_all(refdb_ncbi, input_files, **kwargs):
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

    columns_detailed = ["taxon_a", "taxon_b", "taxon_lca_ab", "rank_lca_ab", "genus_lca_ab"] + input_names
    data_detailed = []
    for taxid_a, taxid_b in sorted(all_pairs):
        lca_taxid, lca_name, lca_rank, lca_genus = ncbi.lca_info([taxid_a, taxid_b])
        ids = [str(taxid_a) + "|" + taxid_to_name[taxid_a], str(taxid_b) + "|" + taxid_to_name[taxid_b], str(lca_taxid) + "|" +lca_name, lca_rank, lca_genus ]
        results = [sc(ncbi, all_results, input_name, taxid_a, taxid_b, lca_taxid, lca_name, lca_rank, lca_genus) for input_name in input_names]
        data_detailed.append(ids+results)
    df_detailed = pandas.DataFrame(columns = columns_detailed, data = data_detailed)

    columns_summary = ["Name"] + ["{}: {}".format(header_shortcuts[h], h) for h in header]
    data_summary = [[input_name] + [ all_counts[input_name][h] for h in header ] for input_name in input_names]
    df_summary = pandas.DataFrame(columns = columns_summary, data = data_summary)
    return df_detailed, df_summary
#ABOLG: has A, has B, has others, good LCA, good genus
#ABOlG: has A, has B, has others, no good LCA, good genus
#ABOLg: has A, has B, has others, good LCA, no good genus
#ABOlg: has A, has B, has others, no good LCA, no good genus
#ABOL: has A, has B, has others, good LCA, n/a genus
#ABOl: has A, has B, has others, no good LCA, n/a genus
#ABo: has A, has B, no others
#AbO: has A, misses B, has others
#Abo: has A, misses B, no others
#aBO: misses A, has B, has others
#aBo: misses A, has B, no others
#abO: misses A, misses B, has others
#NR: No results
colors = {
  "ABOLG": "#40E0D0", # turquoise
  "ABOlG": "#A7C7E7", # pastel blue
  "ABOLg": "#40E0D0", # turquoise
  "ABOlg": "#954535", # chestnut
  "ABOL": "#40E0D0",  # turquoise
  "ABOl": "#F88379",  # coral pink
  "ABo": "blue",
  "AbO": "#FAFA33", # lemon yellow
  "Abo": "#FFE5B4", # peach
  "aBO": "#FAFA33", # lemon yellow
  "aBo": "#C9CC3F", # pear
  "abO": "#FAFA33", # lemon yellow
  "NR": "grey",
}
def cat_to_color(x):
    return colors[x] if x in colors else "grey" # 4 XXX datapoints at 0.01 where we failed to draw any reads


def subplot(df, column, label, ax):
    df.plot.scatter("pca_1", "pca_2", c = df[column].map(cat_to_color), label = label, ax = ax)
    ax.axis('off')
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.3, box.width, box.height*0.7])
    ax.legend(loc='upper center', edgecolor = "white", markerscale = 0, bbox_to_anchor=(0.5, -0.05))

#https://stackoverflow.com/a/40535454
def fix_column_width(writer, sheet_name, df):
    worksheet = writer.sheets[sheet_name]  # pull worksheet object
    lengths = [len(str(df[col].name)) for idx, col in enumerate(df)]
    worksheet.set_column(0, len(lengths), max(lengths))

def do(refdb_ncbi, input_files, pca_tsv, output_tsv, output_png, output_xlsx):
    results_df_detailed, results_df_summary = read_all(refdb_ncbi, input_files)
    results_df_detailed.to_csv(output_tsv, sep = "\t", index = False)

    results_df_detailed = results_df_detailed.drop(["taxon_lca_ab", "rank_lca_ab", "genus_lca_ab"], axis=1)
    results_df_detailed['taxon_a'] =  results_df_detailed['taxon_a'].map(lambda x: x.split("|")[0])
    results_df_detailed['taxon_b'] =  results_df_detailed['taxon_b'].map(lambda x: x.split("|")[0])

    pca_df = pandas.read_csv(pca_tsv, sep = "\t").astype({"taxon_a": str, "taxon_b": str})
    df = pca_df.merge(results_df_detailed, on = ["taxon_a", "taxon_b"])

    fig, axs = plt.subplots(3, 2)
    subplot(df, "EukDetectLo", "EukDetect, 0.01 cov.", axs[0][0])
    subplot(df, "EukDetectMid", "EukDetect, 0.05 cov.", axs[1][0])
    subplot(df, "EukDetectHi", "EukDetect, 0.10 cov.", axs[2][0])
    subplot(df, "CORRALLo", "CORRAL, 0.01 cov.", axs[0][1])
    subplot(df, "CORRALMid", "CORRAL, 0.05 cov.", axs[1][1])
    subplot(df, "CORRALHi", "CORRAL, 0.10 cov.", axs[2][1])

    fig.savefig(output_png, bbox_inches='tight', dpi=199)

    df = pca_df.merge(results_df_detailed, on = ["taxon_a", "taxon_b"], how='left')
    df.fillna("")
    cs = df.columns
    df = df.set_index(['taxon_a_name', 'taxon_b_name'])
    
    with pandas.ExcelWriter(output_xlsx, engine="xlsxwriter") as writer:
        sheet_name_summary = "Coding legend and result counts"
        results_df_summary.transpose().to_excel(writer, sheet_name=sheet_name_summary, header = False)
        worksheet_summary = writer.sheets[sheet_name_summary]
        cell_format_summary = writer.book.add_format({'bold': True})
        worksheet_summary.set_column(0, 0, 30, cell_format_summary )
        worksheet_summary.set_column(1, 10, 20)
        worksheet_summary.set_row(0, 20, cell_format_summary )

        sheet_name_detailed = "PCA coordinates and results"
        df.to_excel(writer, sheet_name = sheet_name_detailed,  float_format="%.3f")
        worksheet_detailed = writer.sheets[sheet_name_detailed]
        cell_format_detailed = writer.book.add_format({'bold': True})
        worksheet_detailed.set_column(0, 1, 30,cell_format_detailed )
        worksheet_detailed.set_column(1, 100, 20)


def opts(argv):
    parser = argparse.ArgumentParser(
      description="parse pairs results and plot / spreadsheet them",
      formatter_class = argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--refdb-ncbi", type=str, action="store", dest="refdb_ncbi", help = "argument for ete.NCBITaxa", required = True)
    parser.add_argument("--input", type=str, action="append", dest="input_files", help = "results summary inputs", required = True)
    parser.add_argument("--pca-tsv", type=str, action="store", dest="pca_tsv", required=True)
    parser.add_argument("--output-tsv", type=str, action="store", dest="output_tsv", help = "result tsv", required = True)
    parser.add_argument("--output-png", type=str, action="store", dest="output_png", required=True)
    parser.add_argument("--output-xlsx", type=str, action="store", dest="output_xlsx", required=True)
    return parser.parse_args(argv)

def main(argv=sys.argv[1:]):
    options = opts(argv)
    kwargs = vars(options)
    do(**kwargs)

if __name__ == '__main__':
    main()
