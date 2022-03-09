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
import itertools
import numpy as np

import matplotlib.pyplot as plt
import xlsxwriter
import pandas

#https://stackoverflow.com/a/40535454
def fix_column_width(writer, sheet_name, df):
    worksheet = writer.sheets[sheet_name]  # pull worksheet object
    lengths = [len(str(df[col].name)) for idx, col in enumerate(df)]
    worksheet.set_column(0, len(lengths), max(lengths))

def read_sample_to_run(input_file):
    result = {}
    with open(input_file, 'r') as f:
        for l in f.readlines():
            (sample, run) = l.rstrip().split("\t")
            result[sample] = run
    return result

#3112358	766728|Meyerozyma	3
def read_triples(ncbi, key, value_type, input_file, sr):
    lines = []
    with open(input_file, 'r') as f:
        for l in f.readlines():
            (sample, taxon_string, value_string) = l.rstrip().split("\t")
            value = value_type(value_string)
            taxon = ncbi.get_taxid_from_string(taxon_string) 
            lines.append({"sample": sample, "run": sr[sample], "taxon" : taxon, "our_result": taxon_string, key: value})
    return pandas.DataFrame.from_records(lines)

            
def read_eukdetect(csv_eukdetect):
    c = pandas.read_csv(csv_eukdetect)
    return pandas.DataFrame.from_dict({
        "run": c["Sample ID"],
        "taxon": [int(x) for x in c["TaxID"]],
        "eukdetect_result": c["Name"],
        "eukdetect_markers": c["Observed_markers"],
        "eukdetect_reads": c["Read_counts"]})

def read(refdb_ncbi, triples_our_markers, triples_our_reads, csv_eukdetect, sample_to_run):
    ncbi = NCBITaxa2(refdb_ncbi)
    sr = read_sample_to_run(sample_to_run)
    df_our_markers = read_triples(ncbi, 'our_markers', int, triples_our_markers, sr)
    df_our_reads = read_triples(ncbi, 'our_reads', float, triples_our_reads, sr)
    df_ours =  pandas.merge(df_our_markers, df_our_reads, how = 'inner', on = ["sample", "run", "taxon", "our_result"])
    df_eukdetect = read_eukdetect(csv_eukdetect)
    df_ours.taxon = df_ours.taxon.astype(str)
    df_eukdetect.taxon = df_eukdetect.taxon.astype(str)

    m = pandas.merge(df_ours, df_eukdetect, how = 'outer', on = ["run", "taxon"])

    all_taxids = set(m["taxon"])
    taxid_to_name = ncbi.get_taxid_translator(all_taxids)
    m["taxon_name"] = [taxid_to_name[int(taxid)] for taxid in m["taxon"]]
    cols = m.columns.tolist()
    cols = cols[0:3] + cols[-1:] + cols[3:-1]
    m = m[cols]
    m = m.set_index(["sample", "run", "taxon"])
    return m

#    pandas.merge(df_ours, df_eukdetect, how = 'outer', on = ["run", "taxon"])
    return df_ours, df_eukdetect

def do(refdb_ncbi, triples_our_markers, triples_our_reads, csv_eukdetect, sample_to_run, output_xlsx):
    m = read(refdb_ncbi, triples_our_markers, triples_our_reads, csv_eukdetect, sample_to_run)
    with pandas.ExcelWriter(output_xlsx, engine="xlsxwriter") as writer:
        sheet_name = "DIABIMMUNE results"
        m.to_excel(writer, sheet_name = sheet_name,  float_format="%.2f")
        worksheet = writer.sheets[sheet_name]
        fix_column_width(writer, sheet_name=sheet_name, df = m)
    
            

def opts(argv):
    parser = argparse.ArgumentParser(
      description="make diabimmune comparison spreadsheet",
      formatter_class = argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--refdb-ncbi", type=str, action="store", dest="refdb_ncbi", help = "argument for ete.NCBITaxa", required = True)
    parser.add_argument("--input-tsv-triples-our-num-markers", type=str, action="store", dest="triples_our_markers", required=True)
    parser.add_argument("--input-tsv-triples-our-num-reads", type=str, action="store", dest="triples_our_reads", required=True)
    parser.add_argument("--input-tsv-sample-to-run", type=str, action="store", dest="sample_to_run", required=True)
    parser.add_argument("--input-csv-eukdetect", type=str, action="store", dest="csv_eukdetect", required=True)

    parser.add_argument("--output-xlsx", type=str, action="store", dest="output_xlsx", required=True)

    return parser.parse_args(argv)

def main(argv=sys.argv[1:]):
    options = opts(argv)
    do(**vars(options))

if __name__ == '__main__':
    main()
