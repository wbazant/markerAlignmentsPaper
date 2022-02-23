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

def do(input_tsvs_with_per_species_data, input_jsons_with_stats, stats_tab_name, output_xlsx):
    with pandas.ExcelWriter(output_xlsx, engine="xlsxwriter") as writer:
        for input_file in input_tsvs_with_per_species_data:
            print(input_file)
            if ":" in input_file:
                (input_name, input_path) = input_file.split(":")
            else:
                input_name = input_file
                input_path = input_file
            df = pandas.read_csv(input_path, sep="\t")
            df.to_excel(writer, sheet_name=input_name, index = False)
            fix_column_width(writer, sheet_name=input_name, df = df)

        stats_dfs = []
        for input_file in input_jsons_with_stats:
            print(input_file)
            if ":" in input_file:
                (input_name, input_path) = input_file.split(":")
            else:
                input_name = input_file
                input_path = input_file
            df = pandas.io.json.read_json(input_path)
            ix_cs = ["mutation_rate", "base_error_rate", "read_length"]
            cs = ix_cs + [c for c in df if c not in ix_cs]
            cs = [c for c in cs if not c.startswith("num") and c not in ["base_error_rate", "read_length", "query_lengthsAvg", "query_lengthsFractionAtLeast60"] and not c.endswith("Buscos")]
            cs_rename = {c:c.replace("Queries", "") for c in cs if c.endswith("Queries")}
            df = df[cs]
            df = df.rename(cs_rename, axis = "columns")
            df.insert(loc=0, column='Group', value=input_name)
            stats_dfs.append(df)
        stats_df = pandas.concat(stats_dfs)
        stats_df.to_excel(writer, sheet_name = stats_tab_name, index = False, float_format="%.3f", na_rep=0)
        worksheet = writer.sheets[stats_tab_name]
        cell_format = writer.book.add_format({'bold': True})
        worksheet.set_column(0, 1, 30,cell_format )
        fix_column_width(writer, sheet_name=stats_tab_name, df = stats_df)


            

def opts(argv):
    parser = argparse.ArgumentParser(
      description="make simulated whole reads spreadsheet",
      formatter_class = argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input-tsv-for-per-species-tab", type=str, action="append", dest="input_tsvs_with_per_species_data", help = "per species tab", required = True)
    parser.add_argument("--input-json-for-stats-tab", type=str, action="append", dest="input_jsons_with_stats", help = "per species tab", required = True)
    parser.add_argument("--stats-tab-name", type=str, action="store", dest="stats_tab_name", required=True)
    parser.add_argument("--output-xlsx", type=str, action="store", dest="output_xlsx", required=True)

    return parser.parse_args(argv)

def main(argv=sys.argv[1:]):
    options = opts(argv)
    do(**vars(options))

if __name__ == '__main__':
    main()
