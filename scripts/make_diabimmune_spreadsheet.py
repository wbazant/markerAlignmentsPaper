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



def do(triples_our_markers, triples_our_reads, csv_eukdetect, output_xlsx):
    with pandas.ExcelWriter(output_xlsx, engine="xlsxwriter") as writer:
        pandas.DataFrame(data=[1]).to_excel(writer, sheet_name = "TODO")
    
            

def opts(argv):
    parser = argparse.ArgumentParser(
      description="make diabimmune comparison spreadsheet",
      formatter_class = argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input-tsv-triples-our-num-markers", type=str, action="store", dest="triples_our_markers", required=True)
    parser.add_argument("--input-tsv-triples-our-num-reads", type=str, action="store", dest="triples_our_reads", required=True)
    parser.add_argument("--input-csv-eukdetect", type=str, action="store", dest="csv_eukdetect", required=True)

    parser.add_argument("--output-xlsx", type=str, action="store", dest="output_xlsx", required=True)

    return parser.parse_args(argv)

def main(argv=sys.argv[1:]):
    options = opts(argv)
    do(**vars(options))

if __name__ == '__main__':
    main()
