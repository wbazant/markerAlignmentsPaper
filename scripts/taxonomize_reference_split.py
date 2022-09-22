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

from ete3 import NCBITaxa
import pandas

from pathlib import Path


def main_loop(**kwargs):
    result = []
   
    ncbi = NCBITaxa(kwargs['refdb_ncbi'])
    nine_tenth = pandas.read_csv(kwargs['nine_tenth_path'], names = ["Species", "Taxid"])
    one_tenth = pandas.read_csv(kwargs['one_tenth_path'], names = ["Species", "Taxid"])
    df = get_df(ncbi, one_tenth, nine_tenth)
    df.to_csv(kwargs['output_path'], sep = "\t", index = False)

def get_df(ncbi, one_tenth, nine_tenth): 
    nine_tenth_labels = dict(zip(nine_tenth["Taxid"], nine_tenth["Species"]))
    result = []
    for species, taxid in zip(one_tenth["Species"], one_tenth["Taxid"]):
        t = ncbi.get_topology(set(nine_tenth['Taxid']) | {taxid})
        L = t & taxid
        descendants = L.up.get_descendants()

        nearest_common_rank = L.up.rank if L.up.rank != 'no rank' else ncbi.get_taxid_translator([L.up.taxid])[L.up.taxid] + " (no rank)"

        for x in descendants:
            if x.is_leaf() and x.taxid != taxid:
                result.append([str(taxid), species, nearest_common_rank, str(x.taxid), nine_tenth_labels[x.taxid]])
            else:
                pass
    return pandas.DataFrame.from_records(result, columns = ["Holdout taxid","Holdout name", "rank", "Reference taxid", "Reference name"]) 

def main(argv=sys.argv[1:]):
    parser = argparse.ArgumentParser(
      description="taxonomize_reference_split",
      formatter_class = argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--refdb-ncbi", type=str, action="store", dest="refdb_ncbi", help = "argument for ete.NCBITaxa")
    parser.add_argument("--holdout-species-path", type=str, action="store", dest="one_tenth_path", help = "Species being held out")
    parser.add_argument("--remaining-species-path", type=str, action="store", dest="nine_tenth_path", help = "Species remaining in reference")
    parser.add_argument("--out-tsv-path", type=str, action="store", dest="output_path", help = "Output")

    options=parser.parse_args(argv)


    logging.basicConfig(format='%(asctime)s: %(message)s')
    logger=logging.getLogger(__name__)

    main_loop(**vars(options), logger = logger )
    


if __name__ == '__main__':
    main()
