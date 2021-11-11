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

import matplotlib.pyplot as plt
import pandas

from marker_alignments.store import SqliteStore
from ete3 import NCBITaxa
#pandas.read_sql_query(sql)

sql = '''
         select source_taxon, avg(mapq) as avg_mapq, 1.0 * sum(mapq_at_least_30) / count(*) as fraction_mapq_at_least_30, 1.0 * sum(is_match) / count(*) as precision from (
              select a.query, a.source_taxon, a.matched_taxon, mapq, case when mapq >=30 then 1 else 0 end as mapq_at_least_30, case when sum(c) == count(c) then 1 else 0 end as is_match
              from (
                select query, source_taxon, matched_taxon, mapq, case when match_type in ('true_match', 'species') then 1 else 0 end as c
                from alignment_from_known_source
              ) a
              group by a.query, source_taxon, matched_taxon, mapq
        ) group by source_taxon
        order by avg_mapq
'''
def do(input_db, output_png, refdb_ncbi):

    ncbi = NCBITaxa(refdb_ncbi)
    df = pandas.read_sql_query(sql, sqlite3.connect(input_db))
    ax = df.plot.scatter(x='fraction_mapq_at_least_30', y='precision')
    fig = ax.get_figure()
    fig.savefig(output_png, bbox_inches='tight', dpi=199)
    


def main(argv=sys.argv[1:]):
    parser = argparse.ArgumentParser(
      description="plot precision by species",
      formatter_class = argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input-alignments-sqlite", type=str, action="store", dest="input_db", required=True)
    parser.add_argument("--output-png", type=str, action="store", dest="output_png", required=True)
    parser.add_argument("--refdb-ncbi", type=str, action="store", dest="refdb_ncbi", help = "argument for ete.NCBITaxa")

    options=parser.parse_args(argv)


    do(**vars(options))
    


if __name__ == '__main__':
    main()
