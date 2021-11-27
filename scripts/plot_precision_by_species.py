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
import pandas

from marker_alignments.store import SqliteStore
from ete3 import NCBITaxa
#pandas.read_sql_query(sql)

sql = '''
select 
  source_taxon,
  avg(mapq) as avg_mapq,
  1.0 * sum(mapq_at_least_30) / count(*) as fraction_mapq_at_least_30,
  1.0 * sum(is_match) / count(*) as precision,
  1.0 * sum(is_match_and_mapq_at_least_30) / sum(mapq_at_least_30) as precision_mapq_at_least_30
from (
  select
    a.query,
    a.source_taxon,
    a.matched_taxon,
    mapq,
    case when mapq >=30 then 1 else 0 end as mapq_at_least_30,
    case when sum(c) == count(c) then 1 else 0 end as is_match,
    case when sum(c) == count(c) and mapq >=30 then 1 else 0 end as is_match_and_mapq_at_least_30
  from (
    select query, source_taxon, matched_taxon, mapq, case when {MATCH_TYPE_CLAUSE} then 1 else 0 end as c
    from alignment_from_known_source
  ) a
    group by a.query, source_taxon, matched_taxon, mapq
) group by source_taxon
order by avg_mapq
'''

def match_types(same_what):
    refs = ['true_match', 'species', 'genus', 'family']
    if same_what not in refs:
        raise ValueError("Same what?: " + same_what)
    ix = refs.index(same_what)
    return "match_type in (" + ", ".join(["'{}'".format(x) for x in refs[:ix+1]]) + ")"

def get_data(input_db, same_what, refdb_ncbi, **kwargs):
    ncbi = NCBITaxa(refdb_ncbi)
    df = pandas.read_sql_query(sql.format(MATCH_TYPE_CLAUSE=match_types(same_what)), sqlite3.connect(input_db))
    df['kingdom'] = [get_kingdom(ncbi, x) for x in df['source_taxon']]
    df['delta_precision'] = df['precision_mapq_at_least_30'] - df['precision']

    print("All taxa that map: " + str(len(df)))
    print("Taxa that map perfectly: " + taxa_that_map_perfectly(ncbi, df))
    print("Taxa that map perfectly only with mapq filter: " + taxa_that_map_perfectly_only_with_mapq_filter(ncbi, df))
    print("Taxa that only map incorrectly: " + taxa_that_dont_map(ncbi, df))
    print("Taxa that only map incorrectly with mapq filter: " + taxa_that_dont_map_only_with_mapq_filter(ncbi, df))
    print("Taxa that get worse: " + taxa_that_get_worse(ncbi, df))
    print("Taxa that get wiped: " + taxa_that_get_wiped(ncbi, df))
    return df


def taxa_stat(ncbi, df, ix):
    l = list(itertools.chain.from_iterable([ncbi.get_taxid_translator([taxid]).values() for taxid in df[ix]['source_taxon']]))
    result = str(len(l))
    if len(l) < 10:
        result +=  ". " + ", ".join(l)
    return(str(len(l)))

def taxa_that_dont_map(ncbi, df):
    ix = df['precision'] < 0.0001
    return taxa_stat(ncbi, df, ix)

def taxa_that_dont_map_only_with_mapq_filter(ncbi, df):
    ix = (df['precision'] > 0.0001) & ( df['precision_mapq_at_least_30'] < 0.0001 )
    return taxa_stat(ncbi, df, ix)

def taxa_that_map_perfectly(ncbi, df):
    ix = df['precision'] >= 0.99999
    return taxa_stat(ncbi, df, ix)

def taxa_that_map_perfectly_only_with_mapq_filter(ncbi, df):
    ix = (df['precision'] < 0.99999 ) & (df['precision_mapq_at_least_30'] >= 0.99999 )
    return taxa_stat(ncbi, df, ix)

def taxa_that_get_worse(ncbi, df):
    ix = df['precision_mapq_at_least_30'] - df['precision'] < - 0.01
    return taxa_stat(ncbi, df, ix)


def taxa_that_get_wiped(ncbi, df):
    ix = np.isnan(df['precision_mapq_at_least_30'])
    return taxa_stat(ncbi, df, ix)



def do(input_db, same_what, refdb_ncbi, output_png):

    df = get_data(input_db, same_what, refdb_ncbi)

    fig, ax = plt.subplots()
# no legend!
#    colors = { "Metazoa": "yellow", "Viridiplantae": "blue", "Fungi": "green", "Other": "grey"}
#    df['kingdomColors'] = [colors[get_kingdom(ncbi, x)] for x in df['source_taxon']]
#    df.plot.scatter(x='fraction_mapq_at_least_30', y='precision', alpha=0.3,c='kingdomColors', ax=ax)

# annoyingly overplotted
    for name, group in [t for t in df.groupby("kingdom")]:
      plt.plot(group['precision'], group['precision_mapq_at_least_30'], marker="x", linestyle="", label=name, alpha=0.2)

    ax.legend()
    ax.set_xlabel("Precision")
    ax.set_ylabel("Precision when MAPQ >= 30")
    fig.savefig(output_png, bbox_inches='tight', dpi=199)
    
def get_kingdom(ncbi, taxid):
    lineage = ncbi.get_lineage(taxid)
    names = ncbi.get_taxid_translator(lineage)
    d = ncbi.get_rank(lineage)
    ks = [k for k in d if d[k] == 'kingdom']
    if len(ks) == 1:
      return names[ks[0]]
    if 2759 in d:
        return "Other"
    raise ValueError(d)


def opts(argv):
    parser = argparse.ArgumentParser(
      description="plot precision by species",
      formatter_class = argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input-alignments-sqlite", type=str, action="store", dest="input_db", required=True)
    parser.add_argument("--output-png", type=str, action="store", dest="output_png", required=True)
    parser.add_argument("--refdb-ncbi", type=str, action="store", dest="refdb_ncbi", help = "argument for ete.NCBITaxa")

    parser.add_argument("--aggregation-level", type=str, action="store", dest="same_what", required=True)
    return parser.parse_args(argv)

def main(argv=sys.argv[1:]):


    options = opts(argv)
    do(**vars(options))
    


if __name__ == '__main__':
    main()
