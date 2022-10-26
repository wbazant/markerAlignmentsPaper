import pandas
import sqlite3
from sklearn.decomposition import PCA
import hashlib
from scipy.spatial.distance import euclidean
import os
import sys
sys.path.append(os.path.join(os.path.realpath(os.path.dirname(__file__)), "../"))
from lib.ncbi2 import NCBITaxa2
import argparse
import logging

sql = '''
SELECT pairs.taxon_a,
       pairs.taxon_b,
       sources_a.num_source_reads  AS num_source_reads_total_a,
       self_a.num_self_matches     AS num_self_matches_a,
       matches_a.num_matched_reads AS num_matched_reads_total_a,
       pairs.num_pairs             AS num_reads_a_to_b,
       sources_b.num_source_reads  AS num_source_reads_total_b,
       self_b.num_self_matches     AS num_self_matches_b,
       matches_b.num_matched_reads AS num_matched_reads_total_b,
       ifnull(pairs_inv.num_pairs, 0)         AS num_reads_b_to_a
FROM   
       (SELECT source_taxon as taxon_a,
               matched_taxon as taxon_b,
               Count(*) AS num_pairs
        FROM   alignment_from_known_source
        WHERE  source_taxon != matched_taxon
        GROUP  BY source_taxon,
                  matched_taxon) pairs
join
       (SELECT source_taxon as taxon,
               Count(*) AS num_source_reads
        FROM   alignment_from_known_source
        GROUP  BY source_taxon) sources_a
  on pairs.taxon_a = sources_a.taxon        
join
       (SELECT source_taxon as taxon,
               Count(*) AS num_source_reads
        FROM   alignment_from_known_source
        GROUP  BY source_taxon) sources_b
  on pairs.taxon_b = sources_b.taxon
join
       (SELECT matched_taxon as taxon,
               Count(*) AS num_matched_reads
        FROM   alignment_from_known_source
        GROUP  BY matched_taxon) matches_a
  on pairs.taxon_a = matches_a.taxon
join
       (SELECT matched_taxon as taxon,
               Count(*) AS num_matched_reads
        FROM   alignment_from_known_source
        GROUP  BY matched_taxon) matches_b
  on pairs.taxon_b = matches_b.taxon
join
       (SELECT source_taxon as taxon,
               Count(*) AS num_self_matches
        FROM   alignment_from_known_source
        WHERE  source_taxon = matched_taxon
        GROUP  BY source_taxon) self_a
  on pairs.taxon_a = self_a.taxon       
join
       (SELECT source_taxon as taxon,
               Count(*) AS num_self_matches
        FROM   alignment_from_known_source
        WHERE  source_taxon = matched_taxon
        GROUP  BY source_taxon) self_b
  on pairs.taxon_b = self_b.taxon
left join
       (SELECT source_taxon as taxon_a,
               matched_taxon as taxon_b,
               Count(*) AS num_pairs
        FROM   alignment_from_known_source
        WHERE  source_taxon != matched_taxon
        GROUP  BY source_taxon,
                  matched_taxon) pairs_inv
on (pairs.taxon_a = pairs_inv.taxon_b and pairs.taxon_b = pairs_inv.taxon_a)
'''

sqlite_path = 'tmp/100.0.0.0.0.alignments.sqlite'
def get_df(sqlite_path):
    df = pandas.read_sql_query(sql, sqlite3.connect(sqlite_path));
    df['rate_emit_cross_matches_a'] = df['num_reads_a_to_b'] / df['num_source_reads_total_a']
    df['rate_emit_other_matches_a'] = (df['num_source_reads_total_a']  - df['num_self_matches_a'] - df['num_reads_a_to_b'] )  / df['num_source_reads_total_a']
    df['rate_accept_cross_matches_a'] = df['num_reads_b_to_a'] / df['num_matched_reads_total_a']
    df['rate_accept_other_matches_a'] = (df['num_matched_reads_total_a'] - df['num_self_matches_a'] - df['num_reads_b_to_a']) / df['num_matched_reads_total_a']
    df['rate_emit_cross_matches_b'] = df['num_reads_b_to_a'] / df['num_source_reads_total_b']
    df['rate_emit_other_matches_b'] = (df['num_source_reads_total_b']  - df['num_self_matches_b'] - df['num_reads_b_to_a'] )  / df['num_source_reads_total_b']
    df['rate_accept_cross_matches_b'] = df['num_reads_a_to_b'] / df['num_matched_reads_total_b']
    df['rate_accept_other_matches_b'] = (df['num_matched_reads_total_b'] - df['num_self_matches_b'] - df['num_reads_a_to_b']) / df['num_matched_reads_total_b']

    df = df[['taxon_a', 'taxon_b', 'rate_emit_cross_matches_a', 'rate_emit_other_matches_a','rate_accept_cross_matches_a', 'rate_accept_other_matches_a', 'rate_emit_cross_matches_b', 'rate_emit_other_matches_b', 'rate_accept_cross_matches_b', 'rate_accept_other_matches_b']]
    return df

    
def pca_components(df, logger):
    pca = PCA(n_components=2)
    logger.debug("Fitting PCA")
    pca.fit(df)
    logger.debug("Done. Explained variance: %s", pca.explained_variance_ratio_)
    logger.debug("PCA1 = " + " + ".join([str("%.04f"%x) + " * " + str(y) for x, y in zip(pca.components_[0], df.columns)]))
    logger.debug("PCA2 = " + " + ".join([str("%.04f"%x) + " * " + str(y) for x, y in zip(pca.components_[1], df.columns)]))
    return pca.transform(df)

def pick_subsample(df, DISTANCE_MIN):
    # E. dispar, E. hystolytica
    eukdetect_paper_example = ("370354", "294381")

    eukdetect_x,eukdetect_y = df.loc[eukdetect_paper_example][["pca_1", "pca_2"]]

    # https://gis.stackexchange.com/questions/436908/selecting-n-samples-uniformly-from-a-grid-points
    names_ok = [eukdetect_paper_example]
    list_ok = [(eukdetect_x,eukdetect_y)]

    for index, row in df.iterrows():
        x,y = row[["pca_1", "pca_2"]]
        if any((euclidean((x,y), point_ok) < DISTANCE_MIN for point_ok in list_ok)):
            continue
        else:
            names_ok.append(index)
            list_ok.append((x,y))
    return names_ok

def get_df_for_pca(sqlite_path):
    df = get_df(sqlite_path)
    df = df.set_index(['taxon_a', 'taxon_b'])
    # if both (x,y) and (y,x) are in the dataset, pseudorandomly pick one of them
    xs = [hashlib.md5(x.encode()).hexdigest() < hashlib.md5(y.encode()).hexdigest() or (y,x) not in df.index for x,y in df.index]
    df = df.loc[xs]
    return df

def do(refdb_ncbi, refdb_markers, sqlite_path, subsample_distance, output_tsv, logger):
    df = get_df_for_pca(sqlite_path)
    data_columns = df.columns.to_list()
    logger.debug("Num data points: %s", len(df))
    df[["pca_1", "pca_2"]] = pca_components(df, logger)

    df = df.reset_index()
    x = set(pandas.read_csv(refdb_markers)['Taxonomy_ID'])
    df['is_both_species'] = df.apply(lambda r: int(r['taxon_a']) in x and int(r['taxon_b']) in x, axis = 1)
    df = df.set_index(['taxon_a', 'taxon_b'])

    logger.debug("Num data points where both taxa are species level: %s", len(df[df['is_both_species']]))

    subsample = pick_subsample(df[df['is_both_species']], DISTANCE_MIN = subsample_distance) 
    logger.debug("Num subsampled points: %s", len(subsample))
    df["is_picked"] =  df.index.map(lambda n: n in subsample)
    df = df.reset_index()
    ncbi = NCBITaxa2(refdb_ncbi)
    t = ncbi.get_taxid_translator(set(df["taxon_a"].values) | set(df["taxon_b"].values))
    df["taxon_a_name"] = df["taxon_a"].map(lambda x: t[int(x)])
    df["taxon_b_name"] = df["taxon_b"].map(lambda x: t[int(x)])

    lca_columns = ["lca_taxon", "lca_taxon_name", "lca_rank", "lca_genus"]
    df[lca_columns] = [ncbi.lca_info([taxon_a, taxon_b]) for taxon_a, taxon_b in zip(df["taxon_a"], df["taxon_b"])]

    df = df.sort_values(by="pca_1", ascending = False)
    df = df[["taxon_a", "taxon_a_name", "taxon_b", "taxon_b_name", "is_picked", "is_both_species"] + lca_columns + ["pca_1", "pca_2"] + data_columns ]
    df.to_csv(output_tsv, index = False, sep = "\t", float_format='%.3f')


def opts(argv):
    parser = argparse.ArgumentParser(
      description="make diabimmune comparison spreadsheet",
      formatter_class = argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--refdb-ncbi", type=str, action="store", dest="refdb_ncbi", help = "argument for ete.NCBITaxa", required = True)
    parser.add_argument("--refdb-markers", type=str, action="store", dest="refdb_markers", help = "marker genes per species csv", required = True)
    parser.add_argument("--input-sqlite", type=str, action="store", dest="sqlite_path", required=True)
    parser.add_argument("--subsample-pca-distance", type=float, action = "store", dest="subsample_distance", default = 0.05)
    parser.add_argument("--verbose", action="store_true", dest="verbose")
    parser.add_argument("--output-tsv", type=str, action="store", dest="output_tsv", required=True)
    return parser.parse_args(argv)

def main(argv=sys.argv[1:]):
    options = opts(argv)
    logging.basicConfig(format='%(asctime)s: %(message)s')
    logger=logging.getLogger(__name__)
    if options.verbose:
        logger.setLevel(logging.DEBUG)
    kwargs = vars(options)
    del kwargs['verbose']
    do(**kwargs, logger = logger)

if __name__ == '__main__':
    main()


