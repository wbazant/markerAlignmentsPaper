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

# https://gis.stackexchange.com/questions/436908/selecting-n-samples-uniformly-from-a-grid-points
def pick_subsample(df, initial_taxa, DISTANCE_MIN):

    names_ok = [(a,b) for a, b in initial_taxa]
    values_ok = [(pca_1, pca_2) for pca_1, pca_2 in zip(df.loc[initial_taxa]["pca_1"], df.loc[initial_taxa]["pca_2"])] 

    for index, row in df.iterrows():
        x,y = row[["pca_1", "pca_2"]]
        if any((euclidean((x,y), point_ok) < DISTANCE_MIN for point_ok in values_ok)):
            continue
        else:
            names_ok.append(index)
            values_ok.append((x,y))
    return names_ok

def get_df_for_pca(sqlite_path):
    df = get_df(sqlite_path)
    df = df.set_index(['taxon_a', 'taxon_b'])
    # if both (x,y) and (y,x) are in the dataset, pseudorandomly pick one of them
    xs = [hashlib.md5(x.encode()).hexdigest() < hashlib.md5(y.encode()).hexdigest() or (y,x) not in df.index for x,y in df.index]
    df = df.loc[xs]
    return df

def is_both_species_column(refdb_markers, df):
    x = set(pandas.read_csv(refdb_markers)['Taxonomy_ID'])
    return df.reset_index().apply(lambda r: int(r['taxon_a']) in x and int(r['taxon_b']) in x, axis = 1).to_list()

def do(refdb_ncbi, refdb_markers, sqlite_path, subsample_distance, second_subsample_distance, second_subsample_size, third_subsample_distance, third_subsample_size, fourth_subsample_size, output_tsv, logger):
    df = get_df_for_pca(sqlite_path)
    data_columns = df.columns.to_list()
    logger.debug("Num data points: %s", len(df))
    df[["pca_1", "pca_2"]] = pca_components(df, logger)

    df['is_both_species'] = is_both_species_column(refdb_markers, df)

    logger.debug("Num data points where both taxa are species level: %s", len(df[df['is_both_species']]))


    subsample = pick_subsample(df[df['is_both_species']],
            initial_taxa = [("370354", "294381")], # E. dispar, E. hystolytica
            DISTANCE_MIN = subsample_distance) 

    logger.debug("Num subsampled points: %s", len(subsample))

    second_subsample = pick_subsample(df[df['is_both_species']],
            initial_taxa = subsample,
            DISTANCE_MIN = second_subsample_distance) 
    logger.debug("Num subsampled points for a second time: %s, will keep top %s with lowest PCA1", len(second_subsample), second_subsample_size)
    second_subsample_top_only = df[df.index.map(lambda n: n in second_subsample)].sort_values(by="pca_1").head(second_subsample_size).index.to_list()

    third_subsample = pick_subsample(df[df['is_both_species']],
            initial_taxa = second_subsample,
            DISTANCE_MIN = third_subsample_distance) 
    logger.debug("Num subsampled points for a third time: %s, will keep top %s with lowest PCA1", len(third_subsample), third_subsample_size)

    xs = df[df.index.map(lambda n: n in third_subsample)]

    eukdetect_paper_example = ("370354", "294381")
    eukdetect_x,eukdetect_y = df.loc[eukdetect_paper_example][["pca_1", "pca_2"]]
    logger.debug("Eukdetect example, PCA1 = %s, PCA2 = %s", eukdetect_x,eukdetect_y)

    df['distance_from_example'] = [euclidean((x,y), (eukdetect_x,eukdetect_y)) for x,y in zip(df["pca_1"], df["pca_2"])]
    
    third_subsample_top_only = df[df.index.map(lambda n: n in third_subsample)].sort_values(by="pca_1").head(third_subsample_size).index.to_list()
    third_subsample_nearest_example_only = df[df.index.map(lambda n: n in third_subsample)].sort_values(by="distance_from_example").head(fourth_subsample_size).index.to_list()

    df["is_picked_for_representative_subset"] =  df.index.map(lambda n: n in subsample)
#    df["is_picked_for_focused_subset"] =  df.index.map(lambda n: n in second_subsample_top_only)
#    df["is_picked_for_top_subset"] =  df.index.map(lambda n: n in third_subsample_top_only)
    df["is_picked_for_star_subset"] =  df.index.map(lambda n: n in third_subsample_nearest_example_only)
    df = df.reset_index()
    ncbi = NCBITaxa2(refdb_ncbi)
    t = ncbi.get_taxid_translator(set(df["taxon_a"].values) | set(df["taxon_b"].values))
    df["taxon_a_name"] = df["taxon_a"].map(lambda x: t[int(x)])
    df["taxon_b_name"] = df["taxon_b"].map(lambda x: t[int(x)])

    lca_columns = ["lca_taxon", "lca_taxon_name", "lca_rank", "lca_genus"]
    df[lca_columns] = [ncbi.lca_info([taxon_a, taxon_b]) for taxon_a, taxon_b in zip(df["taxon_a"], df["taxon_b"])]

    df = df.sort_values(by="pca_1", ascending = False)
    df = df[["taxon_a", "taxon_a_name", "taxon_b", "taxon_b_name", "is_picked_for_representative_subset", "is_picked_for_star_subset", "is_both_species"] + lca_columns + ["pca_1", "pca_2"] + data_columns ]
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
    parser.add_argument("--second-subsample-pca-distance", type=float, action = "store", dest="second_subsample_distance", default = 0.05 / 3)
    parser.add_argument("--second-subsample-size", type=int, action = "store", dest="second_subsample_size", default = 100)
    parser.add_argument("--third-subsample-pca-distance", type=float, action = "store", dest="third_subsample_distance", default = 0.05 / 9)
    parser.add_argument("--third-subsample-size", type=int, action = "store", dest="third_subsample_size", default = 100)
    parser.add_argument("--fourth-subsample-size", type=int, action = "store", dest="fourth_subsample_size", default = 50)
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


