import pandas
import sqlite3
from sklearn.decomposition import PCA

sql = '''
SELECT pairs.source_taxon          AS taxon_a,
       pairs.matched_taxon         AS taxon_b,
       sources_a.num_source_reads  AS num_source_reads_total_a,
       self_a.num_self_matches     AS num_self_matches_a,
       matches_a.num_matched_reads AS num_matched_reads_total_a,
       pairs.num_pairs             AS num_reads_a_to_b,
       sources_b.num_source_reads  AS num_source_reads_total_b,
       self_b.num_self_matches     AS num_self_matches_b,
       matches_b.num_matched_reads AS num_matched_reads_total_b,
       pairs_inv.num_pairs         AS num_reads_b_to_a
FROM   (SELECT source_taxon,
               Count(*) AS num_source_reads
        FROM   alignment_from_known_source
        GROUP  BY source_taxon) sources_a,
       (SELECT source_taxon,
               Count(*) AS num_source_reads
        FROM   alignment_from_known_source
        GROUP  BY source_taxon) sources_b,
       (SELECT matched_taxon,
               Count(*) AS num_matched_reads
        FROM   alignment_from_known_source
        GROUP  BY matched_taxon) matches_a,
       (SELECT matched_taxon,
               Count(*) AS num_matched_reads
        FROM   alignment_from_known_source
        GROUP  BY matched_taxon) matches_b,
       (SELECT source_taxon,
               matched_taxon,
               Count(*) AS num_pairs
        FROM   alignment_from_known_source
        WHERE  source_taxon < matched_taxon
        GROUP  BY source_taxon,
                  matched_taxon) pairs,
       (SELECT source_taxon,
               matched_taxon,
               Count(*) AS num_pairs
        FROM   alignment_from_known_source
        WHERE  source_taxon > matched_taxon
        GROUP  BY source_taxon,
                  matched_taxon) pairs_inv,
       (SELECT source_taxon,
               Count(*) AS num_self_matches
        FROM   alignment_from_known_source
        WHERE  source_taxon = matched_taxon
        GROUP  BY source_taxon) self_a,
       (SELECT source_taxon,
               Count(*) AS num_self_matches
        FROM   alignment_from_known_source
        WHERE  source_taxon = matched_taxon
        GROUP  BY source_taxon) self_b
WHERE  sources_a.source_taxon = pairs.source_taxon
       AND matches_a.matched_taxon = pairs.source_taxon
       AND sources_b.source_taxon = pairs.matched_taxon
       AND matches_b.matched_taxon = pairs.matched_taxon
       AND pairs.source_taxon = pairs_inv.matched_taxon
       AND pairs.matched_taxon = pairs_inv.source_taxon
       AND pairs.source_taxon = self_a.source_taxon
       AND pairs.matched_taxon = self_b.source_taxon 

'''

df = pandas.read_sql_query(sql, sqlite3.connect('tmp/100.0.0.0.0.alignments.sqlite'));

df['rate_emit_cross_matches_a'] = df['num_reads_a_to_b'] / df['num_source_reads_total_a']
df['rate_emit_other_matches_a'] = (df['num_source_reads_total_a']  - df['num_self_matches_a'] - df['num_reads_a_to_b'] )  / df['num_source_reads_total_a']
df['rate_accept_cross_matches_a'] = df['num_reads_b_to_a'] / df['num_matched_reads_total_a']
df['rate_accept_other_matches_a'] = (df['num_matched_reads_total_a'] - df['num_self_matches_a'] - df['num_reads_b_to_a']) / df['num_matched_reads_total_a']
df['rate_emit_cross_matches_b'] = df['num_reads_b_to_a'] / df['num_source_reads_total_b']
df['rate_emit_other_matches_b'] = (df['num_source_reads_total_b']  - df['num_self_matches_b'] - df['num_reads_b_to_a'] )  / df['num_source_reads_total_b']
df['rate_accept_cross_matches_b'] = df['num_reads_a_to_b'] / df['num_matched_reads_total_b']
df['rate_accept_other_matches_b'] = (df['num_matched_reads_total_b'] - df['num_self_matches_b'] - df['num_reads_a_to_b']) / df['num_matched_reads_total_b']

df = df[['taxon_a', 'taxon_b', 'rate_emit_cross_matches_a', 'rate_emit_other_matches_a','rate_accept_cross_matches_a', 'rate_accept_other_matches_a', 'rate_emit_cross_matches_b', 'rate_emit_other_matches_b', 'rate_accept_cross_matches_b', 'rate_accept_other_matches_b']]
df = df.set_index(['taxon_a', 'taxon_b'])

pca = PCA(n_components=2)
pca.fit(df)


# first component: how much a and b accept other matches
# second component: how confusable a and b are with each other
"""
>>> pca.explained_variance_ratio_.tolist()
[0.6285058653211216, 0.1808318437373069]
>>> [x for x in zip(df.columns, pca.components_[0] )]
[('rate_emit_cross_matches_a', 0.009032944786850607), ('rate_emit_other_matches_a', 0.5083683425430616), ('rate_accept_cross_matches_a', 0.012117526198515582), ('rate_accept_other_matches_a', 0.5180661055012918), ('rate_emit_cross_matches_b', 0.008617479581191023), ('rate_emit_other_matches_b', 0.48132271985088476), ('rate_accept_cross_matches_b', 0.005649506383429148), ('rate_accept_other_matches_b', 0.49108346701623734)]
>>> [x for x in zip(df.columns, pca.components_[1] )]
[('rate_emit_cross_matches_a', 0.5157742464118173), ('rate_emit_other_matches_a', 0.06198509331246638), ('rate_accept_cross_matches_a', 0.5146651863043344), ('rate_accept_other_matches_a', 0.0773771132609782), ('rate_emit_cross_matches_b', 0.46972839845403724), ('rate_emit_other_matches_b', -0.09512043486213578), ('rate_accept_cross_matches_b', 0.4709117796519074), ('rate_accept_other_matches_b', -0.08841229347553647)]
"""

df[["pca_1", "pca_2"]] =  pca.transform(df)

# 
