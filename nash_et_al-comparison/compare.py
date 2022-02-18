import sys
import os
sys.path.append(os.path.join(os.path.realpath(os.path.dirname(__file__)), "../"))
from lib.ncbi2 import NCBITaxa2
import pandas
import argparse

import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.api import add_constant

def regress_value_nash_and_value_us(df):
    return sm.OLS(df["value_us"], df["value_nash"]).fit()

def read_nash_xlsx_as_series(nash_xlsx):
    df = pandas.read_excel(nash_xlsx, sheet_name = "Table S3 - Fungal Reads", engine = "openpyxl")
    df = df.rename({"SRA_Sample": "Sample", "Name": "Species"}, axis = 1)
    df = df.groupby(["Sample", "Species"]).size()
    df.name = 'value'
    df = df.astype(float)
    return df

def top_hits(df):
    res = pandas.DataFrame()
    for name, group in df.groupby("Sample"):
        m = group["value"].max()
        ms = group.loc[group["value"] == m]
        res = res.append(ms)
    return res



def read_our_triples_as_series(our_triples_tsv):
    df = pandas.read_csv(our_triples_tsv, sep = "\t", names = ["Sample", "Species", "value"])
    df = df.set_index(["Sample", "Species"])
    return df["value"]

def do(nash_xlsx, our_triples_tsv, refdb_ncbi):
    ncbi = NCBITaxa2(refdb_ncbi)
    nash_df = read_nash_xlsx_as_series(nash_xlsx).to_frame().reset_index()
    our_df = read_our_triples_as_series(our_triples_tsv).to_frame().reset_index()

    #"Fungi" in ncbi.get_taxid_translator(ncbi.get_lineage(76773)).values() 
    nash_df["Taxon_id"] = nash_df["Species"].apply(ncbi.get_taxid_from_string)
    our_df["Taxon_id"] = our_df["Species"].apply(ncbi.get_taxid_from_string)

    our_fungi_bv = our_df["Taxon_id"].apply(lambda x: "Fungi" in ncbi.get_taxid_translator(ncbi.get_lineage(x)).values())
    our_df = our_df.loc[our_fungi_bv]

    nash_res = top_hits(nash_df)
    our_res = top_hits(our_df)
    all_res = nash_res.merge(our_res, on = "Sample", how = "inner", suffixes = ("_nash", "_us"))
    all_res["Match_type"] = all_res.apply(lambda x: ncbi.get_match_type([x.Taxon_id_nash, x.Taxon_id_us]), axis = 1)
    all_res["Match_parent"] = all_res.apply(lambda x: ncbi.get_common_parent([x.Taxon_id_nash, x.Taxon_id_us]), axis = 1)


    print("Nash species: " + repr(len(set(nash_df["Species"].values))))
    print("Nash samples: " + repr(len(set(nash_df["Sample"].values))))
    print("Nash data points: " + repr(len(nash_df)))

    print("Our species: " + repr(len(set(our_df["Species"].values))))
    print("Our samples: " + repr(len(set(our_df["Sample"].values))))
    print("Our data points: " + repr(len(our_df)))


    print("Match types of top hits: \n" + repr(all_res[["Match_type","Match_parent"]].value_counts()))

    df = all_res.loc[(all_res["Match_type"] == "species") & (all_res["Match_parent"] != "Schizosaccharomyces pombe") ].sort_values("value_us").reset_index()
    m = regress_value_nash_and_value_us(df)
    print("Regression on matching species - Nash vs us: \n" + repr(m.summary()))
    return nash_df, our_df, nash_res, our_res, ncbi, all_res
# nash_df, our_df, nash_res, our_res, ncbi, all_res = do("nash_et_al.xlsx", "taxon_num_reads.triples.tsv", "/home/wbazant/dev/markerAlignmentsPaper/refdb/taxa.sqlite")
    # top_hits_per_sample = res.sort_values(by=0, ascending = False).value_counts("Name")



def main(argv=sys.argv[1:]):
    parser = argparse.ArgumentParser(
      description="compare to Nash et al.",
      formatter_class = argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--nashXlsx", type=str, action="store", dest="nash_xlsx", help = "Spreadsheet from the paper", required = True)
    parser.add_argument("--ourTriplesTsv", type=str, action="store", dest="our_triples_tsv", help = "our triples tsv", required = True)
    parser.add_argument("--refdb-ncbi", type=str, action="store", dest="refdb_ncbi", help = "argument for ete.NCBITaxa", required = True)



    options=parser.parse_args(argv)
    do(**vars(options))


if __name__ == '__main__':
    main()
