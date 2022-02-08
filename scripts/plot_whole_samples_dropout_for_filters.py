import sys
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet


def main(argv=sys.argv[1:]):
    parser = argparse.ArgumentParser(
      description="plot over mutation rate",
      formatter_class = argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input-tsv", type=str, action="store", dest="data_path", required=True)
    parser.add_argument("--output-png", type=str, action="store", dest="output_path", required=True)
    options=parser.parse_args(argv)
    M = pd.read_csv(options.data_path, sep = "\t", comment = "#", index_col = 0)
    nrows = len(M.index)
    df = M.apply(pd.value_counts).transpose()

    is_m2 = "two markers"
    is_r4 = "four reads"
    is_M30 = "MAPQ >= 30"
    df[is_m2] = df.index.to_series().str.contains("m2") | df.index.to_series().str.contains("EukDetect")
    df[is_r4] = df.index.to_series().str.contains("r4") | df.index.to_series().str.contains("EukDetect")
    df[is_M30] = df.index.to_series().str.contains("M30") | df.index.to_series().str.endswith("EukDetect")

    df = df.head(8)
    df["NRf"] = df["NR"] / nrows

    cs = pd.Series(df["NRf"].array, index=pd.MultiIndex.from_frame(df[[is_m2, is_r4, is_M30]]))

    u = UpSet(cs, sort_by="cardinality")

    fig, ax = plt.subplots()

    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

    d = u.plot(fig)
    fig.delaxes(d['totals'])
    d['intersections'].set_ylim([0,1])
    d['intersections'].set_ylabel("Samples with no results")
    fig.savefig(options.output_path, bbox_inches='tight', dpi=199)

if __name__ == '__main__':
    main()
