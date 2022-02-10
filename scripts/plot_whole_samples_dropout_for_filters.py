import sys
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet
from matplotlib import cm

report_classes = {
  "NR": "No results",
  "OC": "Single result",
  "MC": "Many results",
  "OI": "Single result",
  "MI": "Many results",
}

def make_series(data_path):
    M = pd.read_csv(data_path, sep = "\t", comment = "#", index_col = 0)
    nrows = len(M.index)
    I = M.stack().index.to_series()
    S = M.stack().apply(lambda x: report_classes[x])
    df = pd.DataFrame(data = {"Results": S})
    is_m2 = "Two markers"
    is_r4 = "Four reads"
    is_M30 = "MAPQ >= 30"
    is_Euk = "EukDetect"
    is_EukMod = "mod. EukDetect"
    is_our = "Our method"
    df[is_m2] = I.apply(lambda x : "m2" in x[1] or "EukDetect" in x[1] or "our method" in x[1])
    df[is_r4] = I.apply(lambda x : "r4" in x[1] or "EukDetect" in x[1])
    df[is_M30] = I.apply(lambda x : "M30" in x[1] or "EukDetect" == x[1])
    df[is_Euk] = I.apply(lambda x : "EukDetect" == x[1])
    df[is_EukMod] = I.apply(lambda x: "EukDetect " in x[1])
    df[is_our] = I.apply(lambda x: "our method" == x[1])
    return df.set_index([is_m2, is_r4, is_M30, is_Euk, is_EukMod, is_our])

def main(argv=sys.argv[1:]):
    parser = argparse.ArgumentParser(
      description="plot over mutation rate",
      formatter_class = argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input-tsv", type=str, action="store", dest="data_path", required=True)
    parser.add_argument("--output-png", type=str, action="store", dest="output_path", required=True)
    options=parser.parse_args(argv)
    df = make_series(options.data_path)

    u = UpSet(df, intersection_plot_elements=0, sort_by = None)
    u.add_stacked_bars(by="Results", colors=cm.Pastel1,
                   title="Count", elements=10)

    fig, ax = plt.subplots()

    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

    d = u.plot(fig)
    fig.delaxes(d['totals'])
    d['extra0'].legend(bbox_to_anchor=(-0.1, 1.0))
#    d['intersections'].set_ylim([0,1])
#    d['intersections'].set_ylabel("Samples with no results")
    fig.savefig(options.output_path, bbox_inches='tight', dpi=199)
    #from PIL import Image
    #im = Image.open(options.output_path)
    #h, w  = im.size
    #im_crop = im.crop((70, 0, w - 70, h))
    #im_crop.save(options.output_path)


if __name__ == '__main__':
    main()
