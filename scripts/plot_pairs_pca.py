import pandas
import argparse
import sys
import matplotlib.pyplot as plt

#ABOLG: has A, has B, has others, good LCA, good genus
#ABOlG: has A, has B, has others, no good LCA, good genus
#ABOLg: has A, has B, has others, good LCA, no good genus
#ABOlg: has A, has B, has others, no good LCA, no good genus
#ABOL: has A, has B, has others, good LCA, n/a genus
#ABOl: has A, has B, has others, no good LCA, n/a genus
#ABo: has A, has B, no others
#AbO: has A, misses B, has others
#Abo: has A, misses B, no others
#aBO: misses A, has B, has others
#aBo: misses A, has B, no others
#abO: misses A, misses B, has others
#NR: No results
colors = {
  "ABOLG": "#40E0D0", # turquoise
  "ABOlG": "#A7C7E7", # pastel blue
  "ABOLg": "#40E0D0", # turquoise
  "ABOlg": "#954535", # chestnut
  "ABOL": "#40E0D0",  # turquoise
  "ABOl": "#F88379",  # coral pink
  "ABo": "blue",
  "AbO": "#FAFA33", # lemon yellow
  "Abo": "#FFE5B4", # peach
  "aBO": "#FAFA33", # lemon yellow
  "aBo": "#C9CC3F", # pear
  "abO": "#FAFA33", # lemon yellow
  "NR": "grey",
}
def cat_to_color(x):
    return colors[x] if x in colors else "grey" # 4 XXX datapoints at 0.01 where we failed to draw any reads


def subplot(df, column, label, ax):
    df.plot.scatter("pca_1", "pca_2", c = df[column].map(cat_to_color), label = label, ax = ax)
    ax.axis('off')
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.3, box.width, box.height*0.7])
    ax.legend(loc='upper center', edgecolor = "white", markerscale = 0, bbox_to_anchor=(0.5, -0.05))

def do(pca_tsv, results_tsv, output_png, output_xlsx):
    pca_df = pandas.read_csv(pca_tsv, sep = "\t").astype({"taxon_a": str, "taxon_b": str})
    results_df = pandas.read_csv(results_tsv, sep = "\t")

    results_df['taxon_a'] =  results_df['taxon_a'].map(lambda x: x.split("|")[0])
    results_df['taxon_b'] =  results_df['taxon_b'].map(lambda x: x.split("|")[0])

    df = pca_df.merge(results_df, on = ["taxon_a", "taxon_b"])

    fig, axs = plt.subplots(3, 2)
    subplot(df, "EukDetectLo", "EukDetect, 0.01 cov.", axs[0][0])
    subplot(df, "EukDetectMid", "EukDetect, 0.05 cov.", axs[1][0])
    subplot(df, "EukDetectHi", "EukDetect, 0.10 cov.", axs[2][0])
    subplot(df, "CORRALLo", "CORRAL, 0.01 cov.", axs[0][1])
    subplot(df, "CORRALMid", "CORRAL, 0.05 cov.", axs[1][1])
    subplot(df, "CORRALHi", "CORRAL, 0.10 cov.", axs[2][1])

    fig.savefig(output_png, bbox_inches='tight', dpi=199)

    df = pca_df.merge(results_df, on = ["taxon_a", "taxon_b"], how='left')
    df.fillna("")
    cs = df.columns
    df = df.set_index(['taxon_a_name', 'taxon_b_name'])
    
    df.to_excel(output_xlsx, sheet_name = "Pairs PCA and results",  float_format="%.3f")


def opts(argv):
    parser = argparse.ArgumentParser(
      description="plot pairs pca",
      formatter_class = argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--pca-tsv", type=str, action="store", dest="pca_tsv", required=True)
    parser.add_argument("--results-tsv", type=str, action="store", dest="results_tsv", required=True)
    parser.add_argument("--output-png", type=str, action="store", dest="output_png", required=True)
    parser.add_argument("--output-xlsx", type=str, action="store", dest="output_xlsx", required=True)
    return parser.parse_args(argv)

def main(argv=sys.argv[1:]):
    options = opts(argv)
    kwargs = vars(options)
    do(**kwargs)

if __name__ == '__main__':
    main()
