import sys
import argparse
import matplotlib.pyplot as plt
import pandas

def main(argv=sys.argv[1:]):
    parser = argparse.ArgumentParser(
      description="plot over mutation rate",
      formatter_class = argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input-json", type=str, action="store", dest="data_path", required=True)
    parser.add_argument("--aggregation-level", type=str, action="store", dest="same_what", required=True)
    parser.add_argument("--output-png", type=str, action="store", dest="output_path", required=True)
    options=parser.parse_args(argv)


    df = pandas.io.json.read_json(options.data_path)

#    df['mutation_rate'] = ["{:.3f}".format(x) for x in df['mutation_rate']]
    
    ax = df.plot(x="mutation_rate", y="precision" + options.same_what + "WhenMapqAtLeast30Queries", label = "Precision, MAPQ>=30 only", linestyle='-', marker='o')
    ax = df.plot(x="mutation_rate", y="precision" + options.same_what + "Queries", ax = ax, label = "Precision", linestyle='-', marker='o')
    ax = df.plot(x="mutation_rate", y="recall" + options.same_what + "Queries", ax = ax, label = "Recall", linestyle='-', marker='o')
    ax = df.plot(x="mutation_rate", y="recall" + options.same_what + "WhenMapqAtLeast30Queries", ax = ax, label = "Recall, MAPQ>=30 only", linestyle='-', marker='o')
    ax = df.plot(x="mutation_rate", y="mapqsFractionAtLeast30", ax = ax, label = "Fraction MAPQ>=30", linestyle='-', marker='o')
    import matplotlib.ticker as mtick
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
#    df.plot(x="mutation_rate", y="mapqsFractionAtLeast30", xlim=(0,0.2))

    ax.set_xlabel("Mutation rate")

    fig = ax.get_figure()

    fig.savefig(options.output_path, bbox_inches='tight', dpi=199)

if __name__ == '__main__':
    main()
