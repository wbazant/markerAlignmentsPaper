import sys
import argparse
import matplotlib.pyplot as plt
import pandas

def main(argv=sys.argv[1:]):
    parser = argparse.ArgumentParser(
      description="plot_simulations",
      formatter_class = argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input-json", type=str, action="store", dest="data_path", required=True)
    parser.add_argument("--output-png", type=str, action="store", dest="output_path", required=True)
    options=parser.parse_args(argv)


    df = pandas.io.json.read_json(options.data_path)

    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.tight_layout()

    df.plot(x="mutation_rate", y="recallQueries", ylim=(0,1), xlim=(0,0.2), kind ="scatter", ax=ax1)
    df.plot(x="mutation_rate", y="mapqsFractionAtLeast30", ylim=(0,1), xlim=(0,0.2), kind ="scatter", ax=ax2)



    plt.savefig(options.output_path, bbox_inches='tight')

if __name__ == '__main__':
    main()
