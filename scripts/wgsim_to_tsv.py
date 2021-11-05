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
    parser.add_argument("--output-tsv", type=str, action="store", dest="output_path", required=True)
    options=parser.parse_args(argv)


    df = pandas.io.json.read_json(options.data_path)

    df['mutation_rate'] = ["{:.3f}".format(x) for x in df['mutation_rate']]
    
    df.set_index('mutation_rate', inplace=True)
    df.to_csv(options.output_path, sep = "\t")

if __name__ == '__main__':
    main()
