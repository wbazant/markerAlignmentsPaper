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

    df['mutation_rate'] = ["{:.3f}".format(x) for x in df['mutation_rate']]
    

    df = pandas.DataFrame({'numQueriesTotal': df['numQueriesTotal'], 'numQueriesMapped': df['numQueriesMapped'], 'numQueriesMappedOnlyAsSameFamily': df['numQueriesMappedOnlyAsSameFamily'], 'numQueriesMappedOnlyAsSameGenus': df['numQueriesMappedOnlyAsSameGenus'], 'numQueriesMappedOnlyAsSameSpecies': df['numQueriesMappedOnlyAsSameSpecies'],
        'numQueriesMappedWhenMapqAtLeast30': df['numQueriesMappedWhenMapqAtLeast30'],
        'numQueriesMappedOnlyAsSameFamilyWhenMapqAtLeast30': df['numQueriesMappedOnlyAsSameFamilyWhenMapqAtLeast30'], 'numQueriesMappedOnlyAsSameGenusWhenMapqAtLeast30': df['numQueriesMappedOnlyAsSameGenusWhenMapqAtLeast30'], 'numQueriesMappedOnlyAsSameSpeciesWhenMapqAtLeast30': df['numQueriesMappedOnlyAsSameSpeciesWhenMapqAtLeast30'],
        'mutation_rate' : df['mutation_rate'] })

    # skip two points so the distances are even in mutation rate
    #df = df.drop(axis="index", labels=[1,2])

    ax = df.plot.bar(x='mutation_rate')

    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.3),
          ncol=2, fancybox=True, shadow=True)

    fig = ax.get_figure()

    fig.savefig(options.output_path, bbox_inches='tight', dpi=199)

if __name__ == '__main__':
    main()
