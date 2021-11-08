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
    

    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.tight_layout()

    df1 = pandas.DataFrame({'numQueriesTotal': df['numQueriesTotal'], 'numQueriesMapped': df['numQueriesMapped'], 'numQueriesMappedOnlyAsSameFamily': df['numQueriesMappedOnlyAsSameFamily'], 'numQueriesMappedOnlyAsSameGenus': df['numQueriesMappedOnlyAsSameGenus'], 'numQueriesMappedOnlyAsSameSpecies': df['numQueriesMappedOnlyAsSameSpecies'], 'mutation_rate' : df['mutation_rate'] })

    # skip two points so the distances are even in mutation rate
    #df = df.drop(axis="index", labels=[1,2])
    
    df1.plot.bar(x='mutation_rate', ax=ax1)

    df2 = pandas.DataFrame({'numQueriesTotal': df['numQueriesTotal'],
        'numQueriesMapped': df['numQueriesMapped'],
        'numQueriesMappedWhenMapqAtLeast30': df['numQueriesMappedWhenMapqAtLeast30'],
        'numQueriesMappedOnlyAsSameFamilyWhenMapqAtLeast30': df['numQueriesMappedOnlyAsSameFamilyWhenMapqAtLeast30'], 'numQueriesMappedOnlyAsSameGenusWhenMapqAtLeast30': df['numQueriesMappedOnlyAsSameGenusWhenMapqAtLeast30'], 'numQueriesMappedOnlyAsSameSpeciesWhenMapqAtLeast30': df['numQueriesMappedOnlyAsSameSpeciesWhenMapqAtLeast30'],
        'mutation_rate' : df['mutation_rate'] })
    df2.plot.bar(x='mutation_rate', ax=ax2)

    fig.savefig(options.output_path, bbox_inches='tight')

if __name__ == '__main__':
    main()
