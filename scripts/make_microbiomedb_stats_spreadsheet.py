import sys
import os
import argparse
import xlsxwriter
import pandas


#https://stackoverflow.com/a/40535454
def write_and_fix_column_width(writer, df, sheet_name, *args, **kwargs):
    df = df.reset_index()
    df.to_excel(writer, sheet_name = sheet_name, index = None, *args, **kwargs)
    worksheet = writer.sheets[sheet_name]
    #Iterate through each column and set the width == the max length in that column. A padding length of 2 is also added.
    for i, col in enumerate(df.columns):
        # find length of column i
        column_len = df[col].astype(str).str.len().max()
        # Setting the length if the column header is larger
        # than the max column value length
        column_len = max(column_len, len(col)) + 2
        # set the column length
        worksheet.set_column(i, i, column_len)

def do(results_dir, name_postfix, output_xlsx):
    dfs = read_dfs(results_dir, name_postfix)
    df, stats_df = read_dfs(results_dir, name_postfix)
    with pandas.ExcelWriter(output_xlsx, engine="xlsxwriter") as writer:
        write_and_fix_column_width(writer, stats_df, sheet_name = "Summary statistics")

        df_vc_genus = df[df["genus"] != ""][["study", "sample", "kingdom", "genus"]].drop_duplicates().value_counts(["kingdom", "genus"]).to_frame('counts_genus')
        write_and_fix_column_width(writer, df_vc_genus, sheet_name = "Value counts - genus detected")

        
        df_vc_taxon = df.value_counts(["kingdom", "genus", "taxon"]).to_frame('counts_taxon')
        write_and_fix_column_width(writer, df_vc_taxon, sheet_name = "Value counts - taxon detected")

        df_vc_genus_study = df[df["genus"] != ""][["study", "sample", "kingdom", "genus"]].drop_duplicates().value_counts(["study", "kingdom", "genus"]).to_frame('counts_genus_study')
        write_and_fix_column_width(writer, df_vc_genus_study, sheet_name = "Value counts - genus, by study") 

        df_vc_taxon_study = df.value_counts(["study", "kingdom", "genus", "taxon"]).to_frame('counts_taxon_study')
        write_and_fix_column_width(writer, df_vc_taxon_study, sheet_name = "Value counts - taxon, by study")

        df_all = df[["study", "sample", "kingdom", "genus", "taxon", "value"]]
        df_all = df_all.set_index(["study", "sample", "kingdom", "genus", "taxon",])
        write_and_fix_column_width(writer, df_all, sheet_name = "All values - CPM", float_format="%.6f")


def read_dfs(results_dir, name_postfix):
    result_files = [ f for f in os.listdir(results_dir) if f.endswith(name_postfix)] 
    dfs = []
    total_num_studies = len(result_files)
    total_num_samples = 0
    for result_file in result_files:
        study = result_file.replace(name_postfix, "")
        df = pandas.read_csv(results_dir + "/" + result_file, sep = "\t")
        total_num_samples += len([c for c in df.columns if c != "lineage"]) 
        df = df.melt(id_vars=["lineage"])
        df = df.loc[df.notna()["value"]]
        df["study"] = study
        dfs.append(df)
    df = pandas.concat(dfs).rename(columns = {"variable": "sample"})
    df["kingdom"] = [x.split(";")[0] for x in df["lineage"]]
    df["taxon"] = ["".join(reversed([ y for y in "".join(reversed(x)).split(";") if y][0])) for x in df["lineage"]]
    df["genus"] = ["".join(reversed("".join(reversed(x)).split(";")[1])) for x in df["lineage"]]

    stats = []
    stats.append(("Total num studies", total_num_studies))
    stats.append(("Total num samples", total_num_samples))
    stats.append(("Num different taxa", len(set(df["taxon"]))))
    stats.append(("Num different genera", len(set(df["genus"]))))
    stats.append(("Total num data points", len(df)))
    stats.append(("Total num data points for Fungi", len(df.loc[df["kingdom"] == "Fungi"])))
    stats.append(("Total num data points for Blastocystis", len(df.loc[df["genus"] == "Blastocystis"])))
    stats_df = pandas.DataFrame.from_records(stats, columns=['statistic', 'value'])
    stats_df = stats_df.set_index('statistic')

    return df, stats_df


def main(argv=sys.argv[1:]):
    parser = argparse.ArgumentParser(
      description="make diabimmune comparison spreadsheet",
      formatter_class = argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input-results-dir", type=str, action="store", dest="results_dir", required=True)
    parser.add_argument("--input-name-postfix", type=str, action="store", dest="name_postfix", required=True)
    parser.add_argument("--output-xlsx", type=str, action="store", dest="output_xlsx", required=True)
    options=parser.parse_args(argv)
    do(**vars(options))

if __name__ == '__main__':
    main()
