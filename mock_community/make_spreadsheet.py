import pandas
from os import listdir

df = pandas.read_excel('Supplementary Tables S1-9.xlsx',  engine='openpyxl', sheet_name="ST2", skiprows=[0])
df = df[[type(s) == type("") and "D6300" in s for s in df['SampleID'].to_list()]]
d_to_protocol = dict(zip(df['SampleID'], df['Protocol']))

errs = pandas.read_csv("errs.tsv", sep="\t", names=["err", "d"], header=None)
err_to_d = dict(zip(errs["err"], errs["d"]))

def read_counts(f):
   counts = pandas.read_csv(f, sep="\t", names=["err", "c"], header=None)
   return dict(zip(counts["err"], counts["c"]))

err_to_al_all = read_counts("alignment_counts_all.tsv")
err_to_al_best= read_counts("alignment_counts_best.tsv")


read_counts = pandas.read_csv("read_counts.tsv", sep="\t")
err_to_read_count = dict(zip(read_counts["run_accession"], read_counts["read_count"]))

def read_results_folder(sheet_name, err_to_d, d_to_protocol, err_to_read_count, err_to_al, results_d, sep, suffix):
    res_dfs = []
    summaries = []
    for fn in listdir(results_d):
        if fn.endswith(suffix):
            err = fn.split(sep)[0]
            sample_id = err_to_d[err].replace("_", "-")
            protocol =  d_to_protocol[sample_id]
            f = results_d + "/" + fn
            res_df = pandas.read_csv(f, sep="\t")
            res_df["Protocol"] = protocol
            res_df["SampleID"] = sample_id
            res_df["Accession"] = err
            num_added_cols=3
            if 'taxon' in res_df:
                num_cns = res_df[["Cryptococcus_neoformans" in x for x in res_df['taxon'].to_list()]]["taxon_num_reads"].to_list()
                num_scs = res_df[["Saccharomyces_cerevisiae" in x for x in res_df['taxon'].to_list()]]["taxon_num_reads"].to_list()
                num_als = round(sum(res_df["taxon_num_reads"]))
            else:
                # EukDetect output
                num_cns = res_df[["Cryptococcus neoformans" in x for x in res_df['Name'].to_list()]]["Read_counts"].to_list()
                num_scs = res_df[["Saccharomyces cerevisiae" in x for x in res_df['Name'].to_list()]]["Read_counts"].to_list()
                num_als = round(sum(res_df["Read_counts"]))

            # Should be unique, really, but ERR4097232 sees two C. neoformans-like species
            num_cn = sum(num_cns) / len(num_cns)
            num_sc = sum(num_scs) / len(num_scs)
            summaries.append({
                    "Protocol": protocol,
                    "SampleID": sample_id,
                    "Accession": err,
                    "Sheet name": sheet_name,
                    "Number of results": len(res_df),
                    "Total reads": err_to_read_count[err],
                    "Total aligned reads": err_to_al[err],
                    "Total reads passing filters": num_als,
                    "S. cerevisiae reads (shared alignments as fractional)": num_scs[0],
                    "C. neoformans reads (shared alignments as fractional)": num_cns[0],
                    "S. cerevisiae / C. neoformans read ratio": 1.0 * num_scs[0] / num_cns[0]
            })
            res_dfs.append(res_df)
    res = pandas.concat(res_dfs)
    cs = res.columns.to_list()
    res = res[cs[-num_added_cols:] + cs[:-num_added_cols]]
    if "taxon_num_reads" in cs:
        res = res.sort_values(by=["taxon_num_reads"], ascending = False)
    res = res.sort_values(by=["Protocol", "SampleID"])
    summary = pandas.DataFrame.from_dict(summaries)
    return res, summary

legend_columns = ["Sheet name", "All aligments or best alignment used", "Filtering description", "Filtering parameters"]
legend = []
summaries = []
results = {}

def add_results(results_folder, sheet_name, all_or_best, filtering_description, filtering_parameters, sep, ext):
    legend.append([sheet_name, all_or_best, filtering_description, filtering_parameters])
    err_to_al = err_to_al_all if all_or_best == "all" else err_to_al_best if all_or_best == "best" else None
    res, summary = read_results_folder(sheet_name, err_to_d, d_to_protocol, err_to_read_count, err_to_al,  results_folder, sep, ext)
    summaries.append(summary)
    results[sheet_name] = res

def add_our_results(**kwargs):
    return add_results(sep=".", ext = ".tsv", **kwargs)

add_our_results(
        results_folder = "results-CORRAL",
        sheet_name = "CORRAL",
        all_or_best = "all",
        filtering_description = "CORRAL canonical",
        filtering_parameters = "--min-read-query-length 60 --min-taxon-num-markers 2 --min-taxon-num-reads 2 --min-taxon-better-marker-cluster-averages-ratio 1.01 --threshold-avg-match-identity-to-call-known-taxon 0.97  --threshold-num-taxa-to-call-unknown-taxon 1 --threshold-num-markers-to-call-unknown-taxon 4     --threshold-num-reads-to-call-unknown-taxon 8",
        )
add_our_results(
        results_folder = "results-all",
        sheet_name = "All_als",
        all_or_best = "all",
        filtering_description = "All alignments, no filtering",
        filtering_parameters = ""
        )
add_our_results(
        results_folder = "results-all-rl",
        sheet_name = "All_als+RL",
        all_or_best = "all",
        filtering_description = "All alignments, at least 60 bases",
        filtering_parameters = "--min-read-query-length 60"
        )
add_our_results(
        results_folder = "results-all-rl-c",
        sheet_name = "All_als+RL+C",
        all_or_best = "all",
        filtering_description = "All alignments, at least 60 bases, only taxa with enough read + marker counts",
        filtering_parameters = "--min-read-query-length 60 --min-taxon-num-markers 2 --min-taxon-num-reads 2"
        )
add_our_results(
        results_folder = "results-best",
        sheet_name = "Best_al",
        all_or_best = "best",
        filtering_description = "Best alignment only, no filtering",
        filtering_parameters = ""
        )
add_our_results(
        results_folder = "results-best-rl",
        sheet_name = "Best_al+RL",
        all_or_best = "best",
        filtering_description = "Best alignment only, at least 60 bases",
        filtering_parameters = "--min-read-query-length 60"
        )
add_our_results(
        results_folder = "results-best-rl-c",
        sheet_name = "Best_al+RL+C",
        all_or_best = "best",
        filtering_description = "Best alignment only, at least 60 bases, only taxa with enough read + marker counts",
        filtering_parameters = "--min-read-query-length 60 --min-taxon-num-markers 2 --min-taxon-num-alignments 4"
        )
add_our_results(
        results_folder = "results-best-rl-mapq",
        sheet_name = "Best_al+RL+MAPQ",
        all_or_best = "best",
        filtering_description = "Best alignment only, at least 60 bases, MAPQ>=30",
        filtering_parameters = "--min-read-query-length 60 --min-taxon-num-markers 2 --min-read-mapq 30"
        )
add_our_results(
        results_folder = "results-best-rl-mapq-c",
        sheet_name = "Best_al+RL+MAPQ+C",
        all_or_best = "best",
        filtering_description = "Best alignment only, at least 60 bases, MAPQ>=30, only taxa with enough read + marker counts",
        filtering_parameters = "--min-read-query-length 60 --min-taxon-num-markers 2 --min-read-mapq 30 --min-taxon-num-markers 2 --min-taxon-num-alignments 4"
        )
add_results(
        results_folder = "results-EukDetect",
        sheet_name = "EukDetect",
        all_or_best = "best",
        filtering_description = "EukDetect final results",
        filtering_parameters = "_filtered_hits.table.txt output",
        sep = "_",
        ext = "_filtered_hits_table.txt",
        )
add_results(
        results_folder = "results-EukDetect/filtering",
        sheet_name = "EukDetect all hits",
        all_or_best = "best",
        filtering_description = "EukDetect intermediate results",
        filtering_parameters = "_all_hits.table.txt output",
        sep = "_",
        ext = "_all_hits_table.txt",
        )

summary = pandas.concat(summaries).sort_values(by=["Protocol", "SampleID", "Sheet name"])
#results_rl_only = read_results_folder(err_to_d, d_to_protocol, err_to_read_count, err_to_al_best,  "results-rl-only", ".", ".tsv")
#results_c_only = read_results_folder(err_to_d, d_to_protocol, err_to_read_count, err_to_al_best,  "results-c-only", ".", ".tsv")
#results_mapq_only = read_results_folder(err_to_d, d_to_protocol, err_to_read_count, err_to_al_best,  "results-mapq-only", ".", ".tsv")
#results_eukdetect = read_results_folder(err_to_d, d_to_protocol, err_to_read_count, err_to_al_best,  "results-EukDetect", "_", "_filtered_hits_table.txt")

def save_df(writer, sheet_name, df, float_format="%.4f", width = None):
    df.to_excel(writer, sheet_name = sheet_name, float_format=float_format, index=False)
    worksheet = writer.sheets[sheet_name]  # pull worksheet object
    lengths = [len(str(df[col].name)) for idx, col in enumerate(df)]
    max_length = width or max(lengths)
    worksheet.set_column(0, len(lengths), max_length)

#legend.append()

# TODO legend
#legend.append(["best+RL+C+MAPQ", "best", "read length, counts per taxon, MAPQ", "--min-read-query-length 60 --min-taxon-num-markers 2 --min-taxon-num-alignments 4 --min-read-mapq 30"])

# TODO EukDetect all hits
# results-EukDetect/filtering/ERR4097268_all_hits_table.txt

with pandas.ExcelWriter("all-results.xlsx", engine="xlsxwriter") as writer:
    save_df(writer, "Legend", pandas.DataFrame(columns = legend_columns, data = legend))
    save_df(writer, "Aggregated values", summary, float_format="%.2f", width = 30)
    for sheet_name in results:
        save_df(writer, sheet_name, results[sheet_name])
