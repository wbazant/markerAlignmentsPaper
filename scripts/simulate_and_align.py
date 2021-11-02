import sys
import re
import os
import argparse
import pysam
import re
import json
import statistics
import logging
import subprocess
from pathlib import Path


from marker_alignments.pysam2 import compute_alignment_identity

def read_marker_to_taxon(path):
    result = {}
    with open(path, 'r') as f:
        for line in f:
            (marker, taxon) = line.rstrip().split("\t")
            result[marker] = taxon
    return result

def main_loop(**kwargs):
    result = []
   
    marker_to_taxon = read_marker_to_taxon(kwargs['refdb_marker_to_taxon_path'])

    for mutation_rate in kwargs['mutation_rates']:
        for base_error_rate in kwargs['base_error_rates']:
            for read_length in kwargs['read_lengths']:
                d = do_one(marker_to_taxon, read_length, base_error_rate, mutation_rate, **kwargs)
                d['mutation_rate'] = mutation_rate
                d['base_error_rate'] = base_error_rate
                d['read_length'] = read_length
                result.append(d)
    print(json.dumps(result, indent = 4))

def do_one(marker_to_taxon, read_length, base_error_rate, mutation_rate, wd, logger, ref_db, sim_source, seed, **kwargs):
    Path(wd).mkdir(parents=True, exist_ok=True)

    prefix = f"{read_length}.{base_error_rate}.{mutation_rate}"
    sim_path_1 = f"{wd}/{prefix}.wgsim_1.fq"
    sim_path_2 = f"{wd}/{prefix}.wgsim_2.fq"
    if not (os.path.isfile(sim_path_1) and os.path.isfile(sim_path_2)):
        wgsim_cmd = ["wgsim",
            "-S", str(seed),
            "-1", str(read_length), "-2", str(read_length),
            "-e", str(base_error_rate),
            "-r", str(mutation_rate),
            sim_source,
            sim_path_1,
            sim_path_2,
        ]
        logger.info("Running: " + " ".join(wgsim_cmd))
        subprocess.run(wgsim_cmd, check=True, stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL)

    sam_path= f"{wd}/{prefix}.alignments.sam"
    if not os.path.isfile(sam_path):
        sam_cmd = ["bowtie2",
            "-x", ref_db,
            "-1", sim_path_1,
            "-2", sim_path_2,
            "--omit-sec-seq",
            "--no-discordant",
            "--no-unal",
            "--seed", str(seed),
        ]
        logger.info("Running: " + " ".join(sam_cmd))
        with open(sam_path + ".tmp", 'w') as f:
            subprocess.run(sam_cmd, check=True, stdout = f, stderr = sys.stderr)
        subprocess.run(["mv", "-v", sam_path + ".tmp", sam_path])

    summary_path= f"{wd}/{prefix}.summary.json"

    if not os.path.isfile(summary_path):
        logger.info("Summarizing " + " ".join([sam_path, sim_path_1]))
        summary = summarize_sim_alignment(marker_to_taxon, sam_path)
        grep_cmd = ['grep', '-c', '^@', sim_path_1]
        logger.info("Running: " + " ".join(grep_cmd))
        summary["numQueriesTotal"] = int(subprocess.run(grep_cmd, stdout = subprocess.PIPE).stdout.decode(encoding="utf-8").replace("\n", ""))
        grep_cmd_2 = [
            'grep', '-o', '^@.*at2759', sim_path_1, '|', "uniq", '|', "wc", "-l" ]
        logger.info("Running: " + " ".join(grep_cmd_2))
        summary["numBuscosTotal"] = int(subprocess.run([" ".join(grep_cmd_2)], shell = True, stdout = subprocess.PIPE).stdout.decode(encoding="utf-8").replace("\n", ""))
        summary["precisionQueries"] = 1.0 * summary["numQueriesMappedOnlyAsMatch"] / summary["numQueriesMapped"]
        summary["precisionBuscos"] = 1.0 * summary["numBuscosMappedOnlyAsMatch"] / summary["numBuscosMapped"]
        summary["recallQueries"] = 1.0 * summary["numQueriesMappedOnlyAsMatch"] / summary["numQueriesTotal"]
        summary["recallBuscos"] = 1.0 * summary["numBuscosMappedOnlyAsMatch"] / summary["numBuscosTotal"]


        with open(summary_path, 'w') as f:
            json.dump(summary, f)

    
    with open(summary_path, 'r') as f:
        return json.load(f)

# removed $ from the end compared to marker_alignments version
#other_MGCollapse_EPSP3-12_382_886_0:0:0_0:0:0_0
eukprot_refdb_regex_taxon = "^[a-z]+-(.*_[a-z-]+(?:_.*)?)-[0-9]+at2759-[A-Z]\d|^[a-z]+-(.*)...Collapse_[^_]*|^()other_MGCollapse_[^_]*"
eukprot_refdb_regex_marker = "^[a-z]+-.*_[a-z-]+(?:_.*)?-([0-9]+at2759)-[A-Z]\d|^[a-z]+-.*(..Collapse_[^_]*)|^(other_MGCollapse_[^_]*)"
eukprot_refdb_regex_busco = "^([a-z]+-.*_[a-z-]+(?:_.*)?-([0-9]+at2759)-[A-Z]\d|^[a-z]+-.*..Collapse_[^_]*|^other_MGCollapse_[^_]*)"

pattern_taxon = re.compile(eukprot_refdb_regex_taxon)
pattern_marker = re.compile(eukprot_refdb_regex_marker)
pattern_busco = re.compile(eukprot_refdb_regex_busco)

def next_g(search):
    return next(g for g in search.groups() if g is not None)

def match(marker_to_taxon, query_name, reference_name):
    # not a cross-match if the read aligns to the other copy of the bu(m)sco

    source_taxon_search = pattern_taxon.search(query_name) 
    if not source_taxon_search:
        raise ValueError(query_name)
    source_taxon = next_g(source_taxon_search)

    source_marker_search = pattern_marker.search(query_name) 
    source_marker = next_g(source_marker_search)

    source_busco_search = pattern_busco.search(query_name) 
    source_busco = next_g(source_busco_search)

    matched_taxon_search = pattern_taxon.search(reference_name) 
    if not matched_taxon_search:
        raise ValueError(reference_name)
    matched_taxon = next_g(matched_taxon_search)

    matched_marker_search = pattern_marker.search(reference_name) 
    matched_marker = next_g(matched_marker_search)

    matched_busco_search = pattern_busco.search(reference_name) 
    matched_busco = next_g(matched_busco_search)



    is_match = re.sub("-D\d$", "", source_busco) == re.sub("-D\d$", "", matched_busco)
    is_match_2 = (source_taxon == matched_taxon and source_marker == matched_marker)

    if is_match != is_match_2:
        raise ValueError(query_name, reference_name, is_match, is_match_2, source_taxon, matched_taxon, source_marker, matched_marker)

    source_taxid = marker_to_taxon[source_busco]
    matched_taxid = marker_to_taxon[matched_busco]
    
    return source_taxid, source_busco, matched_taxid, matched_busco, is_match

def summarize_sim_alignment(marker_to_taxon, input_alignment_file):
    alignment_file = pysam.AlignmentFile(input_alignment_file, check_sq=False)

    queries = {}
    matched_buscos = {}
    mapqs = {}
    query_lengths = {}

    identities_match_true = []
    identities_match_false = []

    for read in alignment_file.fetch():
        query = read.query_name
        identity = compute_alignment_identity(read)
        mapq = read.mapq
        query_length = read.infer_query_length()

        source_taxon, source_busco, matched_taxon, matched_busco, is_match = match(marker_to_taxon, query, read.reference_name)

        if matched_busco not in matched_buscos:
            matched_buscos[matched_busco] = {True: 0, False: 0}
        matched_buscos[matched_busco][is_match] += 1

        if query not in queries:
            queries[query] = {True: 0, False: 0}
        queries[query][is_match] += 1

        if mapq not in mapqs:
            mapqs[mapq] = 0
        mapqs[mapq] += 1

        if query_length not in query_lengths:
            query_lengths[query_length] = 0
        query_lengths[query_length] += 1

        if is_match:
            identities_match_true.append(identity)
        else:
            identities_match_false.append(identity)


    mapqs_n = 0.0
    mapqs_n_at_least_30 = 0.0
    mapqs_sum = 0.0

    if 42 not in mapqs:
        mapqs[42] = 0

    for k in mapqs:
        mapqs_n += mapqs[k]
        if k >= 30:
            mapqs_n_at_least_30 += mapqs[k]
        mapqs_sum += k * mapqs[k]

    mapqs_avg = mapqs_sum / mapqs_n if mapqs_n else 0

    query_lengths_n = 0.0
    query_lengths_n_at_least_60 = 0.0
    query_lengths_sum = 0.0

    for k in query_lengths:
        query_lengths_n += query_lengths[k]
        if k >= 60:
            query_lengths_n_at_least_60 += query_lengths[k]
        query_lengths_sum += k * query_lengths[k]

    query_lengths_avg = query_lengths_sum / query_lengths_n if query_lengths_n else 0
    return {
        "numQueriesMapped": len (queries.keys()),
        "numQueriesMappedOnlyAsMatch": len ([k for k in queries if queries[k][True] > 0 and queries[k][False] == 0 ]),
        "numBuscosMapped": len (matched_buscos.keys()),
        "numBuscosMappedOnlyAsMatch": len ([k for k in matched_buscos if matched_buscos[k][True] > 0 and matched_buscos[k][False] == 0 ]),
        "meanIdentitiesMatchTrue": statistics.mean(identities_match_true) if identities_match_true else 0.0,
        "meanIdentitiesMatchFalse": statistics.mean(identities_match_false) if identities_match_false else 0.0,
        "mapqsAvg": mapqs_avg,
        "mapqsFractionAtLeast30": mapqs_n_at_least_30 / mapqs_n if mapqs_n else 0,
        "mapqsFraction42": mapqs[42] / mapqs_n if mapqs_n else 0,
        "query_lengthsAvg": query_lengths_avg,
        "query_lengthsFractionAtLeast60": query_lengths_n_at_least_60 / query_lengths_n if query_lengths_n else 0,
    }




def main(argv=sys.argv[1:]):
    parser = argparse.ArgumentParser(
      description="simulate_and_align",
      formatter_class = argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--reference", type=str, action="store", dest="ref_db", required=True)
    parser.add_argument("--sim-source", type=str, action="store", dest="sim_source", required=True)
    parser.add_argument("--dir", type=str, action="store", dest="wd", default = os.getcwd())

    parser.add_argument("--verbose", action="store_true", dest="verbose")
    parser.add_argument("--read-lengths", nargs='+', type=int, action = "store", dest = "read_lengths", required = True)
    parser.add_argument("--base-error-rates", nargs='+', type=float,action = "store",  dest = "base_error_rates", required = True)
    parser.add_argument("--mutation-rates", nargs='+', type=float,action = "store",  dest = "mutation_rates", required = True)
    parser.add_argument("--seed", type=int, default = 1337)
    parser.add_argument("--refdb-marker-to-taxon-path", type=str, action="store", dest="refdb_marker_to_taxon_path", help = "Lookup file, two columns - marker name, taxon name")

    options=parser.parse_args(argv)


    logging.basicConfig(format='%(asctime)s: %(message)s')
    logger=logging.getLogger(__name__)
    if options.verbose:
        logger.setLevel(logging.DEBUG)

    main_loop(**vars(options), logger = logger )
    


if __name__ == '__main__':
    main()
