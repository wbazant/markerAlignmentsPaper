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
import sqlite3

from marker_alignments.store import SqliteStore
from ete3 import NCBITaxa

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

    ncbi = NCBITaxa(kwargs['refdb_ncbi'])
    for mutation_rate in kwargs['mutation_rates']:
        for base_error_rate in kwargs['base_error_rates']:
            for read_length in kwargs['read_lengths']:
                d = do_one(ncbi, marker_to_taxon, read_length, base_error_rate, mutation_rate, **kwargs)
                d['mutation_rate'] = mutation_rate
                d['base_error_rate'] = base_error_rate
                d['read_length'] = read_length
                result.append(d)
    print(json.dumps(result, indent = 4))

def do_one(ncbi, marker_to_taxon, read_length, base_error_rate, mutation_rate, wd, logger, ref_db, sim_source, seed, **kwargs):
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
        subprocess.run(["mv", sam_path + ".tmp", sam_path])

    db_path= f"{wd}/{prefix}.alignments.sqlite"

    if not os.path.isfile(db_path):
        logger.info("Populating store " + " ".join([sam_path, db_path]))
        store = write_sim_alignment_to_store(ncbi, marker_to_taxon, sam_path, db_path)
    else:
        store = SqliteStore(db_path = db_path)
        store.connect()
    
    summary_path= f"{wd}/{prefix}.summary.json"
    if not os.path.isfile(summary_path):
        logger.info("Summarizing " + " ".join([db_path, sim_path_1]))
        summary = stats_from_store(store)
        grep_cmd = ['grep', '-c', '^@', sim_path_1]
        logger.info("Running: " + " ".join(grep_cmd))
        summary["numQueriesTotal"] = int(subprocess.run(grep_cmd, stdout = subprocess.PIPE).stdout.decode(encoding="utf-8").replace("\n", ""))

        summary["numBuscosTotal"] = count_buscos_in_sim_fastq(sim_path_1)
        summary["precisionQueries"] = 1.0 * summary["numQueriesMappedOnlyAsMatch"] / summary["numQueriesMapped"]
        summary["precisionBuscos"] = 1.0 * summary["numBuscosMappedOnlyAsMatch"] / summary["numBuscosMapped"]
        summary["recallQueries"] = 1.0 * summary["numQueriesMappedOnlyAsMatch"] / summary["numQueriesTotal"]
        summary["recallBuscos"] = 1.0 * summary["numBuscosMappedOnlyAsMatch"] / summary["numBuscosTotal"]
        summary["precisionSameSpeciesQueries"] = 1.0 * summary["numQueriesMappedOnlyAsSameSpecies"] / summary["numQueriesMapped"]
#        summary["precisionSameSpeciesBuscos"] = 1.0 * summary["numBuscosMappedOnlyAsSameSpecies"] / summary["numBuscosMapped"]
        summary["recallSameSpeciesQueries"] = 1.0 * summary["numQueriesMappedOnlyAsSameSpecies"] / summary["numQueriesTotal"]
#        summary["recallSameSpeciesBuscos"] = 1.0 * summary["numBuscosMappedOnlyAsSameSpecies"] / summary["numBuscosTotal"]
        summary["precisionSameGenusQueries"] = 1.0 * summary["numQueriesMappedOnlyAsSameGenus"] / summary["numQueriesMapped"]
#        summary["precisionSameGenusBuscos"] = 1.0 * summary["numBuscosMappedOnlyAsSameGenus"] / summary["numBuscosMapped"]
        summary["recallSameGenusQueries"] = 1.0 * summary["numQueriesMappedOnlyAsSameGenus"] / summary["numQueriesTotal"]
#        summary["recallSameGenusBuscos"] = 1.0 * summary["numBuscosMappedOnlyAsSameGenus"] / summary["numBuscosTotal"]
        summary["precisionSameFamilyQueries"] = 1.0 * summary["numQueriesMappedOnlyAsSameFamily"] / summary["numQueriesMapped"]
#        summary["precisionSameFamilyBuscos"] = 1.0 * summary["numBuscosMappedOnlyAsSameFamily"] / summary["numBuscosMapped"]
        summary["recallSameFamilyQueries"] = 1.0 * summary["numQueriesMappedOnlyAsSameFamily"] / summary["numQueriesTotal"]
#        summary["recallSameFamilyBuscos"] = 1.0 * summary["numBuscosMappedOnlyAsSameFamily"] / summary["numBuscosTotal"]
        summary["precisionWhenMapqAtLeast30Queries"] = 1.0 * summary["numQueriesMappedOnlyAsMatchWhenMapqAtLeast30"] / summary["numQueriesMappedWhenMapqAtLeast30"]
#        summary["precisionWhenMapqAtLeast30Buscos"] = 1.0 * summary["numBuscosMappedOnlyAsMatchWhenMapqAtLeast30"] / summary["numBuscosMappedWhenMapqAtLeast30"]
        summary["recallWhenMapqAtLeast30Queries"] = 1.0 * summary["numQueriesMappedOnlyAsMatchWhenMapqAtLeast30"] / summary["numQueriesTotal"]
#        summary["recallWhenMapqAtLeast30Buscos"] = 1.0 * summary["numBuscosMappedOnlyAsMatchWhenMapqAtLeast30"] / summary["numBuscosTotal"]
        summary["precisionSameSpeciesWhenMapqAtLeast30Queries"] = 1.0 * summary["numQueriesMappedOnlyAsSameSpeciesWhenMapqAtLeast30"] / summary["numQueriesMappedWhenMapqAtLeast30"]
#        summary["precisionSameSpeciesWhenMapqAtLeast30Buscos"] = 1.0 * summary["numBuscosMappedOnlyAsSameSpeciesWhenMapqAtLeast30"] / summary["numBuscosMappedWhenMapqAtLeast30"]
        summary["recallSameSpeciesWhenMapqAtLeast30Queries"] = 1.0 * summary["numQueriesMappedOnlyAsSameSpeciesWhenMapqAtLeast30"] / summary["numQueriesTotal"]
#        summary["recallSameSpeciesWhenMapqAtLeast30Buscos"] = 1.0 * summary["numBuscosMappedOnlyAsSameSpeciesWhenMapqAtLeast30"] / summary["numBuscosTotal"]
        summary["precisionSameGenusWhenMapqAtLeast30Queries"] = 1.0 * summary["numQueriesMappedOnlyAsSameGenusWhenMapqAtLeast30"] / summary["numQueriesMappedWhenMapqAtLeast30"]
#        summary["precisionSameGenusWhenMapqAtLeast30Buscos"] = 1.0 * summary["numBuscosMappedOnlyAsSameGenusWhenMapqAtLeast30"] / summary["numBuscosMappedWhenMapqAtLeast30"]
        summary["recallSameGenusWhenMapqAtLeast30Queries"] = 1.0 * summary["numQueriesMappedOnlyAsSameGenusWhenMapqAtLeast30"] / summary["numQueriesTotal"]
#        summary["recallSameGenusWhenMapqAtLeast30Buscos"] = 1.0 * summary["numBuscosMappedOnlyAsSameGenusWhenMapqAtLeast30"] / summary["numBuscosTotal"]
        summary["precisionSameFamilyWhenMapqAtLeast30Queries"] = 1.0 * summary["numQueriesMappedOnlyAsSameFamilyWhenMapqAtLeast30"] / summary["numQueriesMappedWhenMapqAtLeast30"]
#        summary["precisionSameFamilyWhenMapqAtLeast30Buscos"] = 1.0 * summary["numBuscosMappedOnlyAsSameFamilyWhenMapqAtLeast30"] / summary["numBuscosMappedWhenMapqAtLeast30"]
        summary["recallSameFamilyWhenMapqAtLeast30Queries"] = 1.0 * summary["numQueriesMappedOnlyAsSameFamilyWhenMapqAtLeast30"] / summary["numQueriesTotal"]
#        summary["recallSameFamilyWhenMapqAtLeast30Buscos"] = 1.0 * summary["numBuscosMappedOnlyAsSameFamilyWhenMapqAtLeast30"] / summary["numBuscosTotal"]
        summary["precisionWhenMapqAtLeast5Queries"] = 1.0 * summary["numQueriesMappedOnlyAsMatchWhenMapqAtLeast5"] / summary["numQueriesMappedWhenMapqAtLeast5"]
#        summary["precisionWhenMapqAtLeast5Buscos"] = 1.0 * summary["numBuscosMappedOnlyAsMatchWhenMapqAtLeast5"] / summary["numBuscosMappedWhenMapqAtLeast5"]
        summary["recallWhenMapqAtLeast5Queries"] = 1.0 * summary["numQueriesMappedOnlyAsMatchWhenMapqAtLeast5"] / summary["numQueriesTotal"]
#        summary["recallWhenMapqAtLeast5Buscos"] = 1.0 * summary["numBuscosMappedOnlyAsMatchWhenMapqAtLeast5"] / summary["numBuscosTotal"]
        summary["precisionSameSpeciesWhenMapqAtLeast5Queries"] = 1.0 * summary["numQueriesMappedOnlyAsSameSpeciesWhenMapqAtLeast5"] / summary["numQueriesMappedWhenMapqAtLeast5"]
#        summary["precisionSameSpeciesWhenMapqAtLeast5Buscos"] = 1.0 * summary["numBuscosMappedOnlyAsSameSpeciesWhenMapqAtLeast5"] / summary["numBuscosMappedWhenMapqAtLeast5"]
        summary["recallSameSpeciesWhenMapqAtLeast5Queries"] = 1.0 * summary["numQueriesMappedOnlyAsSameSpeciesWhenMapqAtLeast5"] / summary["numQueriesTotal"]
#        summary["recallSameSpeciesWhenMapqAtLeast5Buscos"] = 1.0 * summary["numBuscosMappedOnlyAsSameSpeciesWhenMapqAtLeast5"] / summary["numBuscosTotal"]
        summary["precisionSameGenusWhenMapqAtLeast5Queries"] = 1.0 * summary["numQueriesMappedOnlyAsSameGenusWhenMapqAtLeast5"] / summary["numQueriesMappedWhenMapqAtLeast5"]
#        summary["precisionSameGenusWhenMapqAtLeast5Buscos"] = 1.0 * summary["numBuscosMappedOnlyAsSameGenusWhenMapqAtLeast5"] / summary["numBuscosMappedWhenMapqAtLeast5"]
        summary["recallSameGenusWhenMapqAtLeast5Queries"] = 1.0 * summary["numQueriesMappedOnlyAsSameGenusWhenMapqAtLeast5"] / summary["numQueriesTotal"]
#        summary["recallSameGenusWhenMapqAtLeast5Buscos"] = 1.0 * summary["numBuscosMappedOnlyAsSameGenusWhenMapqAtLeast5"] / summary["numBuscosTotal"]
        summary["precisionSameFamilyWhenMapqAtLeast5Queries"] = 1.0 * summary["numQueriesMappedOnlyAsSameFamilyWhenMapqAtLeast5"] / summary["numQueriesMappedWhenMapqAtLeast5"]
#        summary["precisionSameFamilyWhenMapqAtLeast5Buscos"] = 1.0 * summary["numBuscosMappedOnlyAsSameFamilyWhenMapqAtLeast5"] / summary["numBuscosMappedWhenMapqAtLeast5"]
        summary["recallSameFamilyWhenMapqAtLeast5Queries"] = 1.0 * summary["numQueriesMappedOnlyAsSameFamilyWhenMapqAtLeast5"] / summary["numQueriesTotal"]
#        summary["recallSameFamilyWhenMapqAtLeast5Buscos"] = 1.0 * summary["numBuscosMappedOnlyAsSameFamilyWhenMapqAtLeast5"] / summary["numBuscosTotal"]


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

def count_buscos_in_sim_fastq(sim_path):
    result = set()
    with open(sim_path, 'r') as f:
        while True:
            l = f.readline()
            if not l:
                break
            if l.startswith("@"):
                l = re.sub("^@", "", l)
                search = pattern_busco.search(l)
                if not search:
                    raise ValueError(l)
                busco = next_g(search)
                result.add(busco)
    return len(result)


trad_ranks = {"superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"}

match_type_cache = {}
def get_match_type(ncbi, source_taxid, matched_taxid):
    k = (source_taxid, matched_taxid)
    if k in match_type_cache:
        return match_type_cache[k]
    tree = ncbi.get_topology([source_taxid, matched_taxid])

    r = tree.get_tree_root()
    if hasattr(r, 'rank') and r.rank in trad_ranks:
        result = r.rank
    else:
        lineage = ncbi.get_lineage(r.name)
        ranks = ncbi.get_rank(lineage)
        ranks_increasing = [ranks[x] for x in reversed(lineage) if ranks[x] in trad_ranks]

        if not ranks_increasing:
            raise ValueError(lineage, ranks)
        result = ranks_increasing[0]

    match_type_cache[k] = result
    return result

def match(ncbi, marker_to_taxon, query_name, reference_name):
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

    if is_match:
        match_type = 'true_match'
    elif source_taxid == matched_taxid:
        match_type = 'species'
    else:
        match_type = get_match_type(ncbi, source_taxid, matched_taxid)

    is_corresponding_markers = source_marker == matched_marker

    if match_type == 'species' and is_corresponding_markers:
        raise ValueError(source_busco, matched_busco, source_marker)


    return source_taxid, source_busco, matched_taxid, matched_busco, match_type

def query_one_number(store, sql):
    result = [ k for k in store.query(sql)]
    if len(result) != 1:
        raise ValueError(result, sql)
    return result[0][0]

def query_ratio(store, sql):
    result = [ k for k in store.query(sql)]
    if len(result) != 1:
        raise ValueError(result, sql)
    return 1.0 * result[0][0] / result[0][1]

def write_sim_alignment_to_store(ncbi, marker_to_taxon, input_alignment_file, db_path):
    alignment_file = pysam.AlignmentFile(input_alignment_file, check_sq=False)

    store = SqliteStore(db_path = db_path)
    store.connect()
    store.do('''create table alignment_from_known_source (
              query text not null,
              source_taxon text not null,
              source_busco text not null,
              matched_taxon text not null,
              matched_busco text not null,
              match_type text not null,
              identity real not null,
              mapq real not null,
              query_length real not null
            );''')
    store.start_bulk_write()

    for read in alignment_file.fetch():
        query = read.query_name
        identity = compute_alignment_identity(read)
        mapq = read.mapq
        query_length = read.infer_query_length()

        source_taxon, source_busco, matched_taxon, matched_busco, match_type = match(ncbi, marker_to_taxon, query, read.reference_name)

        store.do('insert into alignment_from_known_source (query,source_taxon,source_busco,matched_taxon,matched_busco,match_type,identity,mapq,query_length) values(?,?,?,?,?,?,?,?,?)', [query,source_taxon,source_busco,matched_taxon,matched_busco,match_type,identity,mapq,query_length])

    store.end_bulk_write()
    return store

def stats_from_store(store):
    return {
        "numQueriesMapped": query_one_number(store, '''
          select count(distinct query)
            from alignment_from_known_source
        '''),
        "numQueriesMappedOnlyAsMatch": query_one_number(store, '''
         select count(*) from (
              select a.query, sum(c) as allowed_matches, count(c) as all_matches
              from (
                select query, case when match_type in ('true_match') then 1 else 0 end as c
                from alignment_from_known_source
              ) a
              group by a.query
              having allowed_matches = all_matches
        )
        '''),
        "numQueriesMappedOnlyAsSameSpecies": query_one_number(store, '''
         select count(*) from (
              select a.query, sum(c) as allowed_matches, count(c) as all_matches
              from (
                select query, case when match_type in ('true_match', 'species') then 1 else 0 end as c
                from alignment_from_known_source
              ) a
              group by a.query
              having allowed_matches = all_matches
        )
        '''),
        "numQueriesMappedOnlyAsSameGenus": query_one_number(store, '''
         select count(*) from (
              select a.query, sum(c) as allowed_matches, count(c) as all_matches
              from (
                select query, case when match_type in ('true_match', 'species', 'genus') then 1 else 0 end as c
                from alignment_from_known_source
              ) a
              group by a.query
              having allowed_matches = all_matches
        )
        '''),
        "numQueriesMappedOnlyAsSameFamily": query_one_number(store, '''
         select count(*) from (
              select a.query, sum(c) as allowed_matches, count(c) as all_matches
              from (
                select query, case when match_type in ('true_match', 'species', 'genus', 'family') then 1 else 0 end as c
                from alignment_from_known_source
              ) a
              group by a.query
              having allowed_matches = all_matches
        )
        '''),
        "numQueriesMappedWhenMapqAtLeast30": query_one_number(store, '''
          select count(distinct query)
            from alignment_from_known_source where mapq >=30
        '''),
        "numQueriesMappedOnlyAsMatchWhenMapqAtLeast30": query_one_number(store, '''
         select count(*) from (
              select a.query, sum(c) as allowed_matches, count(c) as all_matches
              from (
                select query, case when match_type in ('true_match') then 1 else 0 end as c
                from alignment_from_known_source where mapq >=30
              ) a
              group by a.query
              having allowed_matches = all_matches
        )
        '''),
        "numQueriesMappedOnlyAsSameSpeciesWhenMapqAtLeast30": query_one_number(store, '''
         select count(*) from (
              select a.query, sum(c) as allowed_matches, count(c) as all_matches
              from (
                select query, case when match_type in ('true_match', 'species') then 1 else 0 end as c
                from alignment_from_known_source where mapq >=30
              ) a
              group by a.query
              having allowed_matches = all_matches
        )
        '''),
        "numQueriesMappedOnlyAsSameGenusWhenMapqAtLeast30": query_one_number(store, '''
         select count(*) from (
              select a.query, sum(c) as allowed_matches, count(c) as all_matches
              from (
                select query, case when match_type in ('true_match', 'species', 'genus') then 1 else 0 end as c
                from alignment_from_known_source where mapq >=30
              ) a
              group by a.query
              having allowed_matches = all_matches
        )
        '''),
        "numQueriesMappedOnlyAsSameFamilyWhenMapqAtLeast30": query_one_number(store, '''
         select count(*) from (
              select a.query, sum(c) as allowed_matches, count(c) as all_matches
              from (
                select query, case when match_type in ('true_match', 'species', 'genus', 'family') then 1 else 0 end as c
                from alignment_from_known_source where mapq >=30
              ) a
              group by a.query
              having allowed_matches = all_matches
        )
        '''),
        "numQueriesMappedWhenMapqAtLeast5": query_one_number(store, '''
          select count(distinct query)
            from alignment_from_known_source where mapq >=5
        '''),
        "numQueriesMappedOnlyAsMatchWhenMapqAtLeast5": query_one_number(store, '''
         select count(*) from (
              select a.query, sum(c) as allowed_matches, count(c) as all_matches
              from (
                select query, case when match_type in ('true_match') then 1 else 0 end as c
                from alignment_from_known_source where mapq >=5
              ) a
              group by a.query
              having allowed_matches = all_matches
        )
        '''),
        "numQueriesMappedOnlyAsSameSpeciesWhenMapqAtLeast5": query_one_number(store, '''
         select count(*) from (
              select a.query, sum(c) as allowed_matches, count(c) as all_matches
              from (
                select query, case when match_type in ('true_match', 'species') then 1 else 0 end as c
                from alignment_from_known_source where mapq >=5
              ) a
              group by a.query
              having allowed_matches = all_matches
        )
        '''),
        "numQueriesMappedOnlyAsSameGenusWhenMapqAtLeast5": query_one_number(store, '''
         select count(*) from (
              select a.query, sum(c) as allowed_matches, count(c) as all_matches
              from (
                select query, case when match_type in ('true_match', 'species', 'genus') then 1 else 0 end as c
                from alignment_from_known_source where mapq >=5
              ) a
              group by a.query
              having allowed_matches = all_matches
        )
        '''),
        "numQueriesMappedOnlyAsSameFamilyWhenMapqAtLeast5": query_one_number(store, '''
         select count(*) from (
              select a.query, sum(c) as allowed_matches, count(c) as all_matches
              from (
                select query, case when match_type in ('true_match', 'species', 'genus', 'family') then 1 else 0 end as c
                from alignment_from_known_source where mapq >=5
              ) a
              group by a.query
              having allowed_matches = all_matches
        )
        '''),
        "numBuscosMapped": query_one_number(store, '''
          select count(distinct source_busco)
            from alignment_from_known_source
        '''),
        "numBuscosMappedOnlyAsMatch": query_one_number(store, '''
         select count(distinct source_busco) from (
              select a.source_busco, a.query, sum(c) as allowed_matches, count(c) as all_matches
              from (
                select source_busco, query, case when match_type in ('true_match') then 1 else 0 end as c
                from alignment_from_known_source
              ) a
              group by a.source_busco, a.query
              having allowed_matches = all_matches
        )
        '''),
        "numBuscosMappedOnlyAsSameSpecies": query_one_number(store, '''
         select count(distinct source_busco) from (
              select a.source_busco, a.query, sum(c) as allowed_matches, count(c) as all_matches
              from (
                select source_busco, query, case when match_type in ('true_match', 'species') then 1 else 0 end as c
                from alignment_from_known_source
              ) a
              group by a.source_busco, a.query
              having allowed_matches = all_matches
        )
        '''),
        "numBuscosMappedOnlyAsSameGenus": query_one_number(store, '''
         select count(distinct source_busco) from (
              select a.source_busco, a.query, sum(c) as allowed_matches, count(c) as all_matches
              from (
                select source_busco, query, case when match_type in ('true_match', 'species', 'genus') then 1 else 0 end as c
                from alignment_from_known_source
              ) a
              group by a.source_busco, a.query
              having allowed_matches = all_matches
        )
        '''),
        "numBuscosMappedOnlyAsSameFamily": query_one_number(store, '''
         select count(distinct source_busco) from (
              select a.source_busco, a.query, sum(c) as allowed_matches, count(c) as all_matches
              from (
                select source_busco, query, case when match_type in ('true_match', 'species', 'genus', 'family') then 1 else 0 end as c
                from alignment_from_known_source
              ) a
              group by a.source_busco, a.query
              having allowed_matches = all_matches
        )
        '''),
        "meanIdentitiesMappedOnlyAsMatch": query_one_number(store, '''
         select avg(identity) from (
              select a.query, identity, sum(c) as allowed_matches, count(c) as all_matches
              from (
                select query, identity, case when match_type in ('true_match') then 1 else 0 end as c
                from alignment_from_known_source
              ) a
              group by a.query
              having allowed_matches = all_matches
        )
        '''),
        "meanIdentitiesMappedOnlyAsSameSpecies": query_one_number(store, '''
         select avg(identity) from (
              select a.query, identity, sum(c) as allowed_matches, count(c) as all_matches
              from (
                select query, identity, case when match_type in ('true_match', 'species') then 1 else 0 end as c
                from alignment_from_known_source
              ) a
              group by a.query
              having allowed_matches = all_matches
        )
        '''),
        "meanIdentitiesMappedOnlyAsSameGenus": query_one_number(store, '''
         select avg(identity) from (
              select a.query, identity, sum(c) as allowed_matches, count(c) as all_matches
              from (
                select query, identity, case when match_type in ('true_match', 'species', 'genus') then 1 else 0 end as c
                from alignment_from_known_source
              ) a
              group by a.query
              having allowed_matches = all_matches
        )
        '''),
        "meanIdentitiesMappedOnlyAsSameFamily": query_one_number(store, '''
         select avg(identity) from (
              select a.query, identity, sum(c) as allowed_matches, count(c) as all_matches
              from (
                select query, identity, case when match_type in ('true_match', 'species', 'genus', 'family') then 1 else 0 end as c
                from alignment_from_known_source
              ) a
              group by a.query
              having allowed_matches = all_matches
        )
        '''),
        "meanTopIdentityPerQuery": query_one_number(store, '''
         select avg(top_identity) from (
              select a.query, max(identity) as top_identity
              from alignment_from_known_source a
              group by a.query
        )
        '''),
        "mapqsAvg": query_one_number(store, '''
         select avg(mapq) from alignment_from_known_source a
        '''),
        "mapqsFractionAtLeast30": query_ratio(store, '''
          select sum(c) as numerator, count(c) as deliminator
            from (
                select case when mapq > 30 then 1 else 0 end as c
                from alignment_from_known_source
              ) a
        '''),
        "mapqsFraction42": query_ratio(store, '''
          select sum(c) as numerator, count(c) as deliminator
            from (
                select case when mapq >= 42 then 1 else 0 end as c
                from alignment_from_known_source
              ) a
        '''),
        "query_lengthsAvg": query_one_number(store, '''
          select avg(query_length)
            from alignment_from_known_source
        '''),
        "query_lengthsFractionAtLeast60": query_ratio(store, '''
          select sum(c) as numerator, count(c) as deliminator
            from (
                select case when query_length >= 60 then 1 else 0 end as c
                from alignment_from_known_source
              ) a
        '''),

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
    parser.add_argument("--refdb-ncbi", type=str, action="store", dest="refdb_ncbi", help = "argument for ete.NCBITaxa")

    options=parser.parse_args(argv)


    logging.basicConfig(format='%(asctime)s: %(message)s')
    logger=logging.getLogger(__name__)
    if options.verbose:
        logger.setLevel(logging.DEBUG)

    main_loop(**vars(options), logger = logger )
    


if __name__ == '__main__':
    main()
