

def get_merge_aln(wildcards, rule):
    fall_through = False
    if rule in ["sort_query", "collate", "metadmg_damage", "bam_filter_lca"] or fall_through:
        if is_activated("bam_filter/filter"):
            return {"aln": rules.align_filter.output.bam, "idx": rules.align_filter.output.idx}
        fall_through = True
    if rule == "align_filter" or fall_through:
        if is_activated("bam_filter/reassign"):
            return {"aln": rules.align_reassign.output.bam, "idx": rules.align_reassign.output.idx}
        fall_through = True
    if rule == "align_reassign" or fall_through:
        if is_activated("bam_filter/filter"):
            return {"aln": rules.align_sort_coord.output.bam, "idx": rules.align_sort_coord.output.idx}
        else:
            return {"aln": rules.align_merge.output.bam}
        fall_through = True
    raise ValueError(f"Invalid rule specified: {rule}")



#############
### RULES ###
#############

rule shard_calmd:
    input:
        aln = lambda w: get_chunk_aln(w, "calmd"),
        ref = lambda w: config["ref"][w.ref]["path"],
        fai = lambda w: config["ref"][w.ref]["path"] + ".fai",
    output:
        bam = temp("temp/shard/calmd/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam"),
#        idx = temp("temp/shard/calmd/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam.csi"),
    log:
        "logs/shard/calmd/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.log"
    benchmark:
        "benchmarks/shard/calmd/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.jsonl"
    params:
        extra = config["align"]["calmd"]["params"],
    threads: 5
    resources:
        mem = lambda w, attempt: f"{10 * attempt} GiB",
        runtime = lambda w, attempt: f"{12 * attempt} h",
    wrapper:
        f"{wrapper_ver}/bio/samtools/calmd"


rule align_merge:
    input:
        lambda w: expand_pandas(get_chunk_aln(w, "merge"), ref_sets, allow_missing=True),
#        lambda w: expand(get_chunk_aln(w, "merge"), zip, **ref_sets.to_dict("list"), allow_missing=True),
    output:
        bam = temp("temp/align/merge/{sample}_{library}_{read_type_map}.bam"),
#        idx = temp("temp/align/merge/{sample}_{library}_{read_type_map}.bam.csi"),
    log:
        "logs/align/merge/{sample}_{library}_{read_type_map}.log"
    benchmark:
        "benchmarks/align/merge/{sample}_{library}_{read_type_map}.jsonl"
    params:
        extra = "-c -p",
    threads: 3
    resources:
        mem = lambda w, input, attempt: f"{(7e-5 * input.size_mb + 40) * attempt} GiB",
        runtime = lambda w, input, attempt: f"{max(5e-5 * input.size_mb * attempt, 0.5)} h",
        tmpdir = get_tmp(),
    wrapper:
        f"{wrapper_ver}/bio/samtools/merge"



# Re-sorting, since it is required even when merging sorted BAM files (needs similar headers):
# https://www.biostars.org/p/251721/
# https://www.biostars.org/p/9509574/
use rule shard_sort_coord as align_sort_coord with:
    input:
        rules.align_merge.output.bam,
    output:
        bam = temp("temp/align/sort_coord/{sample}_{library}_{read_type_map}.bam"),
        idx = temp("temp/align/sort_coord/{sample}_{library}_{read_type_map}.bam.csi"),
    log:
        "logs/align/sort_coord/{sample}_{library}_{read_type_map}.log"
    benchmark:
        "benchmarks/align/sort_coord/{sample}_{library}_{read_type_map}.jsonl"



##########
### QC ###
##########
rule align_samtools_stats:
    input:
        aln = rules.align_sort_coord.output.bam,
    output:
        txt = "stats/align/samtools_stats/{sample}_{library}_{read_type_map}.txt",
    log:
        "logs/align/samtools_stats/{sample}_{library}_{read_type_map}.log",
    benchmark:
        "benchmarks/align/samtools_stats/{sample}_{library}_{read_type_map}.jsonl",
    threads: 4
    resources:
        mem = lambda w, attempt: f"{10 * attempt} GiB",
        runtime = lambda w, attempt: f"{10 * attempt} h",
    wrapper:
        f"{wrapper_ver}/bio/samtools/stats"
