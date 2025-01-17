

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

rule align_merge:
    input:
        lambda w: expand_pandas(rules.chunk_sort_query.output.bam, ref_sets, allow_missing=True),
#        lambda w: expand(rules.chunk_sort_query.output.bam, zip, **ref_sets.to_dict("list"), allow_missing=True),
    output:
        bam = "results/align/merge/{sample}_{library}_{read_type_map}.bam",
    log:
        "logs/align/merge/{sample}_{library}_{read_type_map}.log"
    benchmark:
        "benchmarks/align/merge/{sample}_{library}_{read_type_map}.jsonl"
    params:
        extra = "-N -c -p",
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
use rule chunk_sort_query as align_sort_coord with:
    input:
        rules.align_merge.output.bam,
    output:
        bam = temp("temp/align/sort_coord/{sample}_{library}_{read_type_map}.bam"),
        idx = temp("temp/align/sort_coord/{sample}_{library}_{read_type_map}.bam.csi"),
    log:
        "logs/align/sort_coord/{sample}_{library}_{read_type_map}.log"
    benchmark:
        "benchmarks/align/sort_coord/{sample}_{library}_{read_type_map}.jsonl"
    params:
        mem_overhead_factor=0.4,
    resources:
        mem = lambda w, attempt, threads, input: f"{15 * threads * attempt} GiB",
        runtime = lambda w, attempt, input: f"{np.clip(0.0001 * input.size_mb + 1, 0.1, 20) * attempt} h",



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
