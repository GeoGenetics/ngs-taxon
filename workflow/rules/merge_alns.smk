def get_merge_aln(wildcards, rule):
    fall_through = False
    if rule in ["metadmg_damage", "metadmg_lca"] or fall_through:
        if is_activated("bam_filter/filter") or is_activated("bam_filter/reassign"):
            return {"aln": rules.align_sort_query.output.bam}
        fall_through = True
    if rule == "align_filter" or fall_through:
        if is_activated("bam_filter/reassign"):
            return {"aln": rules.align_reassign.output.bam}
        fall_through = True
    if rule == "align_reassign" or fall_through:
        if is_activated("bam_filter/filter") or is_activated("bam_filter/reassign"):
            return {
                "aln": rules.align_sort_coord.output.bam,
                "idx": rules.align_sort_coord.output.idx,
            }
        else:
            return {"aln": rules.align_merge.output.bam}
        fall_through = True
    raise ValueError(f"Invalid rule specified: {rule}")


#############
### RULES ###
#############


rule align_merge:
    input:
        lambda w: expand_pandas(
            rules.shard_sort_query.output.bam, ref_sets, allow_missing=True
        ),
    # lambda w: expand(rules.shard_sort_query.output.bam, zip, **ref_sets.to_dict("list"), allow_missing=True),
    output:
        bam="results/aligns/merge/{sample}_{library}_{read_type_map}.bam",
    log:
        "logs/aligns/merge/{sample}_{library}_{read_type_map}.log",
    benchmark:
        "benchmarks/aligns/merge/{sample}_{library}_{read_type_map}.jsonl"
    params:
        extra="-n -c -p",
    threads: 4
    resources:
        mem=lambda w, input, attempt: f"{(0.2* input.size_gb+30)* attempt} GiB",
        runtime=lambda w, input, attempt: f"{(0.03* input.size_gb+0.1)* attempt} h",
    wrapper:
        f"{wrapper_ver}/bio/samtools/merge"


##########
### QC ###
##########


rule align_stats:
    input:
        aln=rules.align_merge.output.bam,
    output:
        txt="stats/aligns/samtools_stats/{sample}_{library}_{read_type_map}.txt",
    log:
        "logs/aligns/samtools_stats/{sample}_{library}_{read_type_map}.log",
    benchmark:
        "benchmarks/aligns/samtools_stats/{sample}_{library}_{read_type_map}.jsonl"
    threads: 4
    resources:
        mem=lambda w, input, attempt: f"{(0.04* input.size_gb+5)* attempt} GiB",
        runtime=lambda w, input, attempt: f"{(0.01* input.size_gb+0.1)* attempt} h",
    wrapper:
        f"{wrapper_ver}/bio/samtools/stats"
