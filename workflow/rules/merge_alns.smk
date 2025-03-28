def get_merge_aln(wildcards, rule):
    fall_through = False
    if rule == "sort_query" or fall_through:
        if is_activated("bam_filter/filter"):
            return {"aln": rules.align_filter.output.bam}
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
            rules.shard_sort_coord.output.bam, ref_sets, allow_missing=True
        ),
    #        lambda w: expand(rules.shard_sort_coord.output.bam, zip, **ref_sets.to_dict("list"), allow_missing=True),
    output:
        bam=temp("temp/aligns/merge/{sample}_{library}_{read_type_map}.bam"),
    log:
        "logs/aligns/merge/{sample}_{library}_{read_type_map}.log",
    benchmark:
        "benchmarks/aligns/merge/{sample}_{library}_{read_type_map}.jsonl"
    params:
        extra="-c -p",
    threads: 4
    resources:
        mem=lambda w, input, attempt: f"{(7e-5* input.size_mb+40)* attempt} GiB",
        runtime=lambda w, input, attempt: f"{max(5e-5* input.size_mb* attempt,0.5)} h",
        tmpdir=get_tmp(),
    wrapper:
        f"{wrapper_ver}/bio/samtools/merge"


# Re-sorting, since it is required even when merging sorted BAM files (needs similar headers):
# https://www.biostars.org/p/251721/
# https://www.biostars.org/p/9509574/
use rule shard_sort_coord as align_sort_coord with:
    input:
        rules.align_merge.output.bam,
    output:
        bam=temp("temp/aligns/sort_coord/{sample}_{library}_{read_type_map}.bam"),
        idx=temp("temp/aligns/sort_coord/{sample}_{library}_{read_type_map}.bam.csi"),
    log:
        "logs/aligns/sort_coord/{sample}_{library}_{read_type_map}.log",
    benchmark:
        "benchmarks/aligns/sort_coord/{sample}_{library}_{read_type_map}.jsonl"
    resources:
        mem=lambda w, attempt, threads, input: f"{15* threads* attempt} GiB",
        runtime=lambda w, attempt, input: f"{np.clip(0.0001* input.size_mb+1,0.1,20)* attempt} h",
