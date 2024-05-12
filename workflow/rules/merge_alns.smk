

def get_merge_aln(wildcards, rule):
    for case in switch(rule):
        if case("sort_name") or case("collate") or case("metadmg_damage") or case("bam_filter_lca"):
            if is_activated("bam_filter/filter"):
                return rules.bam_filter_filter.output.bam
        if case("bam_filter_filter"):
            if is_activated("bam_filter/reassign"):
                return rules.bam_filter_reassign.output.bam
        if case("bam_filter_reassign"):
            return rules.sort_merged.output.bam,
        if case():
            raise ValueError(f"Invalid merge_aln rule specified: {rule}")


#############
### RULES ###
#############

rule calmd:
    input:
        aln = lambda w: get_chunk_aln(w, "calmd"),
        ref = lambda w: config["ref"][w.ref]["path"],
        fai = lambda w: config["ref"][w.ref]["path"] + ".fai",
    output:
        bam = temp("temp/align/calmd/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam"),
#        idx = temp("temp/align/calmd/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam.csi"),
    log:
        "logs/align/calmd/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.log"
    benchmark:
        "benchmarks/align/calmd/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.jsonl"
    params:
        extra = config["align"]["calmd"]["params"],
    threads: 5
    resources:
        mem = lambda w, attempt: f"{10 * attempt} GiB",
        runtime = lambda w, attempt: f"{12 * attempt} h",
    wrapper:
        wrapper_ver + "/bio/samtools/calmd"


rule merge_alns:
    input:
        lambda w: expand_pandas(get_chunk_aln(w, "merge_alns"), bowtie_sets, allow_missing=True),
    output:
        bam = temp("temp/align/merge_alns/{sample}_{library}_{read_type_map}.bam"),
#        idx = temp("temp/align/merge_alns/{sample}_{library}_{read_type_map}.bam.csi"),
    log:
        "logs/align/merge_alns/{sample}_{library}_{read_type_map}.log"
    benchmark:
        "benchmarks/align/merge_alns/{sample}_{library}_{read_type_map}.jsonl"
    params:
        extra = "-c -p",
    threads: 3
    resources:
        mem = lambda w, attempt: f"{75 * attempt} GiB",
        runtime = lambda w, attempt: f"{2 * attempt} d",
        tmpdir = get_tmp(),
    wrapper:
        wrapper_ver + "/bio/samtools/merge"
 

# Re-sorting, since merging of BAM files needs similar headers:
# https://www.biostars.org/p/251721/
# https://www.biostars.org/p/9509574/
use rule sort_coord as sort_merged with:
    input:
        rules.merge_alns.output.bam,
    output:
        bam = "results/align/sort_merged/{sample}_{library}_{read_type_map}.bam",
        idx = "results/align/sort_merged/{sample}_{library}_{read_type_map}.bam.csi",
    log:
        "logs/align/sort_merged/{sample}_{library}_{read_type_map}.log"
    benchmark:
        "benchmarks/align/sort_merged/{sample}_{library}_{read_type_map}.jsonl"



##########
### QC ###
##########
rule samtools_stats:
    input:
        aln = rules.merge_alns.output.bam,
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
        wrapper_ver + "/bio/samtools/stats"
