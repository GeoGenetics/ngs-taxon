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
        bam="<results>/<aligns>/merge/{sample}_{library}_{read_type_map}.bam",
    log:
        "<logs>/<aligns>/merge/{sample}_{library}_{read_type_map}.log",
    benchmark:
        "<benchmarks>/<aligns>/merge/{sample}_{library}_{read_type_map}.jsonl"
    threads: 3
    resources:
        mem=lambda w, input, attempt: f"{(0.2* input.size_gb+30)* attempt} GiB",
        runtime=lambda w, input, attempt: f"{(0.03* input.size_gb+1)* attempt} h",
    params:
        extra="-n -c -p",
    wrapper:
        "v7.9.1/bio/samtools/merge"


# https://bioinformatics.stackexchange.com/questions/18538/samtools-sort-most-efficient-memory-and-thread-settings-for-many-samples-on-a-c
rule align_sort_taxon:
    input:
        bam=rules.align_merge.output.bam,
    output:
        bam=temp("<temp>/<aligns>/sort_taxon/{sample}_{library}_{read_type_map}.bam"),
    log:
        "<logs>/<aligns>/sort_taxon/{sample}_{library}_{read_type_map}.log",
    benchmark:
        "<benchmarks>/<aligns>/sort_taxon/{sample}_{library}_{read_type_map}.jsonl"
    threads: 6
    resources:
        mem=lambda w, input, attempt: f"{(10* input.size_gb+20)* attempt} GiB",
        runtime=lambda w, input, attempt: f"{(0.02* input.size_gb+1)* attempt} h",
    params:
        extra="-t XR",
        mem_overhead_factor=0.2,
    wrapper:
        "v9.4.2/bio/samtools/sort"


##########
### QC ###
##########


rule align_samtools_stats:
    input:
        bam=rules.align_merge.output.bam,
    output:
        "<stats>/<aligns>/samtools/stats/{sample}_{library}_{read_type_map}.txt",
    log:
        "<logs>/<aligns>/samtools/stats/{sample}_{library}_{read_type_map}.log",
    benchmark:
        "<benchmarks>/<aligns>/samtools/stats/{sample}_{library}_{read_type_map}.jsonl"
    threads: 2
    resources:
        mem=lambda w, input, attempt: f"{5* attempt} GiB",
        runtime=lambda w, input, attempt: f"{(0.02* input.size_gb+0.5)* attempt} h",
    wrapper:
        "v8.1.1/bio/samtools/stats"


rule align_unicorn_taxstats:
    input:
        bam=rules.align_sort_taxon.output.bam,
        nodes=config["taxonomy"]["nodes"],
        names=config["taxonomy"]["names"],
    output:
        stats="<stats>/<aligns>/unicorn/taxstats/{sample}_{library}_{read_type_map}.tsv",
    log:
        "<logs>/<aligns>/unicorn/taxstats/{sample}_{library}_{read_type_map}.log",
    benchmark:
        "<benchmarks>/<aligns>/unicorn/taxstats/{sample}_{library}_{read_type_map}.jsonl"
    conda:
        urlunparse(
            baseurl._replace(
                path=str(Path(baseurl.path) / "envs" / "enhjoerning.yaml")
            )
        )
    threads: 4
    resources:
        mem=lambda w, input, attempt: f"{(4* input.size_gb+50)* attempt} GiB",
        runtime=lambda w, input, attempt: f"{(0.05* input.size_gb+0.1)* attempt} h",
    params:
        extra=config["unicorn"]["taxstats"]["params"],
    shell:
        "unicorn taxstats --threads {threads} -b {input.bam} --names {input.names} --nodes {input.nodes} --qsize 10000 {params.extra} --outstat {output.stats} > {log} 2>&1"
