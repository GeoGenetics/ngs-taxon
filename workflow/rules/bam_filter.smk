def get_filter_aln(wildcards, rule):
    fall_through = False
    if rule == "sort_coord" or fall_through:
        if is_activated("bam_filter/filter"):
            return {"aln": rules.align_filter.output.bam}
        fall_through = True
    if rule == "align_filter" or fall_through:
        if is_activated("bam_filter/reassign"):
            return {"aln": rules.align_reassign.output.bam}
        fall_through = True
    if rule == "align_reassign" or fall_through:
        return {
            "aln": rules.align_sort_coord.output.bam,
            "idx": rules.align_sort_coord.output.idx,
        }
        fall_through = True
    raise ValueError(f"Invalid rule specified: {rule}")


#############
### RULES ###
#############


# Re-sorting, since it is required even when merging sorted BAM files (needs similar headers):
# https://www.biostars.org/p/251721/
# https://www.biostars.org/p/9509574/
use rule shard_sort_query as align_sort_coord with:
    input:
        rules.align_merge.output.bam,
    output:
        bam=temp("temp/aligns/sort_coord/{sample}_{library}_{read_type_map}.bam"),
        idx=temp("temp/aligns/sort_coord/{sample}_{library}_{read_type_map}.bam.csi"),
    log:
        "logs/aligns/sort_coord/{sample}_{library}_{read_type_map}.log",
    benchmark:
        "benchmarks/aligns/sort_coord/{sample}_{library}_{read_type_map}.jsonl"
    params:
        extra="",


rule align_reassign:
    input:
        unpack(lambda w: get_filter_aln(w, "align_reassign")),
    output:
        bam=temp("temp/aligns/reassign/{sample}_{library}_{read_type_map}.bam"),
    log:
        "logs/aligns/reassign/{sample}_{library}_{read_type_map}.log",
    benchmark:
        "benchmarks/aligns/reassign/{sample}_{library}_{read_type_map}.jsonl"
    params:
        #mem_overhead=0.5,
        extra=config["bam_filter"]["reassign"]["params"],
    conda:
        urlunparse(
            baseurl._replace(path=str(Path(baseurl.path) / "envs" / "bam_filter.yaml"))
        )
    threads: 10
    resources:
        mem=lambda w, attempt, input, threads: f"{np.clip(2* threads* input.size_mb/1024,10* threads,100* threads)* attempt} GiB",
        mem_mb=lambda w, attempt, input, threads: np.clip(
            2 * threads * input.size_mb / 1024, 10 * threads, 100 * threads
        )
        * attempt
        * 1024,
        runtime=lambda w, attempt: f"{2* attempt} d",
    shell:
        #"MEM_THREAD=`echo '{resources.mem_mb}*(1-{params.mem_overhead})/{threads}' | bc`M; "
        "filterBAM reassign --threads {threads} --max-memory {resources.mem_mb}M --bam {input.aln} {params.extra} --tmp-dir {resources.tmpdir} --out-bam {output.bam}  >{log} 2>&1"


rule align_filter:
    input:
        unpack(lambda w: get_filter_aln(w, "align_filter")),
    output:
        bam=temp("temp/aligns/filter/{sample}_{library}_{read_type_map}.bam"),
        read_len="stats/aligns/filter/{sample}_{library}_{read_type_map}.read-length-freqs.json",
        read_hits="stats/aligns/filter/{sample}_{library}_{read_type_map}.read-hits-count.tsv.gz",
        stats="stats/aligns/filter/{sample}_{library}_{read_type_map}.stats.tsv.gz",
        stats_filt="stats/aligns/filter/{sample}_{library}_{read_type_map}.stats-filtered.tsv.gz",
    log:
        "logs/aligns/filter/{sample}_{library}_{read_type_map}.log",
    benchmark:
        "benchmarks/aligns/filter/{sample}_{library}_{read_type_map}.jsonl"
    params:
        #mem_overhead=0.2,
        extra=config["bam_filter"]["filter"]["params"],
    conda:
        urlunparse(
            baseurl._replace(path=str(Path(baseurl.path) / "envs" / "bam_filter.yaml"))
        )
    threads: 10
    resources:
        mem=lambda w, attempt, input, threads: f"{np.clip(2* threads* input.size_mb/1024,10* threads,70* threads)* attempt} GiB",
        runtime=lambda w, attempt: f"{2* attempt} d",
    shell:
        #"MEM_THREAD=`echo '{resources.mem_mb}*(1-{params.mem_overhead})/{threads}' | bc`M; "
        "filterBAM filter --threads {threads} --bam {input.aln} {params.extra} --tmp-dir {resources.tmpdir} --bam-filtered {output.bam} --stats {output.stats} --stats-filtered {output.stats_filt} --read-length-freqs {output.read_len} --read-hits-count {output.read_hits} >{log} 2>&1"


rule align_lca:
    input:
        aln=rules.align_filter.output.bam,
        stats=rules.align_filter.output.stats_filt,
        nodes=config["taxonomy"]["nodes"],
        names=config["taxonomy"]["names"],
        acc2tax=[
            config["ref"][ref]["acc2taxid"]
            for ref in config["ref"]
            if config["ref"][ref]["acc2taxid"]
        ],
    output:
        stats=temp(
            "temp/aligns/bam_filter_lca/{sample}_{library}_{read_type_map}.summary"
        ),
    log:
        "logs/aligns/bam_filter_lca/{sample}_{library}_{read_type_map}.log",
    benchmark:
        "benchmarks/aligns/bam_filter_lca/{sample}_{library}_{read_type_map}.jsonl"
    params:
        #mem_overhead=0.2,
        extra=config["bam_filter"]["lca"]["params"],
    conda:
        urlunparse(
            baseurl._replace(path=str(Path(baseurl.path) / "envs" / "bam_filter.yaml"))
        )
    threads: 10
    resources:
        mem=lambda w, attempt, input, threads: f"{np.clip(2* threads* input.size_mb/1024,20* threads,70* threads)* attempt} GiB",
        runtime=lambda w, attempt: f"{2* attempt} d",
    shell:
        #"MEM_THREAD=`echo '{resources.mem_mb}*(1-{params.mem_overhead})/{threads}' | bc`M; "
        "filterBAM lca --threads {threads} --sort-memory 10G --bam {input.aln} --stats {input.stats} --names {input.names} --nodes {input.nodes} --acc2taxid <(cat {input.acc2tax}) {params.extra} --lca-summary {output.stats} >{log} 2>&1"


use rule shard_sort_query as align_sort_query with:
    input:
        unpack(lambda w: get_filter_aln(w, "sort_coord")),
    output:
        bam=temp("temp/aligns/sort_query/{sample}_{library}_{read_type_map}.bam"),
        idx=temp("temp/aligns/sort_query/{sample}_{library}_{read_type_map}.bam.csi"),
    log:
        "logs/aligns/sort_query/{sample}_{library}_{read_type_map}.log",
    benchmark:
        "benchmarks/aligns/sort_query/{sample}_{library}_{read_type_map}.jsonl"
