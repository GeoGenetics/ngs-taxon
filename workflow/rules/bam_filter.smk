

#############
### RULES ###
#############

rule bam_filter_reassign:
    input:
        aln = lambda w: get_merge_aln(w, "bam_filter_reassign"),
    output:
        bam = temp("temp/align/bam_filter/reassign/{sample}_{library}_{read_type_map}.bam"),
        csi = temp(touch("temp/align/bam_filter/reassign/{sample}_{library}_{read_type_map}.bam.csi")),
    log:
        "logs/align/bam_filter/reassign/{sample}_{library}_{read_type_map}.log"
    benchmark:
        "benchmarks/align/bam_filter/reassign/{sample}_{library}_{read_type_map}.jsonl"
    params:
        extra = config["bam_filter"]["reassign"]["params"],
    conda:
        base_dir / "envs" / "bam_filter.yaml"
    threads: 10
    resources:
        mem = lambda w, attempt, input, threads: f"{max(1.2 * input.size_mb * threads / 1024, 50) * attempt} GiB",
        mem_thread = lambda w, attempt, input, threads: max(1.2 * input.size_mb / 1024, 50 / threads) * attempt,
        runtime = lambda w, attempt: f"{10 * attempt} h",
#        tmpdir = "temp/large_temp",
    shell:
        "filterBAM reassign --threads {threads} --sort-memory {resources.mem_thread}G --bam {input.aln} {params.extra} --tmp-dir {resources.tmpdir} --out-bam {output.bam}  >{log} 2>&1"



rule bam_filter_filter:
    input:
        aln = lambda w: get_merge_aln(w, "bam_filter_filter"),
    output:
        bam = temp("temp/align/bam_filter/filter/{sample}_{library}_{read_type_map}.bam"),
#        csi = temp(touch("temp/align/bam_filter/filter/{sample}_{library}_{read_type_map}.bam.csi")),
        knee = touch("stats/align/bam_filter/filter/{sample}_{library}_{read_type_map}.knee-plot.png"),
        read_len = "stats/align/bam_filter/filter/{sample}_{library}_{read_type_map}.read-length-freqs.json",
        read_hits = "stats/align/bam_filter/filter/{sample}_{library}_{read_type_map}.read-hits-count.tsv.gz",
        stats = "stats/align/bam_filter/filter/{sample}_{library}_{read_type_map}.stats.tsv.gz",
        stats_filt = "stats/align/bam_filter/filter/{sample}_{library}_{read_type_map}.stats-filtered.tsv.gz",
    log:
        "logs/align/bam_filter/filter/{sample}_{library}_{read_type_map}.log"
    benchmark:
        "benchmarks/align/bam_filter/filter/{sample}_{library}_{read_type_map}.jsonl"
    params:
        extra = config["bam_filter"]["filter"]["params"],
    conda:
        base_dir / "envs" / "bam_filter.yaml"
    threads: 10
    resources:
        mem = lambda w, attempt, input, threads: f"{max(1.2 * input.size_mb * threads / 1024, 50) * attempt} GiB",
        mem_thread = lambda w, attempt, input, threads: max(1.2 * input.size_mb / 1024, 50 / threads) * attempt,
        runtime = lambda w, attempt: f"{1 * attempt} d",
#        tmpdir = "temp/large_temp",
    shell:
        "filterBAM filter --threads {threads} --sort-memory {resources.mem_thread}G --low-memory --bam {input.aln} {params.extra} --tmp-dir {resources.tmpdir} --bam-filtered {output.bam} --stats {output.stats} --stats-filtered {output.stats_filt} --read-length-freqs {output.read_len} --read-hits-count {output.read_hits} --knee-plot {output.knee} >{log} 2>&1"



rule bam_filter_lca:
    input:
        aln = lambda w: get_merge_aln(w, "bam_filter_lca"),
        stats = rules.bam_filter_filter.output.stats_filt,
        nodes = config["taxonomy"]["nodes"],
        names = config["taxonomy"]["names"],
        acc2tax = config["taxonomy"]["acc2taxid"],
    output:
        stats = temp("temp/align/bam_filter/lca/{sample}_{library}_{read_type_map}.summary"),
    log:
        "logs/align/bam_filter/lca/{sample}_{library}_{read_type_map}.log"
    benchmark:
        "benchmarks/align/bam_filter/lca/{sample}_{library}_{read_type_map}.jsonl"
    params:
        extra = config["bam_filter"]["lca"]["params"],
    conda:
        base_dir / "envs" / "bam_filter.yaml"
    threads: 10
    resources:
        mem = lambda w, attempt, input, threads: f"{max(1.2 * input.size_mb * threads / 1024, 50) * attempt} GiB",
        mem_thread = lambda w, attempt, input, threads: max(1.2 * input.size_mb / 1024, 50 / threads) * attempt,
        runtime = lambda w, attempt: f"{1 * attempt} d",
#        tmpdir = "temp/large_temp",
    shell:
        "filterBAM lca --threads {threads} --sort-memory {resources.mem_thread}G --bam {input.aln} --stats {input.stats} --names {input.names} --nodes {input.nodes} --acc2taxid {input.acc2tax} {params.extra} --lca-summary {output.stats} >{log} 2>&1"
