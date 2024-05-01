

#############
### RULES ###
#############

# https://bioinformatics.stackexchange.com/questions/18538/samtools-sort-most-efficient-memory-and-thread-settings-for-many-samples-on-a-c
rule sort_coord:
    input:
        bam = rules.clean_header.output.bam,
    output:
        bam = temp("temp/align/sort_coord/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam"),
        idx = temp("temp/align/sort_coord/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam.csi"),
    log:
        "logs/align/sort_coord/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.log"
    benchmark:
        "benchmarks/align/sort_coord/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.jsonl"
    threads: 8
    resources:
        mem = lambda w, attempt, threads: f"{10 * threads * attempt} GiB",
        runtime = lambda w, attempt: f"{1 * attempt} d",
    wrapper:
        wrapper_ver + "/bio/samtools/sort"


use rule sort_coord as sort_name with:
    input:
        lambda w: get_chunk_aln(w, "sort_name"),
    output:
        bam = temp("temp/align/sort_name/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam"),
    log:
        "logs/align/sort_name/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.log"
    benchmark:
        "benchmarks/align/sort_name/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.jsonl"
    params:
        extra = "-n",



rule reassign:
    input:
        unpack(lambda w: ext_dict(get_chunk_aln(w, "reassign", ext=["bam","bam.csi"]), keys=["aln", "idx"])),
    output:
        bam = temp("temp/align/reassign/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam"),
        csi = temp("temp/align/reassign/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam.csi"),
    log:
        "logs/align/reassign/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.log"
    benchmark:
        "benchmarks/align/reassign/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.jsonl"
    params:
        extra = config["align"]["reassign"]["params"],
    conda:
        base_dir / "envs" / "bam_filter.yaml"
    threads: 10
    resources:
        mem = lambda w, attempt, input, threads: f"{max(1.2 * input.size_mb * threads / 1024, 50) * attempt} GiB",
        runtime = lambda w, attempt: f"{10 * attempt} h",
#        tmpdir = "temp/large_temp",
    shell:
        "filterBAM reassign --threads {threads} --bam {input.aln} {params.extra} --tmp-dir {resources.tmpdir} --out-bam {output.bam}  >{log} 2>&1"



rule filter:
    input:
        unpack(lambda w: ext_dict(get_chunk_aln(w, "filter", ext=["bam", "bam.csi"]), keys=["aln", "idx"])),
    output:
        bam = temp("temp/align/filter/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam"),
#        csi = temp("temp/align/filter/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam.csi"),
        knee = touch("stats/align/filter/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.knee-plot.png"),
        read_len = "stats/align/filter/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.read-length-freqs.json",
        read_hits = "stats/align/filter/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.read-hits-count.tsv.gz",
        stats = "stats/align/filter/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.stats.tsv.gz",
        stats_filt = "stats/align/filter/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.stats-filtered.tsv.gz",
    log:
        "logs/align/filter/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.log"
    benchmark:
        "benchmarks/align/filter/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.jsonl"
    params:
        extra = config["align"]["filter"]["params"],
    conda:
        base_dir / "envs" / "bam_filter.yaml"
    threads: 10
    resources:
        mem = lambda w, attempt, input, threads: f"{max(1.2 * input.size_mb * threads / 1024, 50) * attempt} GiB",
        runtime = lambda w, attempt: f"{1 * attempt} d",
#        tmpdir = "temp/large_temp",
    shell:
        "filterBAM filter --threads {threads} --low-memory --bam {input.aln} {params.extra} --tmp-dir {resources.tmpdir} --bam-filtered {output.bam} --stats {output.stats} --stats-filtered {output.stats_filt} --read-length-freqs {output.read_len} --read-hits-count {output.read_hits} --knee-plot {output.knee} >{log} 2>&1"



rule calmd:
    input:
        aln = lambda w: get_chunk_aln(w, "calmd"),
        ref = lambda w: config["ref"][w.ref]["path"],
        fai = lambda w: config["ref"][w.ref]["path"] + ".fai",
    output:
        bam = temp("temp/align/calmd/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam"),
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
