

#############
### RULES ###
#############

rule sort_coord:
    input:
        bam = rules.clean_header.output.bam,
    output:
        bam = temp("temp/align/sort_coord/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam"),
        idx = temp("temp/align/sort_coord/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam.csi"),
    log:
        "logs/align/sort_coord/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.log"
    benchmark:
        "benchmarks/align/sort_coord/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.tsv"
    threads: 4
    resources:
        mem = lambda w, attempt: f"{150 * attempt} GB",
        runtime = lambda w, attempt: f"{4 * attempt} d",
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
        "benchmarks/align/sort_name/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.tsv"
    params:
        extra = "-n",



rule bam_filter:
    input:
        unpack(lambda w: ext_dict(get_chunk_aln(w, "bam_filter", ext=["bam","bam.csi"]), keys=["aln", "idx"])),
    output:
        bam = temp("temp/align/bam_filter/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam"),
        knee = "stats/align/bam_filter/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.knee-plot.png",
        read_len = "stats/align/bam_filter/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.read-length-freqs.json",
        read_hits = "stats/align/bam_filter/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.read-hits-count.tsv.gz",
        stats = "stats/align/bam_filter/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.stats.tsv.gz",
        stats_filt = "stats/align/bam_filter/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.stats-filtered.tsv.gz",
    log:
        "logs/align/bam_filter/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.log"
    benchmark:
        "benchmarks/align/bam_filter/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.tsv"
    params:
        extra = check_cmd(config["align"]["bam_filter"]["params"], forbidden_args = ["-t", "--threads", "-m", "--low-memory", "--sort-memory", "-N", "--sort-by-name", "--disable-sort", "-r", "--reference-lengths", "--read-length-freqs", "--read-hits-count", "--only-stats", "--only-stats-filtered", "--plot", "-p", "--prefix"]),
    conda:
        base_dir / "envs" / "bam_filter.yaml"
    threads: 12
    resources:
        mem_mb = lambda w, attempt, input, threads: max(input.size_mb * 1.2 * threads, 50 * 1024) * attempt,
        runtime = lambda w, attempt: f"{1 * attempt} d",
    shell:
        "filterBAM --threads {threads} --sort-memory $(({resources.mem_mb}/{threads}))M --low-memory --disable-sort --bam {input.aln} --bam-index {input.idx} {params.extra} --tmp-dir {resources.tmpdir} --bam-filtered {output.bam} --stats {output.stats} --stats-filtered {output.stats_filt} --read-length-freqs {output.read_len} --read-hits-count {output.read_hits} --knee-plot {output.knee} >{log} 2>&1"



rule calmd:
    input:
        lambda w: get_chunk_aln(w, "calmd"),
        ref = lambda w: config["ref"][w.ref]["path"],
        fai = lambda w: config["ref"][w.ref]["path"] + ".fai",
    output:
        bam = temp("temp/align/calmd/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam"),
    log:
        "logs/align/calmd/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.log"
    benchmark:
        "benchmarks/align/calmd/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.tsv"
    params:
        extra = check_cmd(config["align"]["calmd"]["params"], forbidden_args = ["--threads", "-@", "-b"]),
    threads: 5
    resources:
        mem = lambda w, attempt: f"{10 * attempt} GB",
        runtime = lambda w, attempt: f"{12 * attempt} h",
    wrapper:
        wrapper_ver + "/bio/samtools/calmd"
