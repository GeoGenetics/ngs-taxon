
rule bamsormadup:
    input:
        lambda w: get_chunk_aln(w, "mark_duplicates"),
    output:
        bam = temp("temp/align/bamsormadup/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam"),
        index = temp("temp/align/bamsormadup/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam.bai"),
        metrics = "stats/align/bamsormadup/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.metrics.txt",
    log:
        "logs/align/bamsormadup/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.log"
    benchmark:
        "benchmarks/align/bamsormadup/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.jsonl"
    params:
        extra = lambda w: "SO=coordinate " + config["align"]["mark_duplicates"]["params"],
    threads: 10
    resources:
        mem = lambda w, attempt: f"{100 * attempt} GiB",
        runtime = lambda w, attempt: f"{5 * attempt} h",
    wrapper:
        wrapper_ver + "/bio/biobambam2/bamsormadup"



rule markdup:
    input:
        bams = lambda w: get_chunk_aln(w, "mark_duplicates"),
    output:
        bam = temp("temp/align/markdup/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam"),
        idx = temp("temp/align/markdup/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam.csi"),
        metrics = "stats/align/markdup/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.metrics.txt",
    log:
        "logs/align/markdup/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.log"
    benchmark:
        "benchmarks/align/markdup/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.jsonl"
    params:
        extra = config["align"]["mark_duplicates"]["params"],
    threads: 1
    resources:
        mem = lambda w, attempt: f"{40 * attempt} GiB",
        runtime = lambda w, attempt: f"{1 * attempt} d",
    wrapper:
        wrapper_ver + "/bio/samtools/markdup"



rule markduplicates:
    input:
        bams = lambda w: get_chunk_aln(w, "mark_duplicates"),
    output:
        bam = temp("temp/align/markduplicates/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam"),
        idx = temp("temp/align/markduplicates/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam.bai"),
        metrics = "stats/align/markduplicates/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.metrics.txt",
    log:
        "logs/align/markduplicates/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.log"
    benchmark:
        "benchmarks/align/markduplicates/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.jsonl"
    params:
        extra = "--CREATE_INDEX " + config["align"]["mark_duplicates"]["params"],
    threads: 1
    resources:
        mem = lambda w, attempt: f"{40 * attempt} GiB",
        runtime = lambda w, attempt: f"{1 * attempt} d",
    wrapper:
        wrapper_ver + "/bio/picard/markduplicates"



use rule markduplicates as markduplicateswithmatecigar with:
    output:
        bam = temp("temp/align/markduplicateswithmatecigar/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam"),
        idx = temp("temp/align/markduplicateswithmatecigar/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam.bai"),
        metrics = "stats/align/markduplicateswithmatecigar/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.metrics.txt",
    log:
        "logs/align/markduplicateswithmatecigar/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.log"
    benchmark:
        "benchmarks/align/markduplicateswithmatecigar/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.jsonl"
    params:
        withmatecigar = True,
        extra = "--CREATE_INDEX " + config["align"]["mark_duplicates"]["params"],



rule markduplicatesspark:
    input:
        lambda w: get_chunk_aln(w, "mark_duplicates"),
    output:
        bam = temp("temp/align/markduplicatesspark/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam"),
        idx = temp("temp/align/markduplicatesspark/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam.bai"),
        metrics = "stats/align/markduplicatesspark/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.metrics.txt",
    log:
        "logs/align/markduplicatesspark/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.log"
    benchmark:
        "benchmarks/align/markduplicatesspark/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.jsonl"
    params:
        extra = "--create-output-bam-index true --create-output-bam-splitting-index false " + config["align"]["mark_duplicates"]["params"],
        spark_extra = "--conf 'spark.local.dir={}'".format(get_tmp()),
    threads: 10
    resources:
        mem = lambda w, attempt: f"{100 * attempt} GiB",
        runtime = lambda w, attempt: f"{1 * attempt} d",
        tmpdir = get_tmp(),
    wrapper:
        wrapper_ver + "/bio/gatk/markduplicatesspark"
