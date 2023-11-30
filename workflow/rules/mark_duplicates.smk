
rule bamsormadup:
    input:
        lambda w: get_chunk_aln(w, "mark_duplicates"),
    output:
        cram = temp("temp/align/bamsormadup/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam"),
        index = temp("temp/align/bamsormadup/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bai"),
        metrics = "stats/align/bamsormadup/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.metrics.txt",
    log:
        "logs/align/bamsormadup/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.log"
    benchmark:
        "benchmarks/align/bamsormadup/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.tsv"
    params:
        extra = lambda w: "SO=coordinate " +
            check_cmd(config["align"]["mark_duplicates"]["params"], forbidden_args = ["threads", "inputformat", "outputformat", "SO", "tmpfile"]),
    threads: 10
    resources:
        mem = lambda w, attempt: f"{100 * attempt} GB",
        runtime = lambda w, attempt: f"{5 * attempt} h",
    wrapper:
        wrapper_ver + "/bio/biobambam2/bamsormadup"



rule markduplicates:
    input:
        bams = lambda w: get_chunk_aln(w, "mark_duplicates"),
    output:
        bam = temp("temp/align/markduplicates/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam"),
        metrics = "stats/align/markduplicates/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.metrics.txt",
    log:
        "logs/align/markduplicates/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.log"
    benchmark:
        "benchmarks/align/markduplicates/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.tsv"
    params:
        extra = check_cmd(config["align"]["mark_duplicates"]["params"], forbidden_args = ["--INPUT", "--TMP_DIR", "--OUTPUT", "--METRICS_FILE"]),
    threads: 1
    resources:
        mem = lambda w, attempt: f"{40 * attempt} GB",
        runtime = lambda w, attempt: f"{1 * attempt} d",
    wrapper:
        wrapper_ver + "/bio/picard/markduplicates"



use rule markduplicates as markduplicateswithmatecigar with:
    output:
        bam = temp("temp/align/markduplicateswithmatecigar/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam"),
        metrics = "stats/align/markduplicateswithmatecigar/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.metrics.txt",
    log:
        "logs/align/markduplicateswithmatecigar/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.log"
    benchmark:
        "benchmarks/align/markduplicateswithmatecigar/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.tsv"
    params:
        withmatecigar = True,
#        extra = check_cmd(config["align"]["mark_duplicates"]["params"], forbidden_args = ["--INPUT", "--TMP_DIR", "--OUTPUT", "--METRICS_FILE"]),



rule markduplicatesspark:
    input:
        lambda w: get_chunk_aln(w, "mark_duplicates"),
    output:
        bam = temp("temp/align/markduplicatesspark/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam"),
#        bai = temp("temp/align/markduplicatesspark/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.crai"),
        metrics = "stats/align/markduplicatesspark/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.metrics.txt",
    log:
        "logs/align/markduplicatesspark/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.log"
    benchmark:
        "benchmarks/align/markduplicatesspark/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.tsv"
    params:
        extra = "--create-output-bam-index false --create-output-bam-splitting-index false " +
            check_cmd(config["align"]["mark_duplicates"]["params"], forbidden_args = ["--input", "--metrics-file", "--create-output-bam-index", "--create-output-bam-splitting-index", "--tmp-dir", "--output"]),
        spark_extra = "--conf 'spark.local.dir={}'".format(get_tmp()),
    threads: 10
    resources:
        mem = lambda w, attempt: f"{100 * attempt} GB",
        runtime = lambda w, attempt: f"{1 * attempt} d",
        tmpdir = get_tmp(),
    wrapper:
        wrapper_ver + "/bio/gatk/markduplicatesspark"
