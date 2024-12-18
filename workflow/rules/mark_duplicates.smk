

if config["align"]["mark_duplicates"]["tool"] == "bamsormadup":
    rule shard_bamsormadup:
        input:
            lambda w: get_chunk_aln(w, "mark_duplicates"),
        output:
            bam = temp("temp/shard/mark_duplicates/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam"),
            index = temp("temp/shard/mark_duplicates/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam.bai"),
            metrics = "stats/shard/mark_duplicates/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.metrics.txt",
        log:
            "logs/shard/mark_duplicates/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.log"
        benchmark:
            "benchmarks/shard/mark_duplicates/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.jsonl"
        params:
            extra = lambda w: "SO=coordinate " + config["align"]["mark_duplicates"]["params"],
        threads: 10
        resources:
            mem = lambda w, attempt: f"{100 * attempt} GiB",
            runtime = lambda w, attempt: f"{5 * attempt} h",
        wrapper:
            f"{wrapper_ver}/bio/biobambam2/bamsormadup"


elif config["align"]["mark_duplicates"]["tool"] == "markdup":
    rule shard_markdup:
        input:
            bams = lambda w: get_chunk_aln(w, "mark_duplicates"),
        output:
            bam = temp("temp/shard/mark_duplicates/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam"),
            idx = temp("temp/shard/mark_duplicates/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam.csi"),
            metrics = "stats/shard/mark_duplicates/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.metrics.txt",
        log:
            "logs/shard/mark_duplicates/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.log"
        benchmark:
            "benchmarks/shard/mark_duplicates/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.jsonl"
        params:
            extra = config["align"]["mark_duplicates"]["params"],
        threads: 1
        resources:
            mem = lambda w, attempt: f"{40 * attempt} GiB",
            runtime = lambda w, attempt: f"{1 * attempt} d",
        wrapper:
            f"{wrapper_ver}/bio/samtools/markdup"


elif config["align"]["mark_duplicates"]["tool"] == "markduplicates":
    rule shard_markduplicates:
        input:
            bams = lambda w: get_chunk_aln(w, "mark_duplicates"),
        output:
            bam = temp("temp/shard/mark_duplicates/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam"),
            idx = temp("temp/shard/mark_duplicates/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam.bai"),
            metrics = "stats/shard/mark_duplicates/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.metrics.txt",
        log:
            "logs/shard/mark_duplicates/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.log"
        benchmark:
            "benchmarks/shard/mark_duplicates/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.jsonl"
        params:
            extra = "--CREATE_INDEX " + config["align"]["mark_duplicates"]["params"],
        threads: 1
        resources:
            mem = lambda w, attempt: f"{40 * attempt} GiB",
            runtime = lambda w, attempt: f"{1 * attempt} d",
        wrapper:
            f"{wrapper_ver}/bio/picard/markduplicates"


elif config["align"]["mark_duplicates"]["tool"] == "markduplicateswithmatecigar":
    rule shard_markduplicateswithmatecigar:
        input:
            bams = lambda w: get_chunk_aln(w, "mark_duplicates"),
        output:
            bam = temp("temp/shard/mark_duplicates/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam"),
            idx = temp("temp/shard/mark_duplicates/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam.bai"),
            metrics = "stats/shard/mark_duplicates/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.metrics.txt",
        log:
            "logs/shard/mark_duplicates/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.log"
        benchmark:
            "benchmarks/shard/mark_duplicates/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.jsonl"
        params:
            withmatecigar = True,
            extra = "--CREATE_INDEX " + config["align"]["mark_duplicates"]["params"],
        threads: 1
        resources:
            mem = lambda w, attempt: f"{40 * attempt} GiB",
            runtime = lambda w, attempt: f"{1 * attempt} d",
        wrapper:
            f"{wrapper_ver}/bio/picard/markduplicates"


elif config["align"]["mark_duplicates"]["tool"] == "markduplicatesspark":
    rule shard_markduplicatesspark:
        input:
            lambda w: get_chunk_aln(w, "mark_duplicates"),
        output:
            bam = temp("temp/shard/mark_duplicates/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam"),
            idx = temp("temp/shard/mark_duplicates/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam.bai"),
            metrics = "stats/shard/mark_duplicates/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.metrics.txt",
        log:
            "logs/shard/mark_duplicates/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.log"
        benchmark:
            "benchmarks/shard/mark_duplicates/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.jsonl"
        params:
            extra = "--create-output-bam-index true --create-output-bam-splitting-index false " + config["align"]["mark_duplicates"]["params"],
            spark_extra = "--conf 'spark.local.dir={}'".format(get_tmp()),
        threads: 10
        resources:
            mem = lambda w, attempt: f"{100 * attempt} GiB",
            runtime = lambda w, attempt: f"{1 * attempt} d",
            tmpdir = get_tmp(),
        wrapper:
            f"{wrapper_ver}/bio/gatk/markduplicatesspark"
