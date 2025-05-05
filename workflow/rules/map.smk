#################
### FUNCTIONS ###
#################


refs = pd.DataFrame.from_dict(config["ref"]).transpose()
ref_sets = refs.reset_index().rename(columns={"index": "ref", "n_shards": "tot_shards"})
ref_sets["read_type_map"] = "collapsed"
ref_sets["n_shard"] = [range(1, tot + 1) for tot in ref_sets["tot_shards"]]
ref_sets = ref_sets.join(pd.json_normalize(ref_sets["map"])).drop("map", axis="columns")


wildcard_constraints:
    read_type_map="|".join(["pe", "se", "singleton", "collapsed"]),


# https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups
def get_read_group(wildcards):
    """Denote sample name and platform in read group."""

    extra_info = {
        "center": "CN",
        "date": "DT",
        "platform": "PL",
        "platform_model": "PM",
    }

    # Get flowcell
    flowcell = (
        units.loc[(wildcards.sample, wildcards.library, slice(None)), "flowcell"]
        .drop_duplicates()
        .item()
    )
    rg_info = [
        f"{wildcards.sample}_{wildcards.library}_{flowcell}_{wildcards.ref}{wildcards.n_shard}o{wildcards.tot_shards}",
        f"SM:{wildcards.sample}",
        f"LB:{wildcards.library}",
        f"PU:{flowcell}",
    ]
    # Add extra info to RG
    for key, abv in extra_info.items():
        value = units.loc[(wildcards.sample, wildcards.library, slice(None))].get(key)
        if value is not None and value.any():
            rg_info.append(f"{abv}:{value.drop_duplicates().item()}")

    return rg_info


#############
### RULES ###
#############


def get_data(wildcards):
    return expand(
        "results/reads/low_complexity/{sample}_{library}_{read_type_trim}.fastq.gz",
        read_type_trim=get_read_type_trim(wildcards.read_type_map),
        allow_missing=True,
    )


def get_index(wildcards):
    if config["ref"][wildcards.ref]["map"].get("bt2l", True):
        return multiext(
            config["ref"][wildcards.ref]["path"],
            ".1.bt2l",
            ".2.bt2l",
            ".3.bt2l",
            ".4.bt2l",
            ".rev.1.bt2l",
            ".rev.2.bt2l",
        )
    else:
        return multiext(
            config["ref"][wildcards.ref]["path"],
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        )


rule bowtie2:
    input:
        sample=get_data,
        idx=get_index,
    output:
        bam=temp(
            "temp/shards/bowtie2/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.bam"
        ),
    log:
        "logs/shards/bowtie2/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.log",
    benchmark:
        "benchmarks/shards/bowtie2/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.jsonl"
    params:
        extra=lambda w: f"""--time --rg-id '{"' --rg '".join(get_read_group(w))}' --rg 'PG:bowtie2' """
        + config["ref"][w.ref]["map"]["params"],
    threads: 20
    resources:
        mem=lambda w, attempt, input: "{} GiB".format(
            (0.7 * sum(Path(f).stat().st_size for f in input.idx) / 1024**3 + 50)
            * attempt
        ),
        runtime=lambda w, attempt: f"{10* attempt} h",
        slurm_extra="--nice=2000",
    wrapper:
        f"{wrapper_ver}/bio/bowtie2/align"


rule bwa_aln:
    input:
        fastq=get_data,
        idx=lambda w: multiext(
            config["ref"][w.ref]["path"], ".amb", ".ann", ".bwt", ".pac", ".sa"
        ),
    output:
        sai=temp(
            "temp/shards/bwa_aln/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.sai"
        ),
    log:
        "logs/shards/bwa_aln/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.log",
    benchmark:
        "benchmarks/shards/bwa_aln/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.jsonl"
    params:
        extra=lambda w: config["ref"][w.ref]["map"]["params"],
    threads: 10
    resources:
        mem=lambda w, attempt, input: "{} GiB".format(
            (2 * sum(Path(f).stat().st_size for f in input.idx) / 1024**3 + 50)
            * attempt
        ),
        runtime=lambda w, attempt: f"{1* attempt} d",
        slurm_extra="--nice=2000",
    wrapper:
        f"{wrapper_ver}/bio/bwa/aln"


rule bwa_samxe:
    input:
        fastq=get_data,
        sai=rules.bwa_aln.output.sai,
        idx=rules.bwa_aln.input.idx,
    output:
        bam=temp(
            "temp/shards/bwa_aln/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.bam"
        ),
    log:
        "logs/shards/bwa_aln/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.samxe.log",
    benchmark:
        "benchmarks/shards/bwa_aln/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.samxe.jsonl"
    params:
        extra=lambda w: f"""-r '@RG\tID:{"\\t".join(get_read_group(w))}\\tPG:bwa_aln' """,
        sort="samtools",
    threads: 1
    resources:
        mem=lambda w, attempt, input: "{} GiB".format(
            (2 * sum(Path(f).stat().st_size for f in input.idx) / 1024**3 + 50)
            * attempt
        ),
        runtime=lambda w, attempt: f"{1* attempt} d",
    wrapper:
        f"{wrapper_ver}/bio/bwa/samxe"


rule shard_count_alns:
    input:
        bam=lambda w: expand(
            "temp/shards/{tool}/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.bam",
            tool=config["ref"][w.ref]["map"]["tool"],
            allow_missing=True,
        ),
    output:
        counts=temp(
            "temp/shards/count_alns/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.tsv"
        ),
    log:
        "logs/shards/count_alns/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.log",
    benchmark:
        "benchmarks/shards/count_alns/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.jsonl"
    threads: 1
    resources:
        mem=lambda w, attempt: f"{5* attempt} GiB",
        runtime=lambda w, attempt: f"{2* attempt} h",
    shell:
        """
        (samtools view {input.bam} | awk 'BEGIN{{print "read_id\tn_aligns"}} {{x[$1]++}} END{{for(read_id in x){{print read_id"\t"x[read_id]}}}}') > {output.counts} 2> {log}
        """


rule shard_saturated_reads_get:
    input:
        counts=expand_pandas(
            rules.shard_count_alns.output.counts, ref_sets, allow_missing=True
        ),
    output:
        read_id=temp(
            "temp/shards/saturated_reads/get/{sample}_{library}_{read_type_map}.tsv"
        ),
    log:
        "logs/shards/saturated_reads/get/{sample}_{library}_{read_type_map}.log",
    benchmark:
        "benchmarks/shards/saturated_reads/get/{sample}_{library}_{read_type_map}.jsonl"
    params:
        extra="--headerless-tsv-output cat then filter '$n_aligns >= {}' then uniq -g read_id then cut -f read_id".format(
            config["filter"]["saturated_reads"]["n_alns"]
        ),
    threads: 4
    resources:
        mem=lambda w, attempt: f"{1* attempt} GiB",
        runtime=lambda w, attempt: f"{1* attempt} h",
    wrapper:
        f"{wrapper_ver}/utils/miller"


rule shard_saturated_reads_remove:
    input:
        bam=rules.shard_count_alns.input.bam,
        read_id=rules.shard_saturated_reads_get.output.read_id,
    output:
        bam=temp(
            "temp/shards/saturated_reads/remove/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.bam"
        ),
    log:
        "logs/shards/saturated_reads/remove/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.log",
    benchmark:
        "benchmarks/shards/saturated_reads/remove/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.jsonl"
    params:
        extra=lambda w, input: f"--qname-file ^{input.read_id}",
    threads: 2
    resources:
        mem=lambda w, attempt: f"{10* attempt} GiB",
        runtime=lambda w, attempt: f"{30* attempt} m",
    wrapper:
        f"{wrapper_ver}/bio/samtools/view"


rule shard_saturated_reads_extract:
    input:
        fastq=get_data,
        read_id=rules.shard_saturated_reads_get.output.read_id,
    output:
        fq="results/reads/saturated/{sample}_{library}_{read_type_map}.fastq.gz",
    log:
        "logs/shards/saturated_reads/extract/{sample}_{library}_{read_type_map}.log",
    benchmark:
        "benchmarks/shards/saturated_reads/extract/{sample}_{library}_{read_type_map}.jsonl"
    params:
        command="subseq",
        compress_lvl=9,
        extra="",
    threads: 4
    resources:
        mem=lambda w, attempt: f"{1* attempt} GiB",
        runtime=lambda w, attempt: f"{30* attempt} m",
    wrapper:
        f"{wrapper_ver}/bio/seqtk"


rule shard_clean_header:
    input:
        bam=(
            rules.shard_saturated_reads_remove.output.bam
            if is_activated("filter/saturated_reads")
            else rules.shard_count_aln.input.bam
        ),
    output:
        bam=temp(
            "temp/shards/clean_header/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.bam"
        ),
    log:
        "logs/shards/clean_header/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.log",
    benchmark:
        "benchmarks/shards/clean_header/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.jsonl"
    threads: 2
    resources:
        mem=lambda w, attempt: f"{20* attempt} GiB",
        runtime=lambda w, attempt: f"{2* attempt} h",
    shell:
        "/projects/caeg/apps/metaDMG-cpp/misc/compressbam --threads {threads} --input {input.bam} --output {output.bam} >{log} 2>&1"


# https://bioinformatics.stackexchange.com/questions/18538/samtools-sort-most-efficient-memory-and-thread-settings-for-many-samples-on-a-c
rule shard_sort_query:
    input:
        bam=rules.shard_clean_header.output.bam,
    output:
        bam=temp(
            "temp/shards/sort_query/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.bam"
        ),
    log:
        "logs/shards/sort_query/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.log",
    benchmark:
        "benchmarks/shards/sort_query/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.jsonl"
    params:
        extra="-n",
        mem_overhead_factor=0.1,
    threads: 8
    resources:
        mem=lambda w, attempt, threads: f"{10* threads* attempt} GiB",
        runtime=lambda w, attempt, input: f"{max(0.0001* input.size_mb+1,0.1)* attempt} h",
    wrapper:
        f"{wrapper_ver}/bio/samtools/sort"
