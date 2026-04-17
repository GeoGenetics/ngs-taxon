#################
### FUNCTIONS ###
#################
import pandas as pd

refs = pd.DataFrame.from_dict(config["ref"]).transpose()
ref_sets = refs.reset_index().rename(columns={"index": "ref", "n_shards": "tot_shards"})
ref_sets["read_type_map"] = "collapsed"
ref_sets["n_shard"] = [range(1, tot + 1) for tot in ref_sets["tot_shards"]]
ref_sets = ref_sets.join(pd.json_normalize(ref_sets["map"])).drop("map", axis="columns")


wildcard_constraints:
    read_type_map="|".join(["pe", "se", "singleton", "collapsed"]),


# https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups
def read_group_get(wildcards, tool):
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
    rg_info = {
        "ID": f"{wildcards.sample}_{wildcards.library}_{flowcell}_{wildcards.ref}{wildcards.n_shard}o{wildcards.tot_shards}",
        "SM": wildcards.sample,
        "LB": wildcards.library,
        "PU": flowcell,
        "PG": tool,
    }
    # Add extra info to RG
    for key, abv in extra_info.items():
        value = units.loc[(wildcards.sample, wildcards.library, slice(None))].get(key)
        if value is not None and value.any():
            rg_info[abv] = value.drop_duplicates().item()

    return rg_info


def read_group_merge(wildcards, tool):
    rg_info = read_group_get(wildcards, tool)
    # Get RG ID
    rg_fmt = [f"{k}:{{{k}}}" for k, v in rg_info.items() if k != "ID"]
    # Define format string
    if tool == "bowtie2":
        rg_fmt = "{ID}' --rg '" + "' --rg '".join(rg_fmt)
    elif tool == "bwa_aln":
        rg_fmt = "@RG\\tID:{ID}\\t" + "\\t".join(rg_fmt)
    elif tool == "dragen":
        rg_fmt = " ".join(
            [f"--RG{k} '{{{k}}}'" for k, v in rg_info.items() if k != "PG"]
        )
    else:
        raise ValueError(f"Tool not supported: {tool}")

    # Format and return
    return rg_fmt.format(**rg_info)


#############
### RULES ###
#############


def get_data(wildcards):
    read_type_trim = {"pe": ["R1", "R2"], "se": ["R"]}
    return local(
        expand(
            "<results>/reads/low_complexity/{sample}_{library}_{read_type_trim}.fastq.gz",
            read_type_trim=read_type_trim.get(
                wildcards.read_type_map, wildcards.read_type_map
            ),
            allow_missing=True,
        )
    )


def get_bowtie2_index(wildcards):
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


rule shard_bowtie2:
    input:
        sample=get_data,
        idx=get_bowtie2_index,
    output:
        bam=temp(
            "<temp>/<shards>/bowtie2/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.bam"
        ),
    log:
        "<logs>/<shards>/bowtie2/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.log",
    benchmark:
        "<benchmarks>/<shards>/bowtie2/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.jsonl"
    threads: 16
    resources:
        mem=lambda w, input, attempt: "{} GiB".format(
            (0.7 * sum(Path(f).stat().st_size for f in input.idx) / 1024**3 + 40)
            * attempt
        ),
        runtime=lambda w, input, attempt: f"{(Path(input.sample[0]).stat().st_size/1024**3+8)* attempt} h",
        slurm_extra="--extra-node-info 1",
    params:
        extra=lambda w: f"""--time --rg-id '{read_group_merge(w, "bowtie2")}' """
        + config["ref"][w.ref]["map"]["params"],
    wrapper:
        "v7.9.1/bio/bowtie2/align"


rule shard_bwa_aln:
    input:
        fastq=get_data,
        idx=lambda w: multiext(
            config["ref"][w.ref]["path"], ".amb", ".ann", ".bwt", ".pac", ".sa"
        ),
    output:
        sai=temp(
            "<temp>/<shards>/bwa_aln/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.sai"
        ),
    log:
        "<logs>/<shards>/bwa_aln/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.log",
    benchmark:
        "<benchmarks>/<shards>/bwa_aln/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.jsonl"
    threads: 10
    resources:
        mem=lambda w, attempt, input: "{} GiB".format(
            (2 * sum(Path(f).stat().st_size for f in input.idx) / 1024**3 + 50)
            * attempt
        ),
        runtime=lambda w, attempt: f"{10* attempt} h",
    params:
        extra=lambda w: config["ref"][w.ref]["map"]["params"],
    wrapper:
        "v5.10.0/bio/bwa/aln"


rule shard_bwa_samxe:
    input:
        fastq=rules.shard_bwa_aln.input.fastq,
        sai=rules.shard_bwa_aln.output.sai,
        idx=rules.shard_bwa_aln.input.idx,
    output:
        bam=temp(
            "<temp>/<shards>/bwa_aln/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.bam"
        ),
    log:
        "<logs>/<shards>/bwa_aln/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.samxe.log",
    benchmark:
        "<benchmarks>/<shards>/bwa_aln/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.samxe.jsonl"
    threads: 1
    resources:
        mem=lambda w, attempt, input: "{} GiB".format(
            (2 * sum(Path(f).stat().st_size for f in input.idx) / 1024**3 + 50)
            * attempt
        ),
        runtime=lambda w, attempt: f"{10* attempt} h",
    params:
        extra=lambda w: f"""-r '{read_group_merge(w, "bwa_aln")}' """,
        sort="samtools",
    wrapper:
        "v7.9.1/bio/bwa/samxe"


rule shard_dragen:
    input:
        sample=get_data,
        #idx=lambda w: multiext(
        #    config["ref"][w.ref]["path"] + "/",
        #    "dragen-replay.json",
        #    "dragen.metrics.json",
        #    "dragen.time_metrics.csv",
        #    "hash_table.cfg",
        #    "hash_table.cfg.bin",
        #    "hash_table.cmp",
        #    "hash_table_stats.txt",
        #    "ref_index.bin",
        #    "reference.bin",
        #    "str_table.bin",
        #),
    output:
        bam=temp(
            "<temp>/<shards>/dragen/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.bam"
        ),
        md5=temp(
            "<temp>/<shards>/dragen/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.bam.md5sum"
        ),
        replay="<stats>/<shards>/dragen/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}-replay.json",
        stats_json="<stats>/<shards>/dragen/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.json",
        stats_fastqc="<stats>/<shards>/dragen/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.fastqc.csv",
        stats_trimmer="<stats>/<shards>/dragen/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.trim.csv",
        stats_map=touch(
            "<stats>/<shards>/dragen/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.map.csv"
        ),
        stats_ploidy="<stats>/<shards>/dragen/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.ploidy.csv",
        stats_time="<stats>/<shards>/dragen/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.time.csv",
    log:
        "<logs>/<shards>/dragen/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.log",
    benchmark:
        "<benchmarks>/<shards>/dragen/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.jsonl"
    threads: 10
    resources:
        mem=lambda w, attempt: f"{100* attempt} GiB",
        runtime=lambda w, attempt: f"{1* attempt} h",
        constraint=lambda w: f"idx{int(w.n_shard)%2:02d}",
    params:
        version="4.4.6",
        idx_dir=lambda w: expand(config["ref"][w.ref]["path"], **w),
        #idx_dir=lambda w, input: os.path.commonprefix(input.idx),
        output_dir=lambda w, output: str(Path(output.bam).parent),
        output_prefix=lambda w, output: Path(output.bam).with_suffix("").name,
        stats_dir=lambda w, output: str(Path(output.stats_time).parent),
        rg=lambda w: read_group_merge(w, "dragen"),
        extra=lambda w: config["ref"][w.ref]["map"]["params"],
    shell:
        """(
        /opt/dragen/{params.version}/bin/dragen --num-threads {threads} -1 {input.sample} -r {params.idx_dir} {params.rg} {params.extra} --generate-zs-tags=true --enable-vcf-indexing=false --enable-bam-indexing=false --qc-coverage-reports-wgs=false --generate-md-tags=true --enable-sort=false --dump-map-align-registers=true --force --output-directory {params.output_dir} --output-file-prefix {params.output_prefix} &&
        mv --verbose {params.output_dir}/{params.output_prefix}-replay.json {output.replay} &&
        mv --verbose {params.output_dir}/{params.output_prefix}.metrics.json {output.stats_json} &&
        mv --verbose {params.output_dir}/{params.output_prefix}.fastqc_metrics.csv {output.stats_fastqc} &&
        mv --verbose {params.output_dir}/{params.output_prefix}.trimmer_metrics.csv {output.stats_trimmer} &&
        [ ! -f {params.output_dir}/{params.output_prefix}.mapping_metrics.csv ] || mv --verbose {params.output_dir}/{params.output_prefix}.mapping_metrics.csv {output.stats_map} &&
        mv --verbose {params.output_dir}/{params.output_prefix}.ploidy_estimation_metrics.csv {output.stats_ploidy} &&
        mv --verbose {params.output_dir}/{params.output_prefix}.time_metrics.csv {output.stats_time} &&
        rm -f --verbose {params.output_dir}/*_usage.txt
        ) >{log[0]} 2>&1;
        """


rule shard_count_alns:
    input:
        bam=lambda w: expand(
            "<temp>/<shards>/{tool}/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.bam",
            tool=config["ref"][w.ref]["map"]["tool"],
            allow_missing=True,
        ),
    output:
        counts=temp(
            "<temp>/<shards>/count_alns/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.tsv"
        ),
    log:
        "<logs>/<shards>/count_alns/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.log",
    benchmark:
        "<benchmarks>/<shards>/count_alns/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.jsonl"
    conda:
        f"https://github.com/snakemake/snakemake-wrappers/raw/v7.9.1/bio/samtools/view/environment.yaml"
    threads: 2
    resources:
        mem=lambda w, input, attempt: f"{(0.5* input.size_gb+30)* attempt} GiB",
        runtime=lambda w, input, attempt: f"{(0.03* input.size_gb+0.1)* attempt} h",
    shell:
        """
        (samtools view {input.bam} | awk 'BEGIN{{print "read_id\tn_aligns"}} {{x[$1]++}} END{{for(read_id in x){{print read_id"\t"x[read_id]}}}}') > {output.counts} 2> {log}
        """


rule shard_saturated_reads_filter:
    input:
        counts=expand_pandas(
            rules.shard_count_alns.output.counts, ref_sets, allow_missing=True
        ),
    output:
        read_id=temp(
            "<temp>/<shards>/saturated_reads/filter/{sample}_{library}_{read_type_map}.tsv"
        ),
    log:
        "<logs>/<shards>/saturated_reads/filter/{sample}_{library}_{read_type_map}.log",
    benchmark:
        "<benchmarks>/<shards>/saturated_reads/filter/{sample}_{library}_{read_type_map}.jsonl"
    threads: 2
    resources:
        mem=lambda w, input, attempt: f"{(0.5* input.size_gb+8)* attempt} GiB",
        runtime=lambda w, input, attempt: f"{(0.02* input.size_gb+0.3)* attempt} h",
    params:
        extra="--headerless-tsv-output cat then filter '$n_aligns >= {}' then sort -f read_id then uniq -g read_id then cut -f read_id".format(
            config["filter"]["saturated_reads"]["n_alns"]
        ),
    wrapper:
        "v7.9.1/utils/miller"


rule shard_saturated_reads_remove:
    input:
        bam=rules.shard_count_alns.input.bam,
        read_id=rules.shard_saturated_reads_filter.output.read_id,
    output:
        bam=temp(
            "<temp>/<shards>/saturated_reads/remove/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.bam"
        ),
    log:
        "<logs>/<shards>/saturated_reads/remove/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.log",
    benchmark:
        "<benchmarks>/<shards>/saturated_reads/remove/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.jsonl"
    threads: 2
    resources:
        mem=lambda w, attempt: f"{10* attempt} GiB",
        runtime=lambda w, input, attempt: f"{(0.01* input.size_gb+1)* attempt} h",
    params:
        extra=lambda w, input: f"--qname-file ^{input.read_id}",
    wrapper:
        "v7.9.1/bio/samtools/view"


rule shard_saturated_reads_extract:
    input:
        fastq=get_data,
        read_id=rules.shard_saturated_reads_filter.output.read_id,
    output:
        fq="<results>/<reads>/saturated/{sample}_{library}_{read_type_map}.fastq.gz",
    log:
        "<logs>/<shards>/saturated_reads/extract/{sample}_{library}_{read_type_map}.log",
    benchmark:
        "<benchmarks>/<shards>/saturated_reads/extract/{sample}_{library}_{read_type_map}.jsonl"
    threads: 4
    resources:
        mem=lambda w, attempt: f"{1* attempt} GiB",
        runtime=lambda w, attempt: f"{1* attempt} h",
    params:
        command="subseq",
        compress_lvl=9,
        extra="",
    wrapper:
        "v7.0.0/bio/seqtk"


rule shard_unicorn:
    input:
        bam=(
            rules.shard_saturated_reads_remove.output.bam
            if is_activated("filter/saturated_reads")
            else rules.shard_count_alns.input.bam
        ),
    output:
        bam=temp(
            "<temp>/<shards>/unicorn/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.bam"
        ),
        stats="<stats>/<shards>/unicorn/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.tsv",
    log:
        "<logs>/<shards>/unicorn/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.log",
    benchmark:
        "<benchmarks>/<shards>/unicorn/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.jsonl"
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
        extra="--minrefl 0 --minreads 1",
    shell:
        "unicorn refstats --threads {threads} -b {input.bam} {params.extra} --outbam {output.bam} --outstat {output.stats} >{log} 2>&1"


# https://bioinformatics.stackexchange.com/questions/18538/samtools-sort-most-efficient-memory-and-thread-settings-for-many-samples-on-a-c
rule shard_sort_query:
    input:
        bam=rules.shard_unicorn.output.bam,
    output:
        bam=temp(
            "<temp>/<shards>/sort_query/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.bam"
        ),
    log:
        "<logs>/<shards>/sort_query/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.log",
    benchmark:
        "<benchmarks>/<shards>/sort_query/{sample}_{library}_{read_type_map}.{ref}.{n_shard}-of-{tot_shards}.jsonl"
    threads: 6
    resources:
        mem=lambda w, input, attempt: f"{(10* input.size_gb+20)* attempt} GiB",
        runtime=lambda w, input, attempt: f"{(0.02* input.size_gb+1)* attempt} h",
    params:
        extra="-n",
        mem_overhead_factor=0.2,
    wrapper:
        "v7.9.1/bio/samtools/sort"
