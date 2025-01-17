
#################
### FUNCTIONS ###
#################


refs = pd.DataFrame.from_dict(config["ref"]).transpose()
ref_sets = refs.reset_index().rename(columns={"index":"ref", "n_chunks":"tot_chunks"})
ref_sets["read_type_map"] = "collapsed"
ref_sets["n_chunk"] = [range(1, tot+1) for tot in ref_sets["tot_chunks"]]
ref_sets.loc[ref_sets["bt2l"].isna(), "bt2l"] = True

wildcard_constraints:
  read_type_map="|".join(["pe", "se", "singleton", "collapsed"])


# https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups
def get_read_group(wildcards):
    """Denote sample name and platform in read group."""

    extra_info = {"center": "CN",
                  "date": "DT",
                  "platform": "PL",
                  "platform_model": "PM",
                 }

    # Get flowcell
    flowcell = units.loc[(wildcards.sample, wildcards.library, slice(None)), "flowcell"].drop_duplicates().item()
    rg_info = [f"{wildcards.sample}_{wildcards.library}_{flowcell}_{wildcards.ref}{wildcards.n_chunk}o{wildcards.tot_chunks}", f"SM:{wildcards.sample}", f"LB:{wildcards.library}", f"PU:{flowcell}"]
    # Add extra info to RG
    for key, abv in extra_info.items():
        value = units.loc[(wildcards.sample, wildcards.library, slice(None))].get(key)
        if value is not None and value.any():
            rg_info.append(f"{abv}:{value.drop_duplicates().item()}")

    return rg_info


def is_bt2l(wildcards):
    return config["ref"][wildcards.ref].get("bt2l", True)



#############
### RULES ###
#############

def get_data(wildcards):
    return expand("results/reads/low_complexity/{sample}_{library}_{read_type_trim}.fastq.gz", read_type_trim = get_read_type_trim(wildcards.read_type_map), allow_missing = True)


def get_index(wildcards):
    if is_bt2l(wildcards):
        return multiext(config["ref"][wildcards.ref]["path"], ".1.bt2l", ".2.bt2l", ".3.bt2l", ".4.bt2l", ".rev.1.bt2l", ".rev.2.bt2l")
    else:
        return multiext(config["ref"][wildcards.ref]["path"], ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")

rule bowtie2:
    input:
        sample = get_data,
        idx = get_index,
    output:
        bam = temp("temp/chunk/bowtie2/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam"),
    log:
        "logs/chunk/bowtie2/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.log"
    benchmark:
        "benchmarks/chunk/bowtie2/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.jsonl"
    params:
        extra = lambda w: f"""--time --rg-id '{"' --rg '".join(get_read_group(w))}' --rg 'PG:bowtie2' """ + config["align"]["map"]["params"],
    threads: 20
    resources:
        mem = lambda w, attempt, input: "{} GiB".format((0.7 * sum(Path(f).stat().st_size for f in input.idx) / 1024**3 + 50) * attempt),
        runtime = lambda w, attempt: f"{3 * attempt} d",
        slurm_extra = "--nice=2000",
    wrapper:
        f"{wrapper_ver}/bio/bowtie2/align"


rule bwa_aln:
    input:
        fastq = get_data,
        idx = lambda w: multiext(config["ref"][w.ref]["path"], ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        sai = temp("temp/chunk/bwa_aln/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.sai")
    log:
        "logs/chunk/bwa_aln/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.log"
    benchmark:
        "benchmarks/chunk/bwa_aln/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.jsonl"
    params:
        extra = lambda w: f"""-r '@RG\\tID:{"\\t".join(get_read_group(w))}\\tPG:bwa_aln' """ + config["align"]["map"]["params"],
    threads: 10
    resources:
        mem = lambda w, attempt, input: "{} GiB".format((2 * sum(Path(f).stat().st_size for f in input.idx) / 1024**3 + 50) * attempt),
        runtime = lambda w, attempt: f"{1 * attempt} d",
        slurm_extra = "--nice=2000",
    wrapper:
        f"{wrapper_ver}/bio/bwa/aln"


rule bwa_samxe:
    input:
        fastq = get_data,
        sai = rules.bwa_aln.output.sai,
        idx = rules.bwa_aln.input.idx,
    output:
        bam = temp("temp/chunk/bwa_aln/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam"),
    log:
        "logs/chunk/bwa_aln/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.log"
    benchmark:
        "benchmarks/chunk/bwa_aln/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.jsonl"
    params:
        extra = lambda w: f"""-r '@RG\tID:{"\t".join(get_read_group(w))}' """ + config["align"]["map"]["params"],
        sort = "samtools",
    threads: 1
    resources:
        mem = lambda w, attempt, input: "{} GiB".format((2 * sum(Path(f).stat().st_size for f in input.idx) / 1024**3 + 50) * attempt),
        runtime = lambda w, attempt: f"{1 * attempt} d",
    wrapper:
        f"{wrapper_ver}/bio/bwa/samxe"



rule chunk_clean_header:
    input:
        bam = expand("temp/chunk/{tool}/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam", tool=config["align"]["map"]["tool"], allow_missing=True),
    output:
        bam = temp("temp/chunk/clean_header/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam"),
    log:
        "logs/chunk/clean_header/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.log"
    benchmark:
        "benchmarks/chunk/clean_header/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.jsonl"
    threads: 1
    resources:
        mem = lambda w, attempt: f"{20 * attempt} GiB",
        runtime = lambda w, attempt: f"{2 * attempt} h",
    shell:
        "/projects/caeg/apps/metaDMG-bb/misc/reheadbam -b {input.bam} -o {output.bam} >{log} 2>&1"


# https://bioinformatics.stackexchange.com/questions/18538/samtools-sort-most-efficient-memory-and-thread-settings-for-many-samples-on-a-c
rule chunk_sort_query:
    input:
        bam = rules.chunk_clean_header.output.bam,
    output:
        bam = temp("temp/chunk/sort_query/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam"),
    log:
        "logs/chunk/sort_query/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.log"
    benchmark:
        "benchmarks/chunk/sort_query/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.jsonl"
    params:
        extra="-N",
        mem_overhead_factor=0.1,
    threads: 8
    resources:
        mem = lambda w, attempt, threads, input: f"{10 * threads * attempt} GiB",
        runtime = lambda w, attempt, input: f"{np.clip(0.0001 * input.size_mb + 1, 0.1, 20) * attempt} h",
    wrapper:
        f"{wrapper_ver}/bio/samtools/sort"
