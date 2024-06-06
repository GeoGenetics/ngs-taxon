
#################
### FUNCTIONS ###
#################


refs = pd.DataFrame.from_dict(config["ref"]).transpose()
bowtie_sets = refs.reset_index().rename(columns={"index":"ref", "n_chunks":"tot_chunks"})
bowtie_sets["read_type_map"] = "collapsed"
bowtie_sets["n_chunk"] = [range(1, tot+1) for tot in bowtie_sets["tot_chunks"]]



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
    rg_info = f"--rg-id '{wildcards.sample}_{wildcards.library}_{flowcell}' --rg 'SM:{wildcards.sample}' --rg 'LB:{wildcards.library}' --rg 'PU:{flowcell}'"
    # Add extra info to RG
    for key, abv in extra_info.items():
        value = units.loc[(wildcards.sample, wildcards.library, slice(None))].get(key)
        if value is not None and value.any():
            rg_info += f" --rg '{abv}:{value.drop_duplicates().item()}'"

    return rg_info


def is_bt2l(wildcards):
    return config["ref"][wildcards.ref].get("bt2l", True)


def get_chunk_aln(wildcards, rule):
    for case in switch(rule):
        if case("merge_alns") or case("cat_alns"):
            if is_activated("align/calmd"):
                return rules.calmd.output.bam
        if case("calmd"):
            if is_activated("align/mark_duplicates"):
                return "temp/align/mark_duplicates/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam"
        if case("mark_duplicates"):
            return rules.sort_coord.output.bam,
        if case():
            raise ValueError(f"Invalid chunk_aln rule specified: {rule}")



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
        bam = temp("temp/align/bowtie2/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam"),
    log:
        "logs/align/bowtie2/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.log"
    benchmark:
        "benchmarks/align/bowtie2/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.jsonl"
    params:
        extra = lambda w: f"--time {get_read_group(w)} " + config["align"]["map"]["params"],
    #priority: lambda w, input: int(input.size_mb)
    threads: 20
    resources:
        mem = lambda w, attempt, input: "{} GiB".format((sum(Path(f).stat().st_size for f in input.idx) / 1024**3 * 0.9 + 100) * attempt),
        runtime = lambda w, attempt: f"{3 * attempt} d",
    wrapper:
        f"{wrapper_ver}/bio/bowtie2/align"


rule clean_header:
    input:
        bam = rules.bowtie2.output.bam,
    output:
        bam = temp("temp/align/clean_header/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam"),
    log:
        "logs/align/clean_header/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.log"
    benchmark:
        "benchmarks/align/clean_header/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.jsonl"
    threads: 4
    resources:
        mem = lambda w, attempt: f"{20 * attempt} GiB",
        runtime = lambda w, attempt: f"{5 * attempt} h",
    shell:
        "/projects/caeg/apps/metaDMG-cpp/misc/compressbam --threads {threads} --input {input.bam} --output {output.bam} >{log} 2>&1"


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
    params:
        mem_overhead_factor=0.2,
    threads: 8
    resources:
        mem = lambda w, attempt, threads: f"{10 * threads * attempt} GiB",
        runtime = lambda w, attempt: f"{1 * attempt} d",
    wrapper:
        f"{wrapper_ver}/bio/samtools/sort"
