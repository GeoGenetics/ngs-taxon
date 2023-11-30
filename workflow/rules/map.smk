
#################
### FUNCTIONS ###
#################

# Add MAP read type
units["read_type_map"] = [get_read_type_map(u.sample, u.library, u.lane) for u in units.itertuples()]


# https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups
def get_read_group(wildcards):
    """Denote sample name and platform in read group."""

    flowcell = pd.unique(units.loc[(wildcards.sample, wildcards.library), "flowcell"])
    if flowcell.size != 1:
        raise ValueError("no flowcell provided")

    platform = pd.unique(units.loc[(wildcards.sample, wildcards.library), "platform"])
    if platform.size != 1:
        raise ValueError("no platform provided")

    return f"--rg-id '{wildcards.sample}' --rg 'PU:{flowcell[0]}.{wildcards.library}' --rg 'SM:{wildcards.sample}' --rg 'LB:{wildcards.library}' --rg 'PL:{platform[0]}'"


def is_bt2l(wildcards):
    return config["ref"][wildcards.ref].get("bt2l", True)


def get_chunk_aln(wildcards, rule, ext=["bam"]):
    for case in switch(rule):
        if case("sort_name"):
            if is_activated("align/calmd"):
                return expand_ext(rules.calmd.output.bam, ext)
        if case("calmd"):
            if is_activated("align/mark_duplicates"):
                return expand_ext(format_partial("temp/align/{tool}/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam", tool=config["align"]["mark_duplicates"]["tool"]), ext)
        if case("mark_duplicates"):
            if is_activated("align/bam_filter"):
                return expand_ext(rules.bam_filter.output.bam, ext)
        if case("bam_filter"):
            if rule != "sort_name":
                return expand_ext(rules.sort_coord.output.bam, ext)
            else:
                return expand_ext(rules.sort_coord.input.bam, ext)
        if case():
            raise ValueError("Invalid ALN rule specified!")



#############
### RULES ###
#############

wildcard_constraints:
    read_type_map  = "|".join(set(flatten(units.read_type_map))),


def get_data(wildcards):
    if is_activated("reads/extension"):
        reads_file = "temp/reads/represent/grep/{tool}/{sample}_{library}_{read_type_trim}.fastq.gz"
    elif is_activated("reads/derep"):
        reads_file = "temp/reads/derep/{tool}/{sample}_{library}_{read_type_trim}.fastq.gz"
    else:
        reads_file = "temp/reads/derep/merge_lanes/{sample}_{library}_{read_type_trim}.fastq.gz"

    return expand(reads_file, tool=config["reads"]["derep"]["tool"], read_type_trim = get_read_type_trim(wildcards.read_type_map), allow_missing = True)


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
        "benchmarks/align/bowtie2/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.tsv"
    params:
        extra = lambda w: "--time {} ".format(get_read_group(w)) + check_cmd(config["align"]["map"]["params"], forbidden_args = ["--threads", "--mm", "-t", "--time", "-q", "-U", "-1", "-2", "--interleaved", "-x", "-b", "-S", "-rd-id", "-rg"]),
#    priority: lambda w, input: int(input.size_mb)
    threads: 20
    resources:
        mem = lambda w, attempt, input: f"{(sum(f.size for f in input.idx) / 1024**2 * 0.9 + 105000) * attempt} MB",
        runtime = lambda w, attempt: f"{3 * attempt} d",
    wrapper:
        wrapper_ver + "/bio/bowtie2/align"


rule clean_header:
    input:
        bam = rules.bowtie2.output.bam,
    output:
        bam = temp("temp/align/clean_header/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.bam"),
    log:
        "logs/align/clean_header/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.log"
    benchmark:
        "benchmarks/align/clean_header/{sample}_{library}_{read_type_map}.{ref}.{n_chunk}-of-{tot_chunks}.tsv"
    threads: 4
    resources:
        mem = lambda w, attempt: f"{50 * attempt} GB",
        runtime = lambda w, attempt: f"{1 * attempt} d",
    shell:
        "/projects/caeg/apps/metaDMG-cpp/misc/compressbam --threads {threads} --input {input.bam} --output {output.bam} >{log} 2>&1"
