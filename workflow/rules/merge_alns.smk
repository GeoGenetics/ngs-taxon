
refs = pd.DataFrame.from_dict(config["ref"]).transpose()
bowtie_sets = refs.reset_index().rename(columns={"index":"ref", "n_chunks":"tot_chunks"})
bowtie_sets["read_type_map"] = "collapsed"
bowtie_sets["n_chunk"] = [range(1, tot+1) for tot in bowtie_sets["tot_chunks"]]


#############
### RULES ###
#############

rule merge_alns:
    input:
        expand_pandas(rules.sort_name.output.bam, bowtie_sets, allow_missing=True),
    output:
        bam = "results/align/merge_alns/{sample}_{library}_{read_type_map}.bam",
    log:
        "logs/align/merge_alns/{sample}_{library}_{read_type_map}.log"
    benchmark:
        "benchmarks/align/merge_alns/{sample}_{library}_{read_type_map}.tsv"
    params:
        extra = "-n",
    threads: 3
    resources:
        mem = lambda w, attempt: f"{50 * attempt} GB",
        runtime = lambda w, attempt: f"{1 * attempt} d",
        tmpdir = get_tmp(),
    wrapper:
        wrapper_ver + "/bio/samtools/merge"



##########
### QC ###
##########
rule samtools_stats:
    input:
        aln = rules.merge_alns.output.bam,
    output:
        txt = "stats/align/samtools_stats/{sample}_{library}_{read_type_map}.txt",
    log:
        "logs/align/samtools_stats/{sample}_{library}_{read_type_map}.log",
    benchmark:
        "benchmarks/align/samtools_stats/{sample}_{library}_{read_type_map}.log",
    threads: 4
    resources:
        mem = lambda w, attempt: f"{50 * attempt} GB",
        runtime = lambda w, attempt: f"{10 * attempt} h",
    wrapper:
        wrapper_ver + "/bio/samtools/stats"
