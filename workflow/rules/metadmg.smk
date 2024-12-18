
#############
### RULES ###
#############

rule align_collate:
    input:
        unpack(lambda w: get_merge_aln(w, "collate")),
    output:
        bam = temp("temp/align/collate/{sample}_{library}_{read_type_map}.bam"),
    log:
        "logs/align/collate/{sample}_{library}_{read_type_map}.log"
    benchmark:
        "benchmarks/align/collate/{sample}_{library}_{read_type_map}.jsonl"
    params:
        extra = "-l 6 -r 100000",
    threads: 2
    resources:
        mem = lambda w, attempt: f"{40 * attempt} GiB",
        runtime = lambda w, attempt: f"{2 * attempt} d",
    wrapper:
        f"{wrapper_ver}/bio/samtools/collate"


use rule shard_sort_coord as align_sort_query with:
    input:
        unpack(lambda w: get_merge_aln(w, "sort_query")),
    output:
        bam = temp("temp/align/sort_query/{sample}_{library}_{read_type_map}.bam"),
    log:
        "logs/align/sort_query/{sample}_{library}_{read_type_map}.log"
    benchmark:
        "benchmarks/align/sort_query/{sample}_{library}_{read_type_map}.jsonl"
    params:
        extra = "-n",
        mem_overhead_factor=0.3,


rule metadmg_damage:
    input:
        aln = rules.sort_merged_name.output.bam,
#        unpack(lambda w: get_merge_aln(w, "metadmg_damage")),
    output:
        dmg = "results/metadmg/damage/{sample}_{library}_{read_type_map}.bdamage.gz",
        res = "results/metadmg/damage/{sample}_{library}_{read_type_map}.res.gz",
        stats = "stats/metadmg/damage/{sample}_{library}_{read_type_map}.stat.gz",
        rlen = "stats/metadmg/damage/{sample}_{library}_{read_type_map}.rlens.gz",
    log:
        "logs/metadmg/damage/{sample}_{library}_{read_type_map}.log"
    benchmark:
        "benchmarks/metadmg/damage/{sample}_{library}_{read_type_map}.jsonl"
    params:
        out_prefix = lambda w, output: str(Path(output.dmg.removesuffix(".gz")).with_suffix("")),
        extra = config["metadmg"]["damage"]["params"],
    threads: 4
    resources:
        mem = lambda w, attempt: f"{10 * attempt} GiB",
        runtime = lambda w, attempt: f"{2 * attempt} d",
    shell:
        """
        /projects/caeg/apps/metaDMG-cpp/metaDMG-cpp getdamage --threads {threads} --run_mode 0 {params.extra} --out_prefix {params.out_prefix} {input.aln} > {log} 2>&1;
        gzip {params.out_prefix}.stat; mv {params.out_prefix}.stat.gz {output.stats};
        mv {params.out_prefix}.rlens.gz {output.rlen};
        """


rule metadmg_lca:
    input:
        aln = rules.sort_merged_name.output.bam,
        nodes = config["taxonomy"]["nodes"],
        names = config["taxonomy"]["names"],
        acc2tax = config["taxonomy"]["acc2taxid"],
    output:
        dmg = "results/metadmg/lca/{sample}_{library}_{read_type_map}.bdamage.gz",
        lca = "results/metadmg/lca/{sample}_{library}_{read_type_map}.lca.gz",
        stats = "stats/metadmg/lca/{sample}_{library}_{read_type_map}.stat.gz",
        rlen = "stats/metadmg/lca/{sample}_{library}_{read_type_map}.rlens.gz",
    log:
        "logs/metadmg/lca/{sample}_{library}_{read_type_map}.log"
    benchmark:
        "benchmarks/metadmg/lca/{sample}_{library}_{read_type_map}.jsonl"
    params:
        out_prefix = lambda w, output: str(Path(output.dmg.removesuffix(".gz")).with_suffix("")),
        extra = config["metadmg"]["lca"]["params"],
    threads: 4
    resources:
        mem = lambda w, attempt, input: f"{30 * attempt} GiB",
        runtime = lambda w, attempt, input: f"{1e-4 * input.size_mb * attempt} h",
    shell:
        """
        /projects/caeg/apps/metaDMG-cpp/metaDMG-cpp lca --threads {threads} --bam {input.aln} --nodes {input.nodes} --names {input.names} --acc2tax {input.acc2tax} {params.extra} --temp {resources.tmpdir}/ --out_prefix {params.out_prefix} > {log} 2>&1
        mv {params.out_prefix}.stat.gz {output.stats};
        mv {params.out_prefix}.rlens.gz {output.rlen};
        """



def _get_library_type(wildcards):
    """ Get library typ (ss or ds) """
    return units.loc[(wildcards.sample, wildcards.library, slice(None)), "library_type"].drop_duplicates().item()

rule metadmg_dfit:
    input:
        dmg = rules.metadmg_lca.output.dmg,
        nodes = config["taxonomy"]["nodes"],
        names = config["taxonomy"]["names"],
    output:
        dfit = "results/metadmg/dfit/{sample}_{library}_{read_type_map}.dfit.gz",
        boot = "results/metadmg/dfit/{sample}_{library}_{read_type_map}.boot.stat.gz",
    log:
        "logs/metadmg/dfit/{sample}_{library}_{read_type_map}.log"
    benchmark:
        "benchmarks/metadmg/dfit/{sample}_{library}_{read_type_map}.jsonl"
    params:
        out_prefix = lambda w, output: str(Path(output.dfit.removesuffix(".gz")).with_suffix("")),
        extra = lambda w: f"--doboot 1 --lib {_get_library_type(w)} " + config["metadmg"]["dfit"]["params"],
    threads: 4
    resources:
        mem = lambda w, attempt: f"{25 * attempt} GiB",
        runtime = lambda w, attempt, input: f"{7 * attempt} h",
    shell:
        "/projects/caeg/apps/metaDMG-cpp/metaDMG-cpp dfit {input.dmg} --threads {threads} --names {input.names} --nodes {input.nodes} {params.extra} --seed $RANDOM --out_prefix {params.out_prefix} > {log} 2>&1"



rule metadmg_aggregate:
    input:
        dmg = rules.metadmg_lca.output.dmg,
        lca = rules.metadmg_lca.output.stats,
        dfit = rules.metadmg_dfit.output.dfit,
        nodes = config["taxonomy"]["nodes"],
        names = config["taxonomy"]["names"],
    output:
        stats = "stats/metadmg/aggregate/{sample}_{library}_{read_type_map}.stat.gz",
    log:
        "logs/metadmg/aggregate/{sample}_{library}_{read_type_map}.log"
    benchmark:
        "benchmarks/metadmg/aggregate/{sample}_{library}_{read_type_map}.jsonl"
    params:
        out_prefix = lambda w, output: str(Path(output.stats.removesuffix(".gz")).with_suffix("")),
    threads: 1
    resources:
        mem = lambda w, attempt: f"{25 * attempt} GiB",
        runtime = lambda w, attempt: f"{15 * attempt} m",
    shell:
        "/projects/caeg/apps/metaDMG-cpp/metaDMG-cpp aggregate {input.dmg} --nodes {input.nodes} --names {input.names} --lcastat {input.lca} --dfit {input.dfit} --out_prefix {params.out_prefix} > {log} 2>&1"
