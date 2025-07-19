#############
### RULES ###
#############


rule metadmg_damage:
    input:
        unpack(lambda w: get_merge_aln(w, "metadmg_damage")),
    output:
        dmg="results/metadmg/damage/{sample}_{library}_{read_type_map}.bdamage.gz",
        res="results/metadmg/damage/{sample}_{library}_{read_type_map}.res.gz",
        stats="stats/metadmg/damage/{sample}_{library}_{read_type_map}.stat.gz",
        rlen="stats/metadmg/damage/{sample}_{library}_{read_type_map}.rlens.gz",
    log:
        "logs/metadmg/damage/{sample}_{library}_{read_type_map}.log",
    benchmark:
        "benchmarks/metadmg/damage/{sample}_{library}_{read_type_map}.jsonl"
    params:
        out_prefix=lambda w, output: str(
            Path(output.dmg.removesuffix(".gz")).with_suffix("")
        ),
        extra=config["metadmg"]["damage"]["params"],
    conda:
        urlunparse(
            baseurl._replace(path=str(Path(baseurl.path) / "envs" / "metadmg.yaml"))
        )
    threads: 4
    resources:
        mem=lambda w, input, attempt: f"{(0.04* input.size_gb+5)* attempt} GiB",
        runtime=lambda w, input, attempt: f"{(0.03* input.size_gb+0.1)* attempt} h",
    shell:
        """
        metaDMG-cpp getdamage --threads {threads} --run_mode 0 {params.extra} --out_prefix {params.out_prefix} {input.aln} > {log} 2>&1;
        mv {params.out_prefix}.stat.gz {output.stats};
        mv {params.out_prefix}.rlens.gz {output.rlen};
        """


rule metadmg_lca:
    input:
        unpack(lambda w: get_merge_aln(w, "metadmg_lca")),
        nodes=config["taxonomy"]["nodes"],
        names=config["taxonomy"]["names"],
        acc2taxid=[
            config["ref"][ref]["acc2taxid"]
            for ref in config["ref"]
            if config["ref"][ref]["acc2taxid"]
        ],
    output:
        dmg="results/metadmg/lca/{sample}_{library}_{read_type_map}.bdamage.gz",
        lca="results/metadmg/lca/{sample}_{library}_{read_type_map}.lca.gz",
        stats="stats/metadmg/lca/{sample}_{library}_{read_type_map}.stat.gz",
        rlen="stats/metadmg/lca/{sample}_{library}_{read_type_map}.rlens.gz",
    log:
        "logs/metadmg/lca/{sample}_{library}_{read_type_map}.log",
    benchmark:
        "benchmarks/metadmg/lca/{sample}_{library}_{read_type_map}.jsonl"
    params:
        out_prefix=lambda w, output: str(
            Path(output.dmg.removesuffix(".gz")).with_suffix("")
        ),
        extra=config["metadmg"]["lca"]["params"],
    conda:
        urlunparse(
            baseurl._replace(path=str(Path(baseurl.path) / "envs" / "metadmg.yaml"))
        )
    threads: 4
    resources:
        mem=lambda w, input, attempt: f"{(0.1* input.size_gb+13)* attempt} GiB",
        runtime=lambda w, input, attempt: f"{(0.05* input.size_gb+2)* attempt} h",
    shell:
        """
        metaDMG-cpp lca --threads {threads} --bam {input.aln} --nodes {input.nodes} --names {input.names} --acc2tax <(cat {input.acc2taxid}) {params.extra} --temp {resources.tmpdir}/ --reallyDump 1 --out_prefix {params.out_prefix} > {log} 2>&1;
        mv {params.out_prefix}.stat.gz {output.stats};
        mv {params.out_prefix}.rlens.gz {output.rlen};
        """


def _get_library_type(wildcards):
    """Get library type (ss or ds)"""
    return (
        units.loc[(wildcards.sample, wildcards.library, slice(None)), "library_type"]
        .drop_duplicates()
        .item()
    )


rule metadmg_dfit:
    input:
        dmg=rules.metadmg_lca.output.dmg,
        nodes=config["taxonomy"]["nodes"],
        names=config["taxonomy"]["names"],
    output:
        dfit="results/metadmg/dfit/{sample}_{library}_{read_type_map}.dfit.gz",
        boot="results/metadmg/dfit/{sample}_{library}_{read_type_map}.boot.stat.gz",
    log:
        "logs/metadmg/dfit/{sample}_{library}_{read_type_map}.log",
    benchmark:
        "benchmarks/metadmg/dfit/{sample}_{library}_{read_type_map}.jsonl"
    params:
        out_prefix=lambda w, output: str(
            Path(output.dfit.removesuffix(".gz")).with_suffix("")
        ),
        extra=lambda w: f"--doboot 1 --lib {_get_library_type(w)} "
        + config["metadmg"]["dfit"]["params"],
    conda:
        urlunparse(
            baseurl._replace(path=str(Path(baseurl.path) / "envs" / "metadmg.yaml"))
        )
    threads: 4
    resources:
        mem=lambda w, input, attempt: f"{(1.5* Path(input.dmg).stat().st_size/1024**2+10)* attempt} GiB",
        runtime=lambda w, input, attempt: f"{(2* Path(input.dmg).stat().st_size/1024**2+7)* attempt} h",
    shell:
        "metaDMG-cpp dfit {input.dmg} --threads {threads} --names {input.names} --nodes {input.nodes} {params.extra} --seed $RANDOM --out_prefix {params.out_prefix} > {log} 2>&1"


rule metadmg_aggregate:
    input:
        dmg=rules.metadmg_lca.output.dmg,
        lca=rules.metadmg_lca.output.stats,
        dfit=rules.metadmg_dfit.output.dfit,
        nodes=config["taxonomy"]["nodes"],
        names=config["taxonomy"]["names"],
    output:
        stats="results/metadmg/aggregate/{sample}_{library}_{read_type_map}.stat.gz",
    log:
        "logs/metadmg/aggregate/{sample}_{library}_{read_type_map}.log",
    benchmark:
        "benchmarks/metadmg/aggregate/{sample}_{library}_{read_type_map}.jsonl"
    params:
        out_prefix=lambda w, output: str(
            Path(output.stats.removesuffix(".gz")).with_suffix("")
        ),
    conda:
        urlunparse(
            baseurl._replace(path=str(Path(baseurl.path) / "envs" / "metadmg.yaml"))
        )
    threads: 1
    resources:
        mem=lambda w, attempt: f"{25* attempt} GiB",
        runtime=lambda w, attempt: f"{15* attempt} m",
    shell:
        "metaDMG-cpp aggregate {input.dmg} --nodes {input.nodes} --names {input.names} --lcastat {input.lca} --dfit {input.dfit} --out_prefix {params.out_prefix} > {log} 2>&1"
