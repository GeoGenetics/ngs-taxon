

#############
### RULES ###
#############

rule metadmg_damage:
    input:
        bam = rules.merge_alns.output.bam,
    output:
        dmg = "results/metadmg/damage/{sample}_{library}_{read_type_map}.bdamage.gz",
        res = "results/metadmg/damage/{sample}_{library}_{read_type_map}.res.gz",
#        stats = "stats/metadmg/damage/{sample}_{library}_{read_type_map}.stat.tsv",
        stats = "results/metadmg/damage/{sample}_{library}_{read_type_map}.stat",
    log:
        "logs/metadmg/damage/{sample}_{library}_{read_type_map}.log"
    benchmark:
        "benchmarks/metadmg/damage/{sample}_{library}_{read_type_map}.tsv"
    params:
        out_prefix = lambda w, output: Path(output.dmg.removesuffix(".gz")).with_suffix(""),
        extra = check_cmd(config["metadmg"]["damage"]["params"], forbidden_args = ["-n", "--threads", "-r", "--run_mode", "-o", "--out_prefix"]),
    threads: 4
    resources:
        mem = lambda w, attempt: f"{100 * attempt} GB",
        runtime = lambda w, attempt: f"{2 * attempt} d",
    shell:
        "/projects/caeg/apps/metaDMG-cpp/metaDMG-cpp getdamage --threads {threads} --run_mode 0 {params.extra} --out_prefix {params.out_prefix} {input.bam} > {log} 2>&1"


rule metadmg_lca:
    input:
        bam = rules.merge_alns.output.bam,
        nodes = config["taxonomy"]["nodes"],
        names = config["taxonomy"]["names"],
        acc2tax = config["taxonomy"]["acc2taxid"],
    output:
        dmg = "results/metadmg/lca/{sample}_{library}_{read_type_map}.bdamage.gz",
        lca = "results/metadmg/lca/{sample}_{library}_{read_type_map}.lca.gz",
#        stats = "stats/metadmg/lca/{sample}_{library}_{read_type_map}.stat.gz",
        stats = "results/metadmg/lca/{sample}_{library}_{read_type_map}.stat.gz",
    log:
        "logs/metadmg/lca/{sample}_{library}_{read_type_map}.log"
    benchmark:
        "benchmarks/metadmg/lca/{sample}_{library}_{read_type_map}.tsv"
    params:
        out_prefix = lambda w, output: Path(output.stats.removesuffix(".gz")).with_suffix(""),
        extra = check_cmd(config["metadmg"]["lca"]["params"], forbidden_args = ["--bam", "--nodes", "--names", "--acc2tax", "--temp", "--out", "--out_prefix"]),
    threads: 1
    resources:
        mem = lambda w, attempt: f"{100 * attempt} GB",
        runtime = lambda w, attempt: f"{1 * attempt} d",
    shell:
        "/projects/caeg/apps/metaDMG-cpp/metaDMG-cpp lca --bam {input.bam} --nodes {input.nodes} --names {input.names} --acc2tax {input.acc2tax} {params.extra} --temp {resources.tmpdir}/ --out_prefix {params.out_prefix} > {log} 2>&1"



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
        "benchmarks/metadmg/dfit/{sample}_{library}_{read_type_map}.tsv"
    params:
        out_prefix = lambda w, output: Path(output.dfit.removesuffix(".gz")).with_suffix(""),
        extra = check_cmd(config["metadmg"]["dfit"]["params"], forbidden_args = ["--nthreads", "--node", "--names", "--out", "--out_prefix"]),
    threads: 4
    resources:
        mem_mb = lambda w, attempt, input: 60 * input.size_mb * attempt,
        runtime = lambda w, attempt: f"{10 * attempt} h",
    shell:
        "/projects/caeg/apps/metaDMG-cpp/metaDMG-cpp dfit {input.dmg} --nthreads {threads} --names {input.names} --nodes {input.nodes} {params.extra} --seed $RANDOM --out_prefix {params.out_prefix} > {log} 2>&1"



rule metadmg_aggregate:
    input:
        dmg = rules.metadmg_lca.output.dmg,
        stats = rules.metadmg_lca.output.stats,
        nodes = config["taxonomy"]["nodes"],
        names = config["taxonomy"]["names"],
    output:
        stats = "stats/metadmg/aggregate/{sample}_{library}_{read_type_map}.stat.gz",
    log:
        "logs/metadmg/aggregate/{sample}_{library}_{read_type_map}.log"
    benchmark:
        "benchmarks/metadmg/aggregate/{sample}_{library}_{read_type_map}.tsv"
    params:
        out_prefix = lambda w, output: Path(output.stats.removesuffix(".gz")).with_suffix(""),
    threads: 1
    resources:
        mem = lambda w, attempt: f"{30 * attempt} GB",
        runtime = lambda w, attempt: f"{10 * attempt} h",
    shell:
        "/projects/caeg/apps/metaDMG-cpp/metaDMG-cpp aggregate {input.dmg} --nodes {input.nodes} --names {input.names} --lcastat {input.stats} --out_prefix {params.out_prefix} > {log} 2>&1"
