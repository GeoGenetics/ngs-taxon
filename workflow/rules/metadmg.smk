

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
        out_prefix = lambda w, output: Path(output.stats).with_suffix(""),
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
#        stats = "stats/metadmg/lca/{sample}_{library}_{read_type_map}.stat.tsv",
        stats = "results/metadmg/lca/{sample}_{library}_{read_type_map}.stat",
    log:
        "logs/metadmg/lca/{sample}_{library}_{read_type_map}.log"
    benchmark:
        "benchmarks/metadmg/lca/{sample}_{library}_{read_type_map}.tsv"
    params:
        out_prefix = lambda w, output: Path(output.stats).with_suffix(""),
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
        dfit = "results/metadmg/dfit/{sample}_{library}_{read_type_map}.dfit.txt.gz",
#        stats_boot = "stats/metadmg/dfit/{sample}_{library}_{read_type_map}.boot.stats.tsv.gz",
        stats_boot = "results/metadmg/dfit/{sample}_{library}_{read_type_map}.boot.stat.txt.gz",
    log:
        "logs/metadmg/dfit/{sample}_{library}_{read_type_map}.log"
    benchmark:
        "benchmarks/metadmg/dfit/{sample}_{library}_{read_type_map}.tsv"
    params:
        out_prefix = lambda w, output: Path(output.dfit).with_suffix("").with_suffix("").with_suffix(""),
        extra = check_cmd(config["metadmg"]["dfit"]["params"], forbidden_args = ["--nthreads", "--node", "--names", "--out", "--out_prefix"]),
    threads: 10
    resources:
        mem_mb = lambda w, attempt, input: 40 * input.size_mb * attempt,
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
#        stats = "stats/metadmg/aggregate/{sample}_{library}_{read_type_map}.stats.tsv.gz",
        stats = "stats/metadmg/aggregate/{sample}_{library}_{read_type_map}.aggregate.stat.txt.gz",
    log:
        "logs/metadmg/aggregate/{sample}_{library}_{read_type_map}.log"
    benchmark:
        "benchmarks/metadmg/aggregate/{sample}_{library}_{read_type_map}.tsv"
    params:
        out_prefix = lambda w, output: Path(output.stats).with_suffix("").with_suffix("").with_suffix("").with_suffix(""),
    threads: 1
    resources:
        mem = lambda w, attempt: f"{30 * attempt} GB",
        runtime = lambda w, attempt: f"{10 * attempt} h",
    shell:
        "/projects/caeg/apps/metaDMG-cpp/metaDMG-cpp aggregate {input.dmg} --nodes {input.nodes} --names {input.names} --lcastat {input.stats} --out_prefix {params.out_prefix} > {log} 2>&1"
