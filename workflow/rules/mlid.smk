

rule mlid_query:
    """
    Creates a query file from a single BAM.

    Alignments to accessions not found in the acc2taxid will be dropped. Header lengths will be
    subtly "wrong" in the output file, and get patched up in `annotate` below.
    """
    input:
        bam=rules.align_merge.output.bam,
        acc2taxid=[
            config["ref"][ref]["acc2taxid"]
            for ref in config["ref"]
            if config["ref"][ref].get("acc2taxid", False)
        ],
    output:
        qry=pipe("<temp>/<mlid>/query/{sample}_{library}_{read_type_map}.bry"),
    log:
        "<logs>/<mlid>/query/{sample}_{library}_{read_type_map}.log",
    benchmark:
        "<benchmarks>/<mlid>/query/{sample}_{library}_{read_type_map}.jsonl"
    threads: 4  # Only used for BAM decompression.
    resources:
        mem=lambda w, attempt: f"{20* attempt} GiB",
        runtime=lambda w, attempt: f"{1* attempt} h",
    params:
        remap=lambda wildcards, input: " ".join(
            f"--remap-file {p}" for p in input.acc2taxid
        ),
    shell:
        "mlid query --threads {threads} {params.remap} --remap-unmapped-drop --output {output.qry} {input.bam} 2> {log}"


rule mlid_annotate:
    """
    Fixes genome lengths in the query file header.

    By default, `mlid query` sums the accession lengths for a taxid when remapping. However,
    production pipeline drops unused accessions (even for used taxids), meaning that the resulting
    summed length is often wrong. This corrects the header by annotating the reference header
    entries with the corrected lengths.
    """
    input:
        qry=rules.mlid_query.output.qry,
        lengths=[
            config["ref"][ref]["lengths"]
            for ref in config["ref"]
            if config["ref"][ref].get("lengths", False)
        ],
    output:
        qry=pipe("<temp>/<mlid>/annotate/{sample}_{library}_{read_type_map}.bry"),
    log:
        "<logs>/<mlid>/annotate/{sample}_{library}_{read_type_map}.log",
    benchmark:
        "<benchmarks>/<mlid>/annotate/{sample}_{library}_{read_type_map}.jsonl"
    threads: 1
    resources:
        mem=lambda w, attempt: f"{20* attempt} GiB",
        runtime=lambda w, attempt: f"{1* attempt} h",
    params:
        # Length file columns: first column has the name (i.e. taxid), second the corrected length.
        format="name,length",
        lengths=lambda wildcards, input: " ".join(
            f"--header-annotations {p}" for p in input.lengths
        ),
    shell:
        "mlid annotate --header-columns reference={params.format} {params.lengths} --output {output.qry} {input.qry} 2> {log}"


rule mlid_abundance:
    """
    Estimates abundances and posterior probabilities.
    """
    input:
        qry=rules.mlid_annotate.output.qry,
    output:
        qry=pipe("<temp>/<mlid>/abundance/{sample}_{library}_{read_type_map}.bry"),
    log:
        "<logs>/<mlid>/abundance/{sample}_{library}_{read_type_map}.log",
    benchmark:
        "<benchmarks>/<mlid>/abundance/{sample}_{library}_{read_type_map}.jsonl"
    threads: 10
    resources:
        mem=lambda w, attempt: f"{20* attempt} GiB",
        runtime=lambda w, attempt: f"{1* attempt} h",
    params:
        extra=config["mlid"]["abundance"]["params"],
    shell:
        "mlid abundance --threads {threads} {params.extra} --temp-dir {resources.tmpdir} --output {output.qry} {input.qry} 2> {log}"


rule mlid_lca:
    """
    Adds posterior LCA assignments based on posterior probabilities.

    A simple LCA is also output for comparison.
    """
    input:
        qry=rules.mlid_abundance.output.qry,
        taxdump=config["taxdump"],
    output:
        qry="<results>/<mlid>/lca/{sample}_{library}_{read_type_map}.bry.gz",
    log:
        "<logs>/<mlid>/lca/{sample}_{library}_{read_type_map}.log",
    benchmark:
        "<benchmarks>/<mlid>/lca/{sample}_{library}_{read_type_map}.jsonl"
    threads: 1
    resources:
        mem=lambda w, attempt: f"{20* attempt} GiB",
        runtime=lambda w, attempt: f"{1* attempt} h",
    params:
        extra=config["mlid"]["lca"]["params"],
    shell:
        "mlid lca --taxdump {input.taxdump} {params.extra} --output {output.qry} {input.qry} 2> {log}"


rule mlid_extract:
    """Extract global abundances to a simple TSV file for convenience."""
    input:
        qry=rules.mlid_lca.output.qry,
    output:
        tsv="<results>/<mlid>/extract/{sample}_{library}_{read_type_map}.tsv.gz",
    log:
        "<logs>/<mlid>/extract/{sample}_{library}_{read_type_map}.log",
    benchmark:
        "<benchmarks>/<mlid>/extract/{sample}_{library}_{read_type_map}.jsonl"
    threads: 1
    resources:
        mem=lambda w, attempt: f"{20* attempt} GiB",
        runtime=lambda w, attempt: f"{1* attempt} h",
    params:
        header_format=lambda _: "#%{name}\t%{abundance}",
        format=lambda _: "%{name}\t%{info/ALCA}\t[%{ref}|%{fmt/PP}]{,}",
    shell:
        """
        mlid view --elect --output-format b {input.qry} | mlid extract --header-format 'reference={params.header_format}' --format '{params.format}' --print-header | gzip > {output.tsv} 2>{log}
        """
