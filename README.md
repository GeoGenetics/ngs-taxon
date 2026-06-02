# NGS Taxon - a generic module for read taxonomic asignment

[![Snakemake](https://img.shields.io/badge/snakemake-≥9.22.0-brightgreen.svg)](https://snakemake.readthedocs.io/en/stable/)
![CI](https://github.com/GeoGenetics/ngs-taxon/actions/workflows/ci.yml/badge.svg)

This module implements read taxonomic assignment steps:
- Map reads agains reference collection
  - [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) / [BWA_aln](https://github.com/lh3/bwa)
- Remove saturated reads
- Reference stats and clean header (for unused references)
  - [unicorn](https://github.com/GeoGenetics/unicorn)
- Taxonomic assignment
  - [metaDMG](https://github.com/metaDMG-dev)

And QC steps:
- [MultiQC](https://multiqc.info/) (aggregates QC from several of the tools above)
  - [Samtools Stats](https://www.htslib.org/doc/samtools-stats.html)

This module can be used directly, but is designed to be used together with other modules.

## Authors

* Filipe G. Vieira

## Usage

For an example on how to use this module, check repo [aeDNA](https://github.com/GeoGenetics/aeDNA).
