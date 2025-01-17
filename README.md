# NGS Taxon - a generic module for read taxonomic asignment

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.11.2-brightgreen.svg)](https://snakemake.readthedocs.io/en/stable/)

This module implements read derelication steps:
- Map reads agains reference collection
  - [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) / [BWA_aln](https://github.com/lh3/bwa)
- Clean header for unused references
- Reassign reads
  - [bam-filter](https://github.com/genomewalker/bam-filter)
- Filter BAM references
  - [bam-filter](https://github.com/genomewalker/bam-filter)
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
