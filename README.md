# NGS Taxon - a generic module for read taxonomic asignment

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.25.0-brightgreen.svg)](https://snakemake.readthedocs.io/en/stable/)

This module implements read derelication steps:
- Map reads agains reference collection
  - [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- Generate MD tag
  - [Samtools CalMD](https://www.htslib.org/doc/samtools-calmd.html)
- Remove Duplicates
  - [Picard MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates)
  - [GATK MarkDuplicatesSpark](https://gatk.broadinstitute.org/hc/en-us/articles/13832682540699-MarkDuplicatesSpark)
  - [Picard MarkDuplicatesWithMateCigar](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicatesWithMateCigar)
  - [BamSorMaDup](https://gitlab.com/german.tischler/biobambam2)
  - [DeDup](https://github.com/apeltzer/DeDup) [TODO]
- Filter BAM references
  - [bam-filter](https://github.com/genomewalker/bam-filter)
- Sort reads by name
  - [Samtools Sort](https://www.htslib.org/doc/samtools-sort.html)
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
