#!/bin/bash

set -euxo pipefail

SNAKEMAKE_OPTS="--snakefile ../../workflow/Snakefile --configfile config/config.yaml --software-deployment-method conda --forceall $@"

for TEST in HD827sonic.all HD827sonic.reassign HD827sonic.filter HD827sonic.none
do
    cd $TEST/
    snakemake $SNAKEMAKE_OPTS --dryrun
    snakemake $SNAKEMAKE_OPTS --rulegraph | dot -Tsvg > rulegraph.svg
    snakemake $SNAKEMAKE_OPTS --filegraph | dot -Tsvg > filegraph.svg
    snakemake $SNAKEMAKE_OPTS --dag | dot -Tsvg > dag.svg

    if [ -d .tests/unit/ ]; then
	pytest -p no:cacheprovider .tests/unit/
    fi
    cd ../
done
