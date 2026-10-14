#!/usr/bin/env bash
set -euo pipefail

# Run after the main population workflow and sequence-pi bootstrap are complete.
# The revised spID assignments must first be propagated by rerunning script 01.
nboot=${NBOOT:-1000}

NBOOT="$nboot" Rscript scripts/10_within_lineage_private_alleles.R \
  results/vcf/s215_popgen_one_snp_per_locus.vcf.gz \
  results/metadata/sample_to_population.tsv \
  results/refugia/within_lineage

Rscript scripts/11_summarize_refugia_metrics.R \
  results/diversity \
  results/refugia/within_lineage \
  results/refugia

echo "Refugial-diversity tables are in results/refugia"

