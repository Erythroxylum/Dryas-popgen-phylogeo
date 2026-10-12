#!/usr/bin/env bash
set -euo pipefail

# Run from the repository root after activating the dryas-popgen environment.
metadata=${METADATA:-data/Dryas_sampledata_popgen.csv}
vcf=${VCF:-data/s384_batchclean_ab25_locusFDR05_mac4_dp15-150_miss10.sorted.vcf.gz}
nboot=${NBOOT:-1000}

Rscript scripts/01_validate_metadata.R "$metadata" "$vcf" results/metadata

bash scripts/02_subset_popgen_vcf.sh \
  "$vcf" \
  results/metadata/assigned_raw_ids.txt \
  results/vcf

NBOOT="$nboot" Rscript scripts/03_population_diversity.R \
  results/vcf/s215_popgen_one_snp_per_locus.vcf.gz \
  results/metadata/sample_to_population.tsv \
  results/diversity

Rscript scripts/04_individual_heterozygosity.R \
  results/vcf/s215_popgen_all_snps.vcf.gz \
  results/metadata/sample_to_population.tsv \
  results/diversity/individual_snp_heterozygosity.tsv

NBOOT="$nboot" Rscript scripts/05_pairwise_fst.R \
  results/vcf/s215_popgen_one_snp_per_locus.vcf.gz \
  results/metadata/sample_to_population.tsv \
  results/fst

echo "Core population-genetic analyses complete."

