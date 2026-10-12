#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: 02_subset_popgen_vcf.sh INPUT.vcf.gz ASSIGNED_RAW_IDS.txt OUTPUT_DIR" >&2
  exit 2
fi

input_vcf=$1
sample_file=$2
output_dir=$3
mkdir -p "$output_dir"

for command_name in bcftools python3; do
  command -v "$command_name" >/dev/null 2>&1 || {
    echo "Required command not found: $command_name" >&2
    exit 1
  }
done

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
all_vcf="$output_dir/s215_popgen_all_snps.vcf.gz"
one_vcf="$output_dir/s215_popgen_one_snp_per_locus.vcf.gz"
regions="$output_dir/one_snp_per_ipyrad_locus.regions.tsv"
selection_report="$output_dir/one_snp_per_ipyrad_locus.tsv"

echo "Subsetting assigned samples and removing sites monomorphic after subsetting..."
bcftools view \
  --samples-file "$sample_file" \
  --min-alleles 2 --max-alleles 2 --types snps \
  --min-ac 1:minor \
  --output-type z --output-file "$all_vcf" \
  "$input_vcf"
bcftools index --tbi --force "$all_vcf"

echo "Selecting the highest-call-rate SNP in every ipyrad locus..."
python3 "$script_dir/02_select_one_snp_per_locus.py" \
  "$all_vcf" "$regions" "$selection_report"

bcftools view \
  --regions-file "$regions" \
  --output-type z --output-file "$one_vcf" \
  "$all_vcf"
bcftools index --tbi --force "$one_vcf"

all_samples=$(bcftools query -l "$all_vcf" | wc -l | tr -d ' ')
all_snps=$(bcftools index --nrecords "$all_vcf")
one_snps=$(bcftools index --nrecords "$one_vcf")
{
  echo -e "dataset\tsamples\tsnps"
  echo -e "all_snps\t${all_samples}\t${all_snps}"
  echo -e "one_snp_per_ipyrad_locus\t${all_samples}\t${one_snps}"
} > "$output_dir/vcf_summary.tsv"

echo "Wrote $all_vcf"
echo "Wrote $one_vcf"

