#!/usr/bin/env bash
set -euo pipefail

INPUT_VCF="${1:-data/s384_batchclean_ab25_locusFDR05_mac4_dp15-150_miss10.vcf.gz}"
OUTPUT_VCF="${2:-data/s384_batchclean_ab25_locusFDR05_mac4_dp15-150_miss10.sorted.vcf.gz}"

WORK_DIR="results/vcf_header_repair"
mkdir -p "$WORK_DIR"

CONTIG_HEADER="$WORK_DIR/s384_contigs.hdr"
BASE_HEADER="$WORK_DIR/s384_header_without_contigs.txt"
REPAIRED_HEADER="$WORK_DIR/s384_repaired_header.txt"
REHEADERED_VCF="$WORK_DIR/s384_batchclean.reheadered.vcf.gz"

command -v bcftools >/dev/null ||
  { echo "ERROR: bcftools was not found"; exit 1; }

[[ -f "$INPUT_VCF" ]] ||
  { echo "ERROR: input VCF not found: $INPUT_VCF"; exit 1; }

echo "Checking compressed VCF..."
gzip -t "$INPUT_VCF"

echo "Extracting reference contigs in first-appearance order..."
gzip -dc "$INPUT_VCF" |
  awk -F '\t' '
    !/^#/ && !seen[$1]++ {
      print "##contig=<ID=" $1 ">"
    }
  ' > "$CONTIG_HEADER"

N_CONTIGS=$(wc -l < "$CONTIG_HEADER" | tr -d ' ')

if [[ "$N_CONTIGS" -eq 0 ]]; then
  echo "ERROR: no contigs were recovered from the VCF"
  exit 1
fi

echo "Recovered $N_CONTIGS reference contigs"

echo "Extracting the existing VCF header..."
#gzip -dc "$INPUT_VCF" |
#  awk '
#    /^##contig=/ {next}
#    /^#/ {print; next}
#    {exit}
#  ' > "$BASE_HEADER"

gzip -dc "$INPUT_VCF" |
  awk '
    /^##contig=/ {next}
    /^#/ {print}
  ' > "$BASE_HEADER"


echo "Adding contig declarations before the #CHROM line..."
awk '
  FNR == NR {
    contig[++n] = $0
    next
  }
  /^#CHROM/ {
    for (i = 1; i <= n; i++) print contig[i]
  }
  {
    print
  }
' "$CONTIG_HEADER" "$BASE_HEADER" > "$REPAIRED_HEADER"

echo "Applying the repaired header..."
bcftools reheader \
  -h "$REPAIRED_HEADER" \
  -o "$REHEADERED_VCF" \
  "$INPUT_VCF"

echo "Sorting variants by contig and position..."
bcftools sort \
  "$REHEADERED_VCF" \
  -Oz \
  -o "$OUTPUT_VCF"

echo "Creating Tabix index..."
bcftools index -f -t "$OUTPUT_VCF"

echo "Validating output..."
gzip -t "$OUTPUT_VCF"

N_SAMPLES=$(bcftools query -l "$OUTPUT_VCF" | wc -l | tr -d ' ')
N_VARIANTS=$(bcftools view -H "$OUTPUT_VCF" | wc -l | tr -d ' ')
N_HEADER_CONTIGS=$(bcftools view -h "$OUTPUT_VCF" |
  grep -c '^##contig=')

echo
echo "Finished:"
echo "  Output:   $OUTPUT_VCF"
echo "  Index:    ${OUTPUT_VCF}.tbi"
echo "  Samples:  $N_SAMPLES"
echo "  Variants: $N_VARIANTS"
echo "  Contigs:  $N_HEADER_CONTIGS"
