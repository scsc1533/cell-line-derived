#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Extract strand-specific FASTA from GENCODE v43 BED files and
# write FIMO-friendly nameOnly headers in one step.
#
# Input BED name column:
#   ENSG00000119661|DNAL1|three_prime_UTR
#
# bedtools getfasta -name usually writes headers like:
#   >ENSG00000119661|DNAL1|three_prime_UTR::chr14:...(+)
#
# FIMO may parse the ::chr:start-end suffix as genomic coordinates and then
# output sequence_name as only the chromosome. Therefore this script removes
# the suffix and writes headers like:
#   >ENSG00000119661|DNAL1|three_prime_UTR
# ============================================================

bed_dir="./03_GENCODE_v43_UTR_CDS_BED_EMDtested_Hep3B_background"
outdir="./04_GENCODE_v43_UTR_CDS_FASTA_EMDtested_Hep3B_background_FIMO_nameOnly"
genome_fa="./GRCh38_GENCODE_20231021/fa/GRCh38.p13.genome.fa"
bedtools="/opt/software/miniforge/envs/rnaseq/bin/bedtools"

mkdir -p "${outdir}"

if [[ ! -s "${genome_fa}" ]]; then
  echo "ERROR: genome FASTA not found or empty: ${genome_fa}" >&2
  exit 1
fi

if [[ ! -x "${bedtools}" ]]; then
  echo "ERROR: bedtools not found or not executable: ${bedtools}" >&2
  exit 1
fi

if [[ ! -d "${bed_dir}" ]]; then
  echo "ERROR: BED directory not found: ${bed_dir}" >&2
  exit 1
fi

summary_file="${outdir}/extract_region_FASTA_summary.tsv"
chrom_check_file="${outdir}/BED_chromosome_not_in_genome.tsv"

printf "bed_file\tfasta_file\tbed_interval_n\tfasta_record_n\texample_header\n" > "${summary_file}"
printf "bed_file\tchrom\n" > "${chrom_check_file}"

genome_fai="${genome_fa}.fai"
if [[ ! -s "${genome_fai}" ]]; then
  "${bedtools}" faidx "${genome_fa}"
fi

tmp_genome_chrom="$(mktemp)"
trap 'rm -f "${tmp_genome_chrom}"' EXIT

cut -f1 "${genome_fai}" | sort -u > "${tmp_genome_chrom}"

shopt -s nullglob
bed_files=("${bed_dir}"/*.bed)
shopt -u nullglob

if [[ "${#bed_files[@]}" -eq 0 ]]; then
  echo "ERROR: no BED files found in ${bed_dir}" >&2
  exit 1
fi

for bed in "${bed_files[@]}"; do
  base="$(basename "${bed}" .bed)"
  fa="${outdir}/${base}.fa"
  tmp_fa="${outdir}/${base}.with_coordinates.tmp.fa"

  if [[ ! -s "${bed}" ]]; then
    echo "Skip empty BED: ${bed}" >&2
    printf "%s\t%s\t0\t0\t%s\n" "$(basename "${bed}")" "$(basename "${fa}")" "" >> "${summary_file}"
    : > "${fa}"
    continue
  fi

  tmp_bed_chrom="$(mktemp)"
  cut -f1 "${bed}" | sort -u > "${tmp_bed_chrom}"
  comm -23 "${tmp_bed_chrom}" "${tmp_genome_chrom}" | awk -v b="$(basename "${bed}")" 'BEGIN{OFS="\t"} {print b, $1}' >> "${chrom_check_file}"
  rm -f "${tmp_bed_chrom}"

  "${bedtools}" getfasta \
    -fi "${genome_fa}" \
    -bed "${bed}" \
    -s \
    -name \
    -fo "${tmp_fa}"

  awk '
    /^>/ {
      h = $0
      sub(/^>[[:space:]]*/, ">", h)
      sub(/::.*/, "", h)
      print h
      next
    }
    {
      print
    }
  ' "${tmp_fa}" > "${fa}"

  rm -f "${tmp_fa}"

  bed_n="$(wc -l < "${bed}" | tr -d ' ')"
  fa_n="$(grep -c '^>' "${fa}" || true)"
  example_header="$(grep -m 1 '^>' "${fa}" || true)"

  printf "%s\t%s\t%s\t%s\t%s\n" \
    "$(basename "${bed}")" "$(basename "${fa}")" "${bed_n}" "${fa_n}" "${example_header}" >> "${summary_file}"
done

echo "Done."
echo "Output directory: ${outdir}"
echo "Summary: ${summary_file}"
echo "Chromosome check: ${chrom_check_file}"