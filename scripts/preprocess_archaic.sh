#!/usr/bin/env bash
# One-off preprocessing before src/08_alle_chromosomen.py:
# archaic BCFs (Skov et al., Zenodo 13368126) -> one TSV per chromosome
# (POS, REF, ALT, genotype of the 4 archaic samples, in the sample order of the BCF).
#
# NOTE: reconstructed from the project documentation (Appendix A). The command as originally typed
# was not kept. The genotype column order is hard-wired in src/genom_chrom.py; check it with:
#   bcftools query -l individuals_highcov.12.bcf
#
# Usage: GENOM_BASIS=/path/to/data scripts/preprocess_archaic.sh   (default: ./data)
set -euo pipefail

BASIS="${GENOM_BASIS:-data}"
IN="$BASIS/Datensätze/Neandertaler_hg38"
OUT="$BASIS/Zwischen/neandertaler"
mkdir -p "$OUT"

for c in $(seq 1 22); do
  bcftools query -f '%POS\t%REF\t%ALT[\t%GT]\n' "$IN/individuals_highcov.$c.bcf" | gzip > "$OUT/chr$c.tsv.gz"
done
