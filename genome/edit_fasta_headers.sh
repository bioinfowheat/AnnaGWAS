#!/bin/bash
# Edit FASTA headers for the P. aegeria reference genome
# - NC_ headers (chromosomes): replace with chromosome number/letter (e.g., >1, >W, >Z)
# - NW_ headers (scaffolds): keep accession ID only, trim at first space
#
# Usage: bash edit_fasta_headers.sh input.fna output.edited.fna

INPUT="${1:-GCF_905163445.1_ilParAegt1.1_genomic.fna}"
OUTPUT="${2:-GCF_905163445.1_ilParAegt1.1_genomic.edited.fna}"

awk '{
  if (/^>NC_/) {
    # Find the word just before the first comma (chromosome number/letter)
    for (i=2; i<=NF; i++) {
      if ($i ~ /,$/) {
        gsub(/,/, "", $i)
        print ">" $i
        break
      }
    }
  } else if (/^>NW_/) {
    # Keep accession ID only (trim at first space)
    print $1
  } else {
    print
  }
}' "$INPUT" > "$OUTPUT"

echo "Done. Edited headers written to: $OUTPUT"
echo "Header count: $(grep -c '^>' "$OUTPUT")"
