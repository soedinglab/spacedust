#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

SPACEPHARER="$1"
DATA="$2"
BASEDIR="$3"

mkdir -p "${BASEDIR}"

"${SPACEPHARER}" createsetdb ${DATA}/GCA*.fna.gz "${BASEDIR}/targets" "${BASEDIR}/tmp" --tax-mapping-file "${DATA}/genome_taxa.tsv" --ncbi-tax-dump "${DATA}/ncbi_taxdump"
"${SPACEPHARER}" createsetdb ${DATA}/GCA*.fna.gz "${BASEDIR}/targets_rev" "${BASEDIR}/tmp" --reverse-fragments 1
"${SPACEPHARER}" easy-predict ${DATA}/*.fas "${BASEDIR}/targets" "${BASEDIR}/result.tsv" "${BASEDIR}/tmp" --tax-mapping-file "${DATA}/spacer_taxa.tsv" --ncbi-tax-dump "${DATA}/ncbi_taxdump"
"${SPACEPHARER}" parsespacer ${DATA}/*_test "${BASEDIR}/query"


awk '/^>/ && $3 < 1E-03 { cnt++; } END { if (cnt != 6) exit 1; }' "${BASEDIR}/result.tsv" \
  || fail "Check 1 failed"
awk 'BEGIN { other = 0; } $2 == 40521 { listeria++; next } $2 == 244310 { burkholdia++; next; } $2 != 0 { other++; } END { if ((listeria" "burkholdia" "other) != "5 1 0") exit 1; }' "${BASEDIR}/result.tsv_lca.tsv" \
  || fail "Check 2 failed"
awk '$1 == "GCA_000836905.1_ViralProj14035_genomic.fna.gz" && $2 == 0 { c++ } $1 == "GCA_000845445.1_ViralProj14409_genomic.fna.gz" && $2 == 28216 { c++ } $1 == "GCA_000849645.1_ViralProj14589_genomic.fna.gz" && $2 == 1639 { c++ } END { if (c != 3) exit 1; }' "${BASEDIR}/result.tsv_lca_per_target.tsv" \
  || fail "Check 3 failed"
awk 'BEGIN{ i = 0 }{ i++; } END{ if (i != 127) exit 1;}' "${BASEDIR}/query.index" \
  || fail "Check 4 failed"