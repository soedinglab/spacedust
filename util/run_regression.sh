#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

SPACEDUST="$1"
DATA="$2"
BASEDIR="$3"

mkdir -p "${BASEDIR}"

"${SPACEDUST}" createsetdb ${DATA}/*.faa "${BASEDIR}/genome" "${BASEDIR}/tmp"
"${SPACEDUST}" clustersearch "${BASEDIR}/genome" "${BASEDIR}/genome" "${BASEDIR}/result.tsv" "${BASEDIR}/tmp" --filter-self-match

awk '/^>/ { cnt++; } END { if (cnt != 308) exit 1; }' "${BASEDIR}/result.tsv" \
  || fail "Check 1 failed. Expected: 308 Actual: $(awk '/^>/ { cnt++; } END { print cnt; }' "${BASEDIR}/result.tsv")"
awk '/^#/ { if ($4 < 1E-20) cnt++; } END { if (cnt != 2) exit 1; }' "${BASEDIR}/result.tsv" \
  || fail "Check 2 failed. Expected: 2 Actual: $(awk '/^#/ { if ($4 < 1E-20) cnt++; } END { print cnt; }' "${BASEDIR}/result.tsv")"