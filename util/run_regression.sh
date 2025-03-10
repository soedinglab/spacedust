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

#foldseek related regression
"${SPACEDUST}" createsetdb ${DATA}/foldseek_testdb/foldseek_test "${BASEDIR}/foldseek_testdb" "${BASEDIR}/tmp"
"${SPACEDUST}" clustersearch "${BASEDIR}/foldseek_testdb" "${BASEDIR}/foldseek_testdb" "${BASEDIR}/result_foldseek.tsv" "${BASEDIR}/tmp" --filter-self-match --search-mode 2 

awk '/^>/ { cnt++; } END { if (cnt != 308) exit 1; }' "${BASEDIR}/result.tsv" \
  || fail "Check 1 failed. Expected: 308 Actual: $(awk '/^>/ { cnt++; } END { print cnt; }' "${BASEDIR}/result.tsv")"
awk '/^#/ { if ($4 < 1E-20) cnt++; } END { if (cnt != 2) exit 1; }' "${BASEDIR}/result.tsv" \
  || fail "Check 2 failed. Expected: 2 Actual: $(awk '/^#/ { if ($4 < 1E-20) cnt++; } END { print cnt; }' "${BASEDIR}/result.tsv")"

#foldseek related checks
awk '/^>/ { cnt++; } END { if (cnt != 568) exit 1; }' "${BASEDIR}/result_foldseek.tsv" \
  || fail "Check 3 failed. Expected: 568 Actual: $(awk '/^>/ { cnt++; } END { print cnt; }' "${BASEDIR}/result_foldseek.tsv")"
