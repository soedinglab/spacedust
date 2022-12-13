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
"${SPACEDUST}" clustersearch "${BASEDIR}/genome" "${BASEDIR}/genome" "${BASEDIR}/result" "${BASEDIR}/tmp" --filter-self-match


awk 'END{ if (NR != 139) exit 1; }' "${BASEDIR}/result" \
  || fail "Check 1 failed"
tr -d '\000' < "${BASEDIR}/result_h" |awk '$3 < 1E-20 { cnt++; } END { if (cnt != 2) exit 1; }' \
  || fail "Check 2 failed"
