#!/bin/sh -e
# shellcheck disable=SC2086
[ -z "$MMSEQS" ] && echo "Please set the environment variable \$MMSEQS to your MMSEQS binary." && exit 1;
[ "$#" -ne 4 ] && echo "Please provide <inputDB> <targetFoldSeekDB> <outDB> <tmpDir>" && exit 1

notExists() {
	[ ! -f "$1" ]
}

IN="$1"
TARGET="$2"
OUT="$3"
TMP_PATH="$4"

[ ! -f "${IN}.dbtype" ] && echo "${IN}.dbtype not found!" && exit 1;
[ ! -f "${TARGET}.dbtype" ] && echo "${TARGET}.dbtype not found!" && exit 1;
[ ! -f "${TARGET}_h.dbtype" ] && echo "${TARGET}_h.dbtype not found!" && exit 1;
[ ! -f "${TARGET}_ss.dbtype" ] && echo "${TARGET}_ss.dbtype not found!" && exit 1;


if notExists "${TMP_PATH}/alnDB.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" search "${IN}" "${TARGET}" "${TMP_PATH}/alnDB" "${TMP_PATH}/search" ${SEARCH_PAR} \
        || fail "search failed"
    fi

if notExists "${TMP_PATH}/topHitalnDB.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" filterdb "${TMP_PATH}/alnDB" "${TMP_PATH}/topHitalnDB" --extract-lines 1 ${THREADS_PAR}\
        || fail "filterdb failed"
fi

if notExists "${TMP_PATH}/topHitalnSwapDB.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" swapdb "${TMP_PATH}/topHitalnDB" "${TMP_PATH}/topHitalnSwapDB" ${THREADS_PAR}\
        || fail "swapdb failed"
fi

if notExists "${TMP_PATH}/topHitalnSwapDB_pref"; then
    # shellcheck disable=SC2086
    "$MMSEQS" prefixid "${TMP_PATH}/topHitalnSwapDB" "${TMP_PATH}/topHitalnSwapDB_pref" --tsv ${THREADS_PAR}\
        || fail "prefixid failed"
fi

if notExists "${TMP_PATH}/outDB.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createsubdb "${TMP_PATH}/topHitalnSwapDB" "${TARGET}" "${TMP_PATH}/outDB" ${VERBOSITY}\
        || fail "createsubdb failed"
fi

cut -f1,2 "${TMP_PATH}/topHitalnSwapDB_pref" > "${TMP_PATH}/mapping"

if notExists "${OUT}.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" renamedbkeys "${TMP_PATH}/mapping" "${TMP_PATH}/outDB" "${OUT}" ${VERBOSITY}\
        || fail "createsubdb failed"
fi

if notExists "${TMP_PATH}/outDB_ss.dbtype"; then
    # shellcheck disable=SC2086 
    "$MMSEQS" createsubdb "${TMP_PATH}/topHitalnSwapDB" "${TARGET}_ss" "${TMP_PATH}/outDB_ss" ${VERBOSITY}\
        || fail "createsubdb failed"
fi

if notExists "${OUT}_ss.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" renamedbkeys "${TMP_PATH}/mapping" "${TMP_PATH}/outDB_ss" "${OUT}_ss" ${VERBOSITY}\
        || fail "createsubdb failed"
fi

"$MMSEQS" lndb "${IN}_h" "${OUT}_h" ${VERBOSITY}

cut -f2 "${TMP_PATH}/mapping" > "${TMP_PATH}/mapping_id"
awk 'NR==FNR{a[$1];next} ($1 in a)' "${TMP_PATH}/mapping_id" "${IN}.lookup" > "${TMP_PATH}/lookup"

if notExists "${OUT}_member_to_set.index"; then
    awk '{ print $1"\t"$3; }' "${TMP_PATH}/lookup" | sort -k1,1n -k2,2n > "${OUT}_member_to_set.tsv"
    "${MMSEQS}" tsv2db "${OUT}_member_to_set.tsv" "${OUT}_member_to_set" --output-dbtype 5 \
        || fail "tsv2db failed"
fi

if notExists "${OUT}_set_to_member.index"; then
    awk '{ print $3"\t"$1; }' "${TMP_PATH}/lookup" | sort -k1,1n -k2,2n > "${OUT}_set_to_member.tsv"
    "${MMSEQS}" tsv2db "${OUT}_set_to_member.tsv" "${OUT}_set_to_member" --output-dbtype 5 \
        || fail "tsv2db failed"
fi

if notExists "${OUT}_set_size.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" result2stats "${OUT}" "${OUT}" "${OUT}_set_to_member" "${OUT}_set_size" ${RESULT2STATS_PAR} \
        || fail "result2stats failed"
fi

#use the original lookup file to ensure fixed indeces
"$MMSEQS" lndb "${IN}.lookup" "${OUT}.lookup" ${VERBOSITY}

"$MMSEQS" lndb "${IN}.source" "${OUT}.source" ${VERBOSITY}

TOTAL_NUM_SEQS=$(wc -l < "${IN}.index")
NUM_SEQS_MAPPED="$(wc -l < "${OUT}.index")"
PERCENTAGE=$(echo "scale=2; $NUM_SEQS_MAPPED / $TOTAL_NUM_SEQS * 100" | bc)
echo "${NUM_SEQS_MAPPED} out of ${TOTAL_NUM_SEQS} sequences (${PERCENTAGE}%) were mapped to target DB (${TARGET})."