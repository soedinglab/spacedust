#!/bin/sh -e
# shellcheck disable=SC2086
[ -z "$MMSEQS" ] && echo "Please set the environment variable \$MMSEQS to your MMSEQS binary." && exit 1;
[ "$#" -ne 3 ] && echo "Please provide <inputDB> <targetFoldSeekDB> <tmpDir>" && exit 1

notExists() {
	[ ! -f "$1" ]
}

IN="$1"
TARGET="$2"
OUT="$1_foldseek"
TMP_PATH="$3"

[ ! -f "${IN}.dbtype" ] && echo "${IN}.dbtype not found!" && exit 1;
[ ! -f "${TARGET}.dbtype" ] && echo "${TARGET}.dbtype not found!" && exit 1;
[ ! -f "${TARGET}_h.dbtype" ] && echo "${TARGET}_h.dbtype not found!" && exit 1;
[ ! -f "${TARGET}_ss.dbtype" ] && echo "${TARGET}_ss.dbtype not found!" && exit 1;
[ ! -f "${TARGET}_ca.dbtype" ] && echo "${TARGET}_ca.dbtype not found!" && exit 1;


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

#fake a dbtype that MMseqs accepts
mv -f "${TARGET}_ca.dbtype" "${TARGET}_ca.dbtype.tmp"
cp -f "${TARGET}_ss.dbtype" "${TARGET}_ca.dbtype"

if notExists "${TMP_PATH}/outDB_ca.dbtype"; then
    # shellcheck disable=SC2086 
    "$MMSEQS" createsubdb "${TMP_PATH}/topHitalnSwapDB" "${TARGET}_ca" "${TMP_PATH}/outDB_ca" ${VERBOSITY}\
        || fail "createsubdb failed"
fi

if notExists "${OUT}_ca.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" renamedbkeys "${TMP_PATH}/mapping" "${TMP_PATH}/outDB_ca" "${OUT}_ca" ${VERBOSITY}\
        || fail "createsubdb failed"
fi
mv -f "${TARGET}_ca.dbtype.tmp" "${TARGET}_ca.dbtype"
cp -f "${TARGET}_ca.dbtype" "${OUT}_ca.dbtype"

"$MMSEQS" lndb "${IN}_h" "${OUT}_h" ${VERBOSITY}

#use the original lookup file to ensure fixed indeces
"$MMSEQS" lndb "${IN}.lookup" "${OUT}.lookup" ${VERBOSITY}

"$MMSEQS" lndb "${IN}.source" "${OUT}.source" ${VERBOSITY}

TOTAL_NUM_SEQS=$(wc -l < "${IN}.index")
NUM_SEQS_MAPPED="$(wc -l < "${OUT}.index")"
PERCENTAGE=$(echo "scale=2; $NUM_SEQS_MAPPED / $TOTAL_NUM_SEQS * 100" | bc)
echo "${NUM_SEQS_MAPPED} out of ${TOTAL_NUM_SEQS} sequences (${PERCENTAGE}%) were mapped to target DB (${TARGET}) and stored in ${OUT}"

if notExists "${OUT}_unmapped.index"; then
    awk 'BEGIN{FS="\t"}NR==FNR{a[$1];next} !($1 in a)' "${OUT}.index" "${IN}.index" > "${TMP_PATH}/unmapped_id"
    # shellcheck disable=SC2086 
    "$MMSEQS" createsubdb "${TMP_PATH}/unmapped_id" "${IN}" "${IN}_unmapped" ${VERBOSITY}\
        || fail "createsubdb failed"
    rm -f "${TMP_PATH}/unmapped_id"
    echo "The remaining unmapped sequences are stored in ${IN}_unmapped."
fi

if [ -n "${REMOVE_TMP}" ]; then
    echo "Remove temporary files"
    rm -rf "${TMP_PATH}/search"
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/alnDB" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/topHitalnDB" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/topHitalnSwapDB" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/topHitalnSwapDB_pref" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/outDB" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/outDB_h" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/outDB_ss" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/outDB_ca" ${VERBOSITY}
    rm -f "${TMP_PATH}/mapping"
    rm -f "${TMP_PATH}/aa2foldseek.sh"
fi