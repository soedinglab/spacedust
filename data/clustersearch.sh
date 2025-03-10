#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

#pre processing
[ -z "$MMSEQS" ] && echo "Please set the environment variable \$MMSEQS to your MMSEQS binary." && exit 1;
# check number of input variables
[ "$#" -ne 4 ] && echo "Please provide <queryDB> <targetDB> <outputDB> <tmpDir>" && exit 1;
# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[ ! -f "$2.dbtype" ] && echo "$2.dbtype not found!" && exit 1;
# TO DO??? add check if $3.dbtype already exists before entire workfolw ???

QUERY="$1"
TARGET="$2"
OUTPUT="$3"
TMP_PATH="$4"

if [ -n "${USE_PROSTT5}" ]; then 
    [ -n "${USE_PROFILE}" ] && [ ! -f "${TARGET}_clu_seq.dbtype" ] && echo "${TARGET}_foldseek_clu_seq.dbtype not found! Please make sure the ${TARGET} is clustered with clusterdb ${TARGET} tmp --search-mode 1" && exit 1;
    [ ! -f "${TARGET}_ss.dbtype" ] && echo "${TARGET}_ss.dbtype not found! Please make sure the ${TARGET} is created using ProstT5. " && exit 1;
elif [ -n "${USE_FOLDSEEK}" ]; then 
    [ -n "${USE_PROFILE}" ] && [ ! -f "${TARGET}_foldseek_clu_seq.dbtype" ] && echo "${TARGET}_foldseek_clu_seq.dbtype not found! Please make sure the ${TARGET}_foldseek is clustered with clusterdb ${TARGET}_foldseek tmp --search-mode 1" && exit 1;
    [ ! -f "${TARGET}_foldseek.dbtype" ] && echo "${TARGET}_foldseek.dbtype not found! Please make sure the ${TARGET}_foldseek is created with aa2foldseek. If ${TARGET} is created with ProstT5 please use --search-mode 2" && exit 1;
fi

if [ -n "${USE_PROFILE}" ]; then
    if [ -n "${USE_FOLDSEEK}" ]; then
        if [ -n "${USE_PROSTT5}" ] ; then
            if notExists "${TMP_PATH}/result.index"; then
                # shellcheck disable=SC2086
                "${FOLDSEEK}" search "${QUERY}" "${TARGET}_clu" "${TMP_PATH}/result" "${TMP_PATH}/search" --cluster-search 1 ${FOLDSEEKSEARCH_PAR}\
                    || fail "foldseek search failed"
            fi
        else
            if notExists "${TMP_PATH}/result_foldseek.index"; then
                # shellcheck disable=SC2086
                "${FOLDSEEK}" search "${QUERY}_foldseek" "${TARGET}_foldseek_clu" "${TMP_PATH}/result_foldseek" "${TMP_PATH}/search" --cluster-search 1 ${FOLDSEEKSEARCH_PAR}\
                    || fail "foldseek search failed"
            fi
            if notExists "${TMP_PATH}/result_clu.index"; then
                # shellcheck disable=SC2086
                "${MMSEQS}" search "${QUERY}_unmapped" "${TARGET}_clu" "${TMP_PATH}/result_clu" "${TMP_PATH}/search" ${SEARCH_PAR} \
                    || fail "mmseqs search failed"
            fi
            if notExists "${TMP_PATH}/result_exp.index"; then
                # shellcheck disable=SC2086
                "${MMSEQS}" expandaln "${QUERY}_unmapped" "${TARGET}_clu" "${TMP_PATH}/result_clu" "${TARGET}_clu_aln" "${TMP_PATH}/result_exp" ${THREADS_PAR} \
                    || fail "expandaln failed"
            fi
            if notExists "${TMP_PATH}/result_mmseqs.index"; then
                # shellcheck disable=SC2086
                "${MMSEQS}" align "${QUERY}_unmapped" "${TARGET}" "${TMP_PATH}/result_exp" "${TMP_PATH}/result_mmseqs" -a --alt-ali 10 ${THREADS_PAR} \
                    || fail "realign failed"
            fi
            if notExists "${TMP_PATH}/result.index"; then
                # shellcheck disable=SC2086
                "${MMSEQS}" concatdbs "${TMP_PATH}/result_foldseek" "${TMP_PATH}/result_mmseqs" "${TMP_PATH}/result" --preserve-keys ${THREADS_PAR} \
                    || fail "concatdbs failed"
            fi
        fi
    else
        if notExists "${TMP_PATH}/result_clu.index"; then
            # shellcheck disable=SC2086
            "${MMSEQS}" search "${QUERY}" "${TARGET}_clu_rep_profile" "${TMP_PATH}/result_clu" "${TMP_PATH}/search" ${SEARCH_PAR} \
                || fail "search failed"
        fi

        #realignment?
        if notExists "${TMP_PATH}/result.index"; then
            # shellcheck disable=SC2086
            "${MMSEQS}" expandaln "${QUERY}" "${TARGET}_clu_rep_profile" "${TMP_PATH}/result_clu" "${TARGET}_clu_aln" "${TMP_PATH}/result" ${THREADS_PAR} \
                || fail "expandaln failed"
        fi
    fi

else
    if notExists "${TMP_PATH}/result.index"; then
        if [ -n "${USE_FOLDSEEK}" ]; then
            if [ -n "${USE_PROSTT5}" ] ; then
                if notExists "${TMP_PATH}/result.index"; then
                    # shellcheck disable=SC2086
                    "${FOLDSEEK}" search "${QUERY}" "${TARGET}" "${TMP_PATH}/result" "${TMP_PATH}/search" ${FOLDSEEKSEARCH_PAR}\
                        || fail "foldseek search failed"
                fi
            else
                if notExists "${TMP_PATH}/result_foldseek.index"; then
                    # shellcheck disable=SC2086
                    "${FOLDSEEK}" search "${QUERY}_foldseek" "${TARGET}_foldseek" "${TMP_PATH}/result_foldseek" "${TMP_PATH}/search" ${FOLDSEEKSEARCH_PAR}\
                        || fail "foldseek search failed"
                fi
                if notExists "${TMP_PATH}/result_mmseqs.index"; then
                    # shellcheck disable=SC2086
                    "${MMSEQS}" search "${QUERY}_unmapped" "${TARGET}" "${TMP_PATH}/result_mmseqs" "${TMP_PATH}/search" ${SEARCH_PAR} \
                        || fail "mmseqs search failed"
                fi
                if notExists "${TMP_PATH}/result.index"; then
                    # shellcheck disable=SC2086
                    "${MMSEQS}" concatdbs "${TMP_PATH}/result_foldseek" "${TMP_PATH}/result_mmseqs" "${TMP_PATH}/result" --preserve-keys ${THREADS_PAR} \
                        || fail "concatdbs failed"
                fi
            fi
        else
            if notExists "${TMP_PATH}/result.index"; then
            # shellcheck disable=SC2086
            "${MMSEQS}" search "${QUERY}" "${TARGET}" "${TMP_PATH}/result" "${TMP_PATH}/search" ${SEARCH_PAR} \
                || fail "mmseqs search failed"
            fi
        fi
    fi
fi

if notExists "${TMP_PATH}/result_prefixed.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" prefixid  "${TMP_PATH}/result" "${TMP_PATH}/result_prefixed" ${THREADS_PAR} \
        || fail "prefixid failed"
fi

if notExists "${TMP_PATH}/aggregate.index"; then
    # aggregation: take for each target set the best hit
    # shellcheck disable=SC2086
    "${MMSEQS}" besthitbyset "${QUERY}" "${TARGET}" "${TMP_PATH}/result_prefixed" "${TMP_PATH}/aggregate" ${BESTHITBYSET_PAR} \
        || fail "aggregate best hit failed"
fi

if notExists "${TMP_PATH}/aggregate_merged.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" mergeresultsbyset "${QUERY}_set_to_member" "${TMP_PATH}/aggregate" "${TMP_PATH}/aggregate_merged" ${THREADS_PAR} \
        || fail "mergesetresults failed"
fi

if notExists "${TMP_PATH}/matches.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" combinehits "${QUERY}" "${TARGET}" "${TMP_PATH}/aggregate_merged" "${TMP_PATH}/matches" "${TMP_PATH}" ${COMBINEHITS_PAR} \
        || fail "combinepvalperset failed"
fi

if notExists "${TMP_PATH}/clusters.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" clusterhits "${QUERY}" "${TARGET}" "${TMP_PATH}/matches" "${TMP_PATH}/clusters" ${CLUSTERHITS_PAR} \
        || fail "clusterhits failed"
fi

# shellcheck disable=SC2086
"${MMSEQS}" summarizeresults "${QUERY}" "${TARGET}" "${TMP_PATH}/clusters" "${OUTPUT}" ${THREADS_PAR} \
    || fail "summarizeresults failed"

#postprocessing
if notExists "${TMP_PATH}/clu_to_seq.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" filterdb "${TMP_PATH}/clusters" "${TMP_PATH}/clu_to_seq" --trim-to-one-column ${THREADS_PAR} \
        || fail "filterdb failed"
fi

if notExists "${TMP_PATH}/seq_to_clu.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" swapdb "${TMP_PATH}/clu_to_seq" "${OUTPUT}_seq_to_clu" ${THREADS_PAR} \
        || fail "swapdb failed"
fi

if [ -n "${REMOVE_TMP}" ]; then
    echo "Remove temporary files"
    rm -rf "${TMP_PATH}/search"
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/result" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/result_prefixed" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/aggregate" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/aggregate_merged" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/matches" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/clusters" ${VERBOSITY}
    rm -f "${TMP_PATH}/clustersearch.sh"
fi

