#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

notDone() {
    [ ! -f "${TMP_PATH}/$1.done" ]
}

markDone() {
    touch "${TMP_PATH}/$1.done"
}

freeDb() {
    [ -n "${REMOVE_TMP}" ] || return 0
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/$1" ${VERBOSITY}
    if [ -f "${TMP_PATH}/$1_h.dbtype" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/$1_h" ${VERBOSITY}
    fi
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
            if notDone "result"; then
                # shellcheck disable=SC2086
                "${FOLDSEEK}" search "${QUERY}" "${TARGET}_clu" "${TMP_PATH}/result" "${TMP_PATH}/search" --cluster-search 1 ${FOLDSEEKSEARCH_PAR}\
                    || fail "foldseek search failed"
                markDone "result"
            fi
        else
            if notDone "result_foldseek"; then
                # shellcheck disable=SC2086
                "${FOLDSEEK}" search "${QUERY}_foldseek" "${TARGET}_foldseek_clu" "${TMP_PATH}/result_foldseek" "${TMP_PATH}/search" --cluster-search 1 ${FOLDSEEKSEARCH_PAR}\
                    || fail "foldseek search failed"
                markDone "result_foldseek"
            fi
            if notDone "result_clu"; then
                # shellcheck disable=SC2086
                "${MMSEQS}" search "${QUERY}_unmapped" "${TARGET}_clu" "${TMP_PATH}/result_clu" "${TMP_PATH}/search" ${SEARCH_PAR} \
                    || fail "mmseqs search failed"
                markDone "result_clu"
            fi
            if notDone "result_exp"; then
                # shellcheck disable=SC2086
                "${MMSEQS}" expandaln "${QUERY}_unmapped" "${TARGET}_clu" "${TMP_PATH}/result_clu" "${TARGET}_clu_aln" "${TMP_PATH}/result_exp" ${THREADS_PAR} \
                    || fail "expandaln failed"
                freeDb "result_clu"
                markDone "result_exp"
            fi
            if notDone "result_mmseqs"; then
                # shellcheck disable=SC2086
                "${MMSEQS}" align "${QUERY}_unmapped" "${TARGET}" "${TMP_PATH}/result_exp" "${TMP_PATH}/result_mmseqs" -a --alt-ali 10 ${THREADS_PAR} \
                    || fail "realign failed"
                freeDb "result_exp"
                markDone "result_mmseqs"
            fi
            if notDone "result"; then
                # shellcheck disable=SC2086
                "${MMSEQS}" concatdbs "${TMP_PATH}/result_foldseek" "${TMP_PATH}/result_mmseqs" "${TMP_PATH}/result" --preserve-keys ${THREADS_PAR} \
                    || fail "concatdbs failed"
                freeDb "result_foldseek"
                freeDb "result_mmseqs"
                markDone "result"
            fi
        fi
    else
        if notDone "result_clu"; then
            # shellcheck disable=SC2086
            "${MMSEQS}" search "${QUERY}" "${TARGET}_clu_rep_profile" "${TMP_PATH}/result_clu" "${TMP_PATH}/search" ${SEARCH_PAR} \
                || fail "search failed"
            markDone "result_clu"
        fi

        #realignment?
        if notDone "result"; then
            # shellcheck disable=SC2086
            "${MMSEQS}" expandaln "${QUERY}" "${TARGET}_clu_rep_profile" "${TMP_PATH}/result_clu" "${TARGET}_clu_aln" "${TMP_PATH}/result" ${THREADS_PAR} \
                || fail "expandaln failed"
            freeDb "result_clu"
            markDone "result"
        fi
    fi

else
    if notDone "result"; then
        if [ -n "${USE_FOLDSEEK}" ]; then
            if [ -n "${USE_PROSTT5}" ] ; then
                if notDone "result"; then
                    # shellcheck disable=SC2086
                    "${FOLDSEEK}" search "${QUERY}" "${TARGET}" "${TMP_PATH}/result" "${TMP_PATH}/search" ${FOLDSEEKSEARCH_PAR}\
                        || fail "foldseek search failed"
                    markDone "result"
                fi
            else
                if notDone "result_foldseek"; then
                    # shellcheck disable=SC2086
                    "${FOLDSEEK}" search "${QUERY}_foldseek" "${TARGET}_foldseek" "${TMP_PATH}/result_foldseek" "${TMP_PATH}/search" ${FOLDSEEKSEARCH_PAR}\
                        || fail "foldseek search failed"
                    markDone "result_foldseek"
                fi
                if notDone "result_mmseqs"; then
                    # shellcheck disable=SC2086
                    "${MMSEQS}" search "${QUERY}_unmapped" "${TARGET}" "${TMP_PATH}/result_mmseqs" "${TMP_PATH}/search" ${SEARCH_PAR} \
                        || fail "mmseqs search failed"
                    markDone "result_mmseqs"
                fi
                if notDone "result"; then
                    # shellcheck disable=SC2086
                    "${MMSEQS}" concatdbs "${TMP_PATH}/result_foldseek" "${TMP_PATH}/result_mmseqs" "${TMP_PATH}/result" --preserve-keys ${THREADS_PAR} \
                        || fail "concatdbs failed"
                    freeDb "result_foldseek"
                    freeDb "result_mmseqs"
                    markDone "result"
                fi
            fi
        else
            if notDone "result"; then
            # shellcheck disable=SC2086
            "${MMSEQS}" search "${QUERY}" "${TARGET}" "${TMP_PATH}/result" "${TMP_PATH}/search" ${SEARCH_PAR} \
                || fail "mmseqs search failed"
                markDone "result"
            fi
        fi
    fi
fi

# the search tmp is dead as soon as search returns
if [ -n "${REMOVE_TMP}" ]; then
    rm -rf "${TMP_PATH}/search"
fi

if notDone "aggregate"; then
    # aggregation: take for each target set the best hit
    # shellcheck disable=SC2086
    "${MMSEQS}" besthitbyset "${QUERY}" "${TARGET}" "${TMP_PATH}/result" "${TMP_PATH}/aggregate" ${BESTHITBYSET_PAR} \
        || fail "aggregate best hit failed"
    markDone "aggregate"
fi

if notDone "aggregate_merged"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" mergeresultsbyset "${QUERY}_set_to_member" "${TMP_PATH}/aggregate" "${TMP_PATH}/aggregate_merged" ${THREADS_PAR} \
        || fail "mergesetresults failed"
    freeDb "aggregate"
    markDone "aggregate_merged"
fi

if notDone "matches"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" combinehits "${QUERY}" "${TARGET}" "${TMP_PATH}/aggregate_merged" "${TMP_PATH}/matches" "${TMP_PATH}" ${COMBINEHITS_PAR} \
        || fail "combinepvalperset failed"
    freeDb "aggregate_merged"
    markDone "matches"
fi

if notDone "clusters"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" clusterhits "${QUERY}" "${TARGET}" "${TMP_PATH}/matches" "${TMP_PATH}/clusters" ${CLUSTERHITS_PAR} \
        || fail "clusterhits failed"
    freeDb "matches"
    markDone "clusters"
fi

# shellcheck disable=SC2086
"${MMSEQS}" summarizeresults "${QUERY}" "${TARGET}" "${TMP_PATH}/clusters" "${OUTPUT}" ${THREADS_PAR} \
    || fail "summarizeresults failed"

#postprocessing
if notDone "clu_to_seq"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" filterdb "${TMP_PATH}/clusters" "${TMP_PATH}/clu_to_seq" --trim-to-one-column ${THREADS_PAR} \
        || fail "filterdb failed"
    markDone "clu_to_seq"
fi

# shellcheck disable=SC2086
"${MMSEQS}" swapdb "${TMP_PATH}/clu_to_seq" "${OUTPUT}_seq_to_clu" ${THREADS_PAR} \
    || fail "swapdb failed"

if [ -n "${REMOVE_TMP}" ]; then
    echo "Remove temporary files"
    rm -rf "${TMP_PATH}/search"
    freeDb "result"
    freeDb "clusters"
    freeDb "result_foldseek"
    freeDb "result_mmseqs"
    freeDb "result_clu"
    freeDb "result_exp"
    freeDb "aggregate"
    freeDb "aggregate_merged"
    freeDb "matches"
    freeDb "clu_to_seq"
    rm -f "${TMP_PATH}"/*.done
    rm -f "${TMP_PATH}/clustersearch.sh"
fi

