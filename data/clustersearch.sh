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
FOLDSEEK="$(pwd)"/foldseek/bin/foldseek

if [ -n "${USE_FOLDSEEK}" ]; then 
    [ -n "${USE_PROFILE}" ] && echo "Profile cluster search with Foldseek is not supported." && exit 1;
    [ ! -f "$FOLDSEEK" ] && echo "Please make sure Foldseek is installed in the working directory." && exit 1;
    [ ! -f "${TARGET}_foldseek.dbtype" ] && echo "${TARGET}_foldseek.dbtype not found! Please make sure the ${TARGET}_foldseek is created with aa2foldseek" && exit 1;
fi

if [ -n "${USE_PROFILE}" ]; then
    if notExists "${TARGET}_clu.index"; then
        # shellcheck disable=SC2086
        "${MMSEQS}" cluster "${TARGET}" "${TARGET}_clu" "${TMP_PATH}/cluster" ${CLUSTER_PAR} \
            || fail "cluster failed"
    fi

    if notExists "${TARGET}_clu_rep.index"; then
        # shellcheck disable=SC2086
        "${MMSEQS}" createsubdb "${TARGET}_clu" "${TARGET}" "${TARGET}_clu_rep" ${VERBOSITY}\
            || fail "createsubdb failed"
    fi

    if notExists "${TARGET}_clu_rep_profile.index"; then
        # shellcheck disable=SC2086
        "${MMSEQS}" result2profile "${TARGET}_clu_rep" "${TARGET}" "${TARGET}_clu" "${TARGET}_clu_rep_profile" ${THREADS_PAR}\
            || fail "result2profile failed"
    fi

    if notExists "${TMP_PATH}/result_clu.index"; then
        # shellcheck disable=SC2086
        "${MMSEQS}" search "${QUERY}" "${TARGET}_clu_rep_profile" "${TMP_PATH}/result_clu" "${TMP_PATH}/search" ${SEARCH_PAR} \
            || fail "search failed"
    fi

    #expandaln?
    if notExists "${TARGET}_clu_aln.index"; then
        # shellcheck disable=SC2086
        "${MMSEQS}" align "${TARGET}" "${TARGET}" "${TARGET}_clu" "${TARGET}_clu_aln" -a ${THREADS_PAR} \
            || fail "align failed"
    fi

    if notExists "${TMP_PATH}/result.index"; then
        # shellcheck disable=SC2086
        "${MMSEQS}" expandaln "${QUERY}" "${TARGET}_clu_rep_profile" "${TMP_PATH}/result_clu" "${TARGET}_clu_aln" "${TMP_PATH}/result" ${THREADS_PAR} \
            || fail "expandaln failed"
    fi

else
    if notExists "${TMP_PATH}/result.index"; then
        if [ -n "${USE_FOLDSEEK}" ]; then
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
        else
            # shellcheck disable=SC2086
            "${MMSEQS}" search "${QUERY}" "${TARGET}" "${TMP_PATH}/result" "${TMP_PATH}/search" ${SEARCH_PAR} \
                || fail "mmseqs search failed"
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

