#!/bin/sh -e
# shellcheck disable=SC2086
[ -z "$MMSEQS" ] && echo "Please set the environment variable \$MMSEQS to your MMSEQS binary." && exit 1;
[ "$#" -ne 2 ] && echo "Please provide <inputDB> <tmpDir>" && exit 1

notExists() {
	[ ! -f "$1" ]
}

IN="$1"
TMP_PATH="$2"
FOLDSEEK="$(pwd)"/foldseek/bin/foldseek

[ ! -f "${IN}.dbtype" ] && echo "${IN}.dbtype not found!" && exit 1;
if [ -n "${USE_FOLDSEEK}" ]; then 
    [ ! -f "$FOLDSEEK" ] && echo "Please make sure Foldseek is installed in the working directory." && exit 1;
    [ ! -f "${IN}_h.dbtype" ] && echo "${IN}_h.dbtype not found!" && exit 1;
    [ ! -f "${IN}_ss.dbtype" ] && echo "${IN}_ss.dbtype not found!" && exit 1;
    [ ! -f "${IN}_ca.dbtype" ] && echo "${IN}_ca.dbtype not found!" && exit 1;
fi

if notExists "${TMP_PATH}/cluster_mmseqs.dbtype"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" cluster "${IN}" "${TMP_PATH}/cluster_mmseqs" "${TMP_PATH}/cluster" ${CLUSTER_PAR} \
        || fail "cluster failed"
fi

if notExists "${TMP_PATH}/db_clu_1.dbtype"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" createsubdb "${TMP_PATH}/cluster_mmseqs" "${IN}" "${TMP_PATH}/db_clu_1" ${VERBOSITY}\
        || fail "createsubdb failed"
fi

if notExists "${TMP_PATH}/db_clu_1_profile.dbtype"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" result2profile "${TMP_PATH}/db_clu_1" "${IN}" "${TMP_PATH}/cluster_mmseqs" "${TMP_PATH}/db_clu_1_profile" ${THREADS_PAR}\
        || fail "result2profile failed"
fi

if notExists "${TMP_PATH}/db_clu_1_consensus.dbtype"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" profile2consensus "${TMP_PATH}/db_clu_1_profile" "${TMP_PATH}/db_clu_1_consensus" ${CONSENSUS_PAR}\
        || fail "profile2consensus failed"
fi

if notExists "${TMP_PATH}/db_clu_mmseqs_aln.dbtype"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" align "${IN}" "${IN}" "${TMP_PATH}/cluster_mmseqs" "${TMP_PATH}/db_clu_mmseqs_aln" -a ${THREADS_PAR} \
        || fail "align failed"
fi

if [ -n "${USE_FOLDSEEK}" ]; then 
    if notExists "${TMP_PATH}/db_clu_1_ss.dbtype"; then
        # shellcheck disable=SC2086
        "${MMSEQS}" createsubdb "${TMP_PATH}/cluster_mmseqs" "${IN}_ss" "${TMP_PATH}/db_clu_1_ss" ${VERBOSITY}\
            || fail "createsubdb failed"
    fi

    if notExists "${TMP_PATH}/db_clu_1_ca.dbtype"; then
        # shellcheck disable=SC2086
        "${MMSEQS}" createsubdb "${TMP_PATH}/cluster_mmseqs" "${IN}_ca" "${TMP_PATH}/db_clu_1_ca" ${VERBOSITY}\
            || fail "createsubdb failed"
    fi

    if notExists "${TMP_PATH}/cluster_foldseek.dbtype"; then
        # shellcheck disable=SC2086
        "${FOLDSEEK}" cluster "${TMP_PATH}/db_clu_1" "${TMP_PATH}/cluster_foldseek" "${TMP_PATH}/cluster" --cov-mode 0 -c 0.8 ${THREADS_PAR} \
            || fail "foldseek cluster failed"
    fi

    if notExists "${TMP_PATH}/cluster_merge.dbtype"; then
        # shellcheck disable=SC2086
        "${FOLDSEEK}" mergeclusters "${IN}" "${TMP_PATH}/cluster_merge" "${TMP_PATH}/cluster_foldseek" "${TMP_PATH}/cluster_mmseqs" $MERGECLU_PAR \
            || fail "mergeclusters died"
    fi

    if notExists "${TMP_PATH}/db_clu_2.dbtype"; then
        # shellcheck disable=SC2086
        "${FOLDSEEK}" createsubdb "${TMP_PATH}/cluster_merge" "${IN}" "${TMP_PATH}/db_clu_2" ${VERBOSITY}\
            || fail "createsubdb failed"
    fi

    if notExists "${TMP_PATH}/db_clu_2_ss.dbtype"; then
        # shellcheck disable=SC2086
        "${FOLDSEEK}" createsubdb "${TMP_PATH}/cluster_merge" "${IN}_ss" "${TMP_PATH}/db_clu_2_ss" ${VERBOSITY}\
            || fail "createsubdb failed"
    fi

    if notExists "${TMP_PATH}/db_clu_2_ca.dbtype"; then
        # shellcheck disable=SC2086
        "${FOLDSEEK}" createsubdb "${TMP_PATH}/cluster_merge" "${IN}_ca" "${TMP_PATH}/db_clu_2_ca" ${VERBOSITY}\
            || fail "createsubdb failed"
    fi

    if notExists "${TMP_PATH}/db_clu_2_profile.dbtype"; then
        # shellcheck disable=SC2086
        "${FOLDSEEK}" result2profile "${TMP_PATH}/db_clu_2" "${IN}" "${TMP_PATH}/cluster_merge" "${TMP_PATH}/db_clu_2_profile" ${THREADS_PAR}\
            || fail "result2profile failed"
    fi

    if notExists "${IN}_clu_foldseek_consensus.dbtype"; then
        # shellcheck disable=SC2086
        "${FOLDSEEK}" profile2consensus "${TMP_PATH}/db_clu_2_profile" "${IN}_clu_foldseek_consensus" ${CONSENSUS_PAR}\
            || fail "profile2consensus failed"
    fi

    if notExists "${TMP_PATH}/db_clu_2_profile_ss.dbtype"; then
        # shellcheck disable=SC2086
        "${FOLDSEEK}" result2profile "${TMP_PATH}/db_clu_2_ss" "${IN}_ss" "${TMP_PATH}/cluster_merge" "${TMP_PATH}/db_clu_2_profile_ss" ${THREADS_PAR}\
            || fail "result2profile failed"
    fi

    if notExists "${IN}_clu_foldseek_consensus_ss.dbtype"; then
        # shellcheck disable=SC2086
        "${FOLDSEEK}" profile2consensus "${TMP_PATH}/db_clu_2_profile_ss" "${IN}_clu_foldseek_consensus_ss" ${THREADS_PAR}\
            || fail "profile2consensus failed"
    fi

    if notExists "${IN}_clu_foldseek_aln.dbtype"; then
        # shellcheck disable=SC2086
        "${MMSEQS}" align "${IN}" "${IN}" "${TMP_PATH}/cluster_merge" "${IN}_clu_foldseek_aln" -a ${THREADS_PAR} \
            || fail "align failed"
    fi
else
    if notExists "${IN}_clu_aln.dbtype"; then
        "${MMSEQS}" mvdb "${TMP_PATH}/db_clu_mmseqs_aln" "${IN}_clu_aln" ${VERBOSITY} \
            || fail "mvdb failed"
    fi

    if notExists "${IN}_clu_aln.dbtype"; then
        "${MMSEQS}" mvdb "${TMP_PATH}/db_clu_1_consensus" "${IN}_clu_consensus" ${VERBOSITY} \
            || fail "mvdb failed"
    fi

    if notExists "${IN}_clu_rep_profile.dbtype"; then
        "${MMSEQS}" mvdb "${TMP_PATH}/db_clu_1_profile" "${IN}_clu_rep_profile" ${VERBOSITY} \
            || fail "mvdb failed"
    fi
fi

if [ -n "${REMOVE_TMP}" ]; then
    echo "Remove temporary files"
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/db_clu_1" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/db_clu_1_profile" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/db_clu_1_consensus" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/db_clu_mmseqs_aln" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/cluster_mmseqs" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/cluster_foldseek" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/cluster_merge" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/db_clu_2" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/db_clu_2_ss" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/db_clu_2_ca" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/db_clu_2_profile" ${VERBOSITY}
    rm -rf "${TMP_PATH}/cluster"
    rm -f "${TMP_PATH}/clusterdb.sh"
fi