#!/bin/sh -e
# shellcheck disable=SC2086
[ -z "$MMSEQS" ] && echo "Please set the environment variable \$MMSEQS to your MMSEQS binary." && exit 1;
[ "$#" -ne 2 ] && echo "Please provide <inputDB> <tmpDir>" && exit 1

fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

IN="$1"
TMP_PATH="$2"

[ ! -f "${IN}.dbtype" ] && echo "${IN}.dbtype not found!" && exit 1;
if [ -n "${USE_FOLDSEEK}" ]; then 
    [ ! -f "$FOLDSEEK" ] && echo "Please make sure Foldseek is installed in the working directory." && exit 1;
    [ ! -f "${IN}_h.dbtype" ] && echo "${IN}_h.dbtype not found!" && exit 1;
    [ ! -f "${IN}_ss.dbtype" ] && echo "${IN}_ss.dbtype not found!" && exit 1;
fi

if [ -n "${USE_FOLDSEEK}" ]; then 
    if notExists "${TMP_PATH}/cluster_foldseek.dbtype"; then
        # shellcheck disable=SC2086
        "${FOLDSEEK}" cluster "${IN}" "${TMP_PATH}/cluster_foldseek" "${TMP_PATH}/structurecluster" --cov-mode 0 -c 0.9 --min-seq-id 0.3 ${THREADS_PAR} \
            || fail "foldseek cluster failed"
    fi

    if notExists "${IN}_clu_seq.dbtype"; then
        # shellcheck disable=SC2086
        "${FOLDSEEK}" cpdb "${IN}" "${IN}_clu_seq" ${VERBOSITY}
    fi

    if notExists "${IN}_clu_seq_h.dbtype"; then
        # shellcheck disable=SC2086
        "${FOLDSEEK}" cpdb "${IN}_h" "${IN}_clu_seq_h" ${VERBOSITY}
    fi

    if notExists "${IN}_clu_seq_ss.dbtype"; then
        # shellcheck disable=SC2086
        "${FOLDSEEK}" cpdb "${IN}_ss" "${IN}_clu_seq_ss" ${VERBOSITY}
    fi

    if [ -f "${IN}_ca.dbtype" ] && notExists "${IN}_clu_seq_ca.dbtype"; then
        # shellcheck disable=SC2086
        "${FOLDSEEK}" cpdb "${IN}_ca" "${IN}_clu_seq_ca" ${VERBOSITY}
    fi

    if notExists "${IN}_clu_clu.dbtype"; then
        # shellcheck disable=SC2086
        "${FOLDSEEK}" cpdb ${TMP_PATH}/cluster_foldseek "${IN}_clu_clu" ${VERBOSITY}
    fi

    if [ -f "${IN}_clu_seq_ca.dbtype" ] && notExists "${IN}_clu_ca.dbtype"; then
        # shellcheck disable=SC2086
        "${FOLDSEEK}" createsubdb ${TMP_PATH}/cluster_foldseek "${IN}_clu_seq_ca" "${IN}_clu_ca" ${VERBOSITY}
    fi

    if notExists "${IN}_clu_aln.dbtype"; then
        # shellcheck disable=SC2086
        "${FOLDSEEK}" structurealign  "${IN}_clu_seq" "${IN}_clu_seq" "${TMP_PATH}/cluster_foldseek" "${IN}_clu_aln" -a -e 0.1 --sort-by-structure-bits 0 ${MERGECLU_PAR}
    fi

    if notExists "${IN}_clu_profile.dbtype"; then
        # shellcheck disable=SC2086
        "${MMSEQS}" result2profile  "${IN}_clu_seq" "${IN}_clu_seq" "${IN}_clu_aln" "${IN}_clu_profile" ${PROFILE_PAR}
    fi

    if notExists "${IN}_clu.dbtype"; then
        # shellcheck disable=SC2086
        "${FOLDSEEK}" profile2consensus "${IN}_clu_profile" "${IN}_clu" ${MERGECLU_PAR}
        if [ -n "$REMOVE_TMP" ]; then
            # shellcheck disable=SC2086
            "${FOLDSEEK}" rmdb "${IN}_clu_profile" ${VERBOSITY}
        fi
    fi

    if notExists "${IN}_clu_profile_ss.dbtype"; then
        # shellcheck disable=SC2086
        "${FOLDSEEK}" result2profile  "${IN}_clu_seq_ss" "${IN}_clu_seq_ss" "${IN}_clu_aln" "${IN}_clu_profile_ss" ${PROFILE_SS_PAR}
    fi

    if notExists "${IN}_clu_ss.dbtype"; then
        # shellcheck disable=SC2086
        "${FOLDSEEK}" profile2consensus "${IN}_clu_profile_ss" "${IN}_clu_ss" ${MERGECLU_PAR}
        if [ -n "$REMOVE_TMP" ]; then
            # shellcheck disable=SC2086
            "${FOLDSEEK}" rmdb "${IN}_clu_profile_ss" ${VERBOSITY}
            # shellcheck disable=SC2086
            "${FOLDSEEK}" rmdb "${IN}_clu_aln" ${VERBOSITY}
        fi
    fi
else
    if notExists "${TMP_PATH}/cluster_mmseqs.dbtype"; then
        # shellcheck disable=SC2086
        "${MMSEQS}" cluster "${IN}" "${TMP_PATH}/cluster_mmseqs" "${TMP_PATH}/cluster" ${CLUSTER_PAR} \
            || fail "cluster failed"
    fi

    if notExists "${TMP_PATH}/db_clu.dbtype"; then
        # shellcheck disable=SC2086
        "${MMSEQS}" createsubdb "${TMP_PATH}/cluster_mmseqs" "${IN}" "${TMP_PATH}/db_clu" ${VERBOSITY}\
            || fail "createsubdb failed"
    fi

    if notExists "${IN}_clu_rep_profile.dbtype"; then
        # shellcheck disable=SC2086
        "${MMSEQS}" result2profile "${TMP_PATH}/db_clu" "${IN}" "${TMP_PATH}/cluster_mmseqs" "${IN}_clu_rep_profile" ${THREADS_PAR}\
            || fail "result2profile failed"
    fi

    if notExists "${IN}_clu.dbtype"; then
        # shellcheck disable=SC2086
        "${MMSEQS}" profile2consensus "${IN}_clu_rep_profile" "${IN}_clu" ${CONSENSUS_PAR}\
            || fail "profile2consensus failed"
    fi

    if notExists "${IN}_clu_aln.dbtype"; then
        # shellcheck disable=SC2086
        "${MMSEQS}" align "${IN}" "${IN}" "${TMP_PATH}/cluster_mmseqs" "${IN}_clu_aln" -a ${THREADS_PAR} \
            || fail "align failed"
    fi
fi


if [ -n "${REMOVE_TMP}" ]; then
    echo "Remove temporary files"
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/cluster_mmseqs" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/cluster_foldseek" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/db_clu_2" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/db_clu_2_ss" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/db_clu_2_profile" ${VERBOSITY}
    rm -rf "${TMP_PATH}/cluster"
    rm -f "${TMP_PATH}/clusterdb.sh"
fi
