#!/bin/sh -e
# Iterative cluster search workflow script
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
[ "$#" -ne 4 ] && echo "Please provide <queryDB> <targetDB> <outDB> <tmp>" && exit 1;
# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[ ! -f "$2.dbtype" ] && echo "$2.dbtype not found!" && exit 1;
[   -f "$3.dbtype" ] && echo "$3.dbtype exists already!" && exit 1;
[ ! -d "$4" ] && echo "tmp directory $4 not found!" && mkdir -p "$4";

QUERY="$1"
#for accessing the metadata
QUERY_STEP_0="$1"
TARGET="$2"
TMP_PATH="$4"

STEP=0
# processing
[ -z "$NUM_IT" ] && NUM_IT=3;
while [ $STEP -lt $NUM_IT ]; do
    # call prefilter module
    if notExists "$TMP_PATH/pref_tmp_${STEP}.done"; then
        PARAM="PREFILTER_PAR_$STEP"
        eval TMP="\$$PARAM"
        if [ $STEP -eq 0 ]; then
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" prefilter "${QUERY}" "${TARGET}" "$TMP_PATH/pref_$STEP" ${TMP} \
                || fail "Prefilter failed"
        else
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" prefilter "${QUERY}" "${TARGET}" "$TMP_PATH/pref_tmp_$STEP" ${TMP} \
                || fail "Prefilter failed"
        fi
        touch "$TMP_PATH/pref_tmp_${STEP}.done"
    fi

    if [ $STEP -ge 1 ]; then
        if notExists "$TMP_PATH/pref_$STEP.done"; then
            STEPONE=$((STEP-1))
            # shellcheck disable=SC2086
            "$MMSEQS" subtractdbs "$TMP_PATH/pref_tmp_$STEP" "$TMP_PATH/aln_$STEPONE" "$TMP_PATH/pref_$STEP" $SUBSTRACT_PAR \
                || fail "Substractdb failed"
            "$MMSEQS" rmdb "$TMP_PATH/pref_tmp_$STEP"
        fi
        touch "$TMP_PATH/pref_$STEP.done"
    fi

	# call alignment module
	if notExists "$TMP_PATH/aln_tmp_$STEP.done"; then
	    PARAM="ALIGNMENT_PAR_$STEP"
        eval TMP="\$$PARAM"
        if [ $STEP -eq 0 ]; then
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" "${ALIGN_MODULE}" "${QUERY}" "${TARGET}" "$TMP_PATH/pref_$STEP" "$TMP_PATH/aln_$STEP" ${TMP} \
                || fail "Alignment failed"
        else
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" "${ALIGN_MODULE}" "${QUERY}" "${TARGET}" "$TMP_PATH/pref_$STEP" "$TMP_PATH/aln_tmp_$STEP" ${TMP} \
                || fail "Alignment failed"
        fi
        touch "$TMP_PATH/aln_tmp_$STEP.done"
    fi



    if [ $STEP -gt 0 ]; then
        if notExists "$TMP_PATH/aln_$STEP.done"; then
            STEPONE=$((STEP-1))
            # shellcheck disable=SC2086
            "$MMSEQS" mergedbs "${QUERY_STEP_0}" "$TMP_PATH/aln_$STEP" "$TMP_PATH/aln_$STEPONE" "$TMP_PATH/aln_tmp_$STEP" \
                || fail "mergedbs failed"
            "$MMSEQS" rmdb "$TMP_PATH/aln_tmp_$STEP"
            "$MMSEQS" rmdb "$TMP_PATH/aln_$STEPONE"
            touch "$TMP_PATH/aln_$STEP.done"
        fi
    fi

    if [ $STEP -eq $((NUM_IT  - 1)) ]; then
        # clustersearch pipeline in the last iteration
        if notExists "$TMP_PATH/cluster_aln.done"; then

            if notExists "${TMP_PATH}/result_prefixed.index"; then
                # shellcheck disable=SC2086
                "${MMSEQS}" prefixid  "${TMP_PATH}/aln_$STEP" "${TMP_PATH}/result_prefixed" ${THREADS_PAR} \
                    || fail "prefixid failed"
            fi

            if notExists "${TMP_PATH}/aggregate.index"; then
            # aggregation: take for each target set the best hit
            # shellcheck disable=SC2086
                "${MMSEQS}" besthitbyset "${QUERY_STEP_0}" "${TARGET}" "${TMP_PATH}/result_prefixed" "${TMP_PATH}/aggregate" ${BESTHITBYSET_PAR} \
                    || fail "aggregate best hit failed"
            fi

            if notExists "${TMP_PATH}/aggregate_merged.index"; then
                # shellcheck disable=SC2086
                "${MMSEQS}" mergeresultsbyset "${QUERY_STEP_0}_set_to_member" "${TMP_PATH}/aggregate" "${TMP_PATH}/aggregate_merged" ${THREADS_PAR} \
                    || fail "mergesetresults failed"
            fi


            if notExists "${TMP_PATH}/matches.index"; then
                # shellcheck disable=SC2086
                "${MMSEQS}" combinehits "${QUERY_STEP_0}" "${TARGET}" "${TMP_PATH}/aggregate_merged" "${TMP_PATH}/matches" "${TMP_PATH}" ${combinehits_PAR} \
                    || fail "combinehits failed"
            fi

            #db-output set to false to not print the null bytes
            if notExists "${TMP_PATH}/cluster.index"; then
                # shellcheck disable=SC2086
                "${MMSEQS}" clusterhits "${QUERY_STEP_0}" "${TARGET}" "${TMP_PATH}/matches" "${TMP_PATH}/cluster" ${CLUSTERHITS_PAR} \
                    || fail "clusterhits failed"
            fi

            #TODO:sort by qid to group by qid?
            if notExists "${TMP_PATH}/cluster_sorted.index"; then
                sort -k1,1n "${TMP_PATH}/cluster" > "${TMP_PATH}/cluster_sorted" \
                || fail "sort failed"
            fi

            #TODO:--output-dbtype 5 which is the alignment dbtype
            # shellcheck disable=SC2086
            "${MMSEQS}" tsv2db "${TMP_PATH}/cluster_sorted" "$3" --output-dbtype 5 ${VERBOSITY} \
                || fail "tsv2db failed"

            #TODO: parameterize create profile or show original clustering
            PARAM="PROFILE_PAR"
            eval TMP="\$$PARAM"
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" result2profile "${QUERY_STEP_0}" "${TARGET}" "$3" "$3_profile" ${TMP} \
                || fail "create profile failed"
            # # shellcheck disable=SC2086
            # "$MMSEQS" cpdb "${QUERY_STEP_0}_h" "$3_profile_h"
            # # shellcheck disable=SC2086
            # "$MMSEQS" cpdb "${QUERY_STEP_0}_nucl" "$3_profile_nucl"
            # # shellcheck disable=SC2086
            # "$MMSEQS" cpdb "${QUERY_STEP_0}_set_size" "$3_profile_set_size"
            # # shellcheck disable=SC2086
            # "$MMSEQS" cpdb "${QUERY_STEP_0}_member_to_set" "$3_profile_member_to_set"
            # # shellcheck disable=SC2086
            # "$MMSEQS" cpdb "${QUERY_STEP_0}_set_to_member" "$3_profile_set_to_member"
            # # shellcheck disable=SC2086
            # "$MMSEQS" cpdb "$3" "$3_profile_aln"


            if [ -n "${REMOVE_TMP}" ]; then
                # shellcheck disable=SC2086
                "$MMSEQS" rmdb "${TMP_PATH}/aggregate" ${VERBOSITY}
                # shellcheck disable=SC2086
                "$MMSEQS" rmdb "${TMP_PATH}/aggregate_merged" ${VERBOSITY}
                # shellcheck disable=SC2086
                "$MMSEQS" rmdb "${TMP_PATH}/cEval" ${VERBOSITY}
                # shellcheck disable=SC2086
                "$MMSEQS" rmdb "${TMP_PATH}/match" ${VERBOSITY}
                # shellcheck disable=SC2086
                "$MMSEQS" rmdb "${TMP_PATH}/aggregate_prefixed" ${VERBOSITY}
                # shellcheck disable=SC2086
                "$MMSEQS" rmdb "${TMP_PATH}/aggregate_prefixed_merged" ${VERBOSITY}
                # shellcheck disable=SC2086
                "$MMSEQS" rmdb "${TMP_PATH}/matches" ${VERBOSITY}
                # shellcheck disable=SC2086
                "$MMSEQS" rmdb "${TMP_PATH}/cluster" ${VERBOSITY}
            fi

        touch "$TMP_PATH/cluster_aln.done"
        fi
    fi

# create profiles
    if [ $STEP -ne $((NUM_IT  - 1)) ]; then
        if notExists "$TMP_PATH/profile_$STEP.dbtype"; then
            PARAM="PROFILE_PAR_$STEP"
            eval TMP="\$$PARAM"
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" result2profile "${QUERY}" "${TARGET}" "$TMP_PATH/aln_$STEP" "$TMP_PATH/profile_$STEP" ${TMP} \
                || fail "create profile failed"
        fi
    fi
	QUERY="$TMP_PATH/profile_$STEP"
	STEP=$((STEP+1))
done

if [ -n "$REMOVE_TMP" ]; then
    STEP=0
    while [ "$STEP" -lt "$NUM_IT" ]; do
        if [ $STEP -gt 0 ]; then
            rm -f -- "$TMP_PATH/aln_$STEP.done" "$TMP_PATH/pref_$STEP.done"
        fi
        rm -f -- "$TMP_PATH/aln_tmp_$STEP.done" "$TMP_PATH/pref_tmp_${STEP}.done"
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/pref_$STEP" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/aln_$STEP" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/profile_$STEP" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/profile_${STEP}_h" ${VERBOSITY}
        STEP=$((STEP+1))
    done
    rm -f "$TMP_PATH/iterativevlustersearch.sh"
fi
