#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

abspath() {
    if [ -d "$1" ]; then
        (cd "$1"; pwd)
    elif [ -f "$1" ]; then
        if [ -z "${1##*/*}" ]; then
            echo "$(cd "${1%/*}"; pwd)/${1##*/}"
        else
            echo "$(pwd)/$1"
        fi
    elif [ -d "$(dirname "$1")" ]; then
        echo "$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"
    fi
}

ARR=""
push_back() {
    # shellcheck disable=SC1003
    CURR="$(printf '%s' "$1" | awk '{ gsub(/'\''/, "'\''\\'\'''\''"); print; }')"
    if [ -z "$ARR" ]; then
        ARR=''\'$CURR\'''
    else
        ARR=$ARR' '\'$CURR\'''
    fi
}

hasCommand () {
    command -v "$1" >/dev/null 2>&1 || { echo "Please make sure that $1 is in \$PATH."; exit 1; }
}

hasCommand awk
hasCommand join
hasCommand sort

[ -z "$MMSEQS" ] && echo "Please set the environment variable \$MMSEQS to your MMSEQS binary." && exit 1;

export MMSEQS_FORCE_MERGE=1

OUTDB="$(abspath "${OUTDB}")"

#check if already created db
if notExists "${1}.dbtype"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" createdb "$@" "${TMP_PATH}/seqDB" ${CREATEDB_PAR} \
        || fail "createdb failed"
        
else
    echo "Input DB already exists. Generating associated metadata."

	if notExists "${TMP_PATH}/seqDB.index"; then
	    # shellcheck disable=SC2086
	    "${MMSEQS}" cpdb "$1" "${TMP_PATH}/seqDB" ${VERBOSITY} \
		|| fail "cpdb failed"
	fi

	if notExists "${TMP_PATH}/seqDB_h.index"; then
	    # shellcheck disable=SC2086
	    "${MMSEQS}" cpdb "$1_h" "${TMP_PATH}/seqDB_h" ${VERBOSITY} \
		|| fail "cpdb failed"
	fi

	if notExists "${OUTDB}_ss.index"; then
	    # shellcheck disable=SC2086
	    "${MMSEQS}" cpdb "$1_ss" "${OUTDB}_ss" ${VERBOSITY} \
		|| fail "cpdb failed"
	fi

    touch "${TMP_PATH}/seqDB.external"
fi

if [ "$("${MMSEQS}" dbtype "${TMP_PATH}/seqDB")" = "Nucleotide" ]; then

    echo "Input DB type is Nucleotide."

    [ -z "$GFFDIR" ] && fail "No GFF directory file is given. Please provide a valid path to GFF directory file with the --gff-dir parameter."
        
    GFFDIR="$(abspath "${GFFDIR}")"

    if notExists "${GFFDIR}"; then 
        fail "Cannot find GFF directory file. Please provide a valid path to GFF directory file with the --gff-dir parameter."
    fi

    if notExists "${OUTDB}_nucl.index"; then
        #read line by line and push back to array
        while IFS= read -r line; do
            push_back "$line"
        done < "${GFFDIR}"
        eval "set -- $ARR"
        # shellcheck disable=SC2086
        "${MMSEQS}" gff2db "${@}" "${TMP_PATH}/seqDB" "${OUTDB}_nucl" ${GFF2DB_PAR} \
            || fail "gff2db failed"
    fi

    if notExists "${OUTDB}.index"; then
        # shellcheck disable=SC2086
        "${MMSEQS}" translatenucs "${OUTDB}_nucl" "${OUTDB}" ${TRANSLATENUCS_PAR} \
            || fail "translatenucs failed"
    fi

elif [ "$("${MMSEQS}" dbtype "${TMP_PATH}/seqDB")" = "Aminoacid" ]; then 

    echo "Input DB type is Aminoacid."

    if notExists "${TMP_PATH}/seqDB_h_pref"; then
        # shellcheck disable=SC2086
        "${MMSEQS}" prefixid "${TMP_PATH}/seqDB_h" "${TMP_PATH}/seqDB_h_pref" --tsv ${THREADS_PAR}\
            || fail "prefixid failed"
    fi

    if notExists "${TMP_PATH}/seqDB.lookup.tmp"; then
        #remove whitespaces -> split prodigal header by "#" and print relavant columns -> swap start and end if gene on minus strand -> merge to one column of header -> sort by seqid -> replace the seq name column in lookup -> recount pos idx by setid
        awk '{ gsub(/ /,""); print }' "${TMP_PATH}/seqDB_h_pref" \
        |awk -F '[\t#]' 'NF{NF-=1};1' OFS='\t' \
        |awk -F'\t' '$5=="-1" { temp = $4; $4 = $3; $3 = temp } 1' OFS='\t' \
        |sort -k1,1n > "${TMP_PATH}/seqDB_h_pref.tmp"
    
        if notExists "${TMP_PATH}/seqDB.external"; then
            join -t "$(printf '\t')" -o '1.1 2.2 2.3 2.4 1.3' "${TMP_PATH}/seqDB.lookup" "${TMP_PATH}/seqDB_h_pref.tmp" \
            |awk -F '[\t]' '{ if (setid == $NF) { counter++ } else { counter = 1; setid = $NF }; print $1"\t"$2"_"counter-1"_"$3"_"$4"\t"$NF }' \
            > "${TMP_PATH}/seqDB.lookup.tmp"
        else
            join -t "$(printf '\t')" -o '1.1 2.2 2.3 2.4 1.3' "${TMP_PATH}/seqDB.lookup" "${TMP_PATH}/seqDB_h_pref.tmp" \
            |sort -k2,2 -V\
            |awk -F '[\t]' '{ if (setid == $NF) { counter++ } else { counter = 1; setid = $NF }; print $1"\t"$2"_"counter-1"_"$3"_"$4"\t"$NF }' \
            |sort -k1,1n\
            > "${TMP_PATH}/seqDB.lookup.tmp"
        fi
    fi
    
    mv -f -- "${TMP_PATH}/seqDB.lookup.tmp" "${TMP_PATH}/seqDB.lookup"
    rm "${TMP_PATH}/seqDB_h_pref"
    rm "${TMP_PATH}/seqDB_h_pref.tmp"

	if notExists "${OUTDB}.index"; then
	    # shellcheck disable=SC2086
	    "${MMSEQS}" mvdb "${TMP_PATH}/seqDB" "${OUTDB}" ${VERBOSITY} \
		|| fail "mvdb failed"
	fi

	if notExists "${OUTDB}_h.index"; then
	    # shellcheck disable=SC2086
	    "${MMSEQS}" mvdb "${TMP_PATH}/seqDB_h" "${OUTDB}_h" ${VERBOSITY} \
		|| fail "mvdb failed"
	fi
    
    # mv only if ${TMP_PATH}/seqDB.source exists
    if [ -f "${TMP_PATH}/seqDB.source" ]; then
        mv -f "${TMP_PATH}/seqDB.source" "${OUTDB}.source"
    fi
else
    fail "Input DB has the wrong type. Allowed input: Nucleotide, Aminoacid";
fi

#TODO: run aa2foldseek instead to generate mapping files
if notExists "${OUTDB}_member_to_set.index"; then
    awk '{ print $1"\t"$3; }' "${OUTDB}.lookup" | sort -k1,1n -k2,2n > "${OUTDB}_member_to_set.tsv"
    "${MMSEQS}" tsv2db "${OUTDB}_member_to_set.tsv" "${OUTDB}_member_to_set" --output-dbtype 5 \
        || fail "tsv2db failed"
fi

if notExists "${OUTDB}_set_to_member.index"; then
    awk '{ print $3"\t"$1; }' "${OUTDB}.lookup" | sort -k1,1n -k2,2n > "${OUTDB}_set_to_member.tsv"
    "${MMSEQS}" tsv2db "${OUTDB}_set_to_member.tsv" "${OUTDB}_set_to_member" --output-dbtype 5 \
        || fail "tsv2db failed"
fi

if notExists "${OUTDB}_set_size.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" result2stats "${OUTDB}" "${OUTDB}" "${OUTDB}_set_to_member" "${OUTDB}_set_size" ${RESULT2STATS_PAR} \
        || fail "result2stats failed"
fi

#add check for lookup format

if [ -n "${REMOVE_TMP}" ]; then
    echo "Remove temporary files"
    rm -f "${TMP_PATH}/createsetdb.sh"
fi
