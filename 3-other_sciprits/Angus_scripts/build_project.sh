#!/usr/bin/env bash

if [ $# -ne 2 ]; then
    echo "Incorrect number of args specified"
    exit 1
fi

indir=$1
outdir=$2

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" > /dev/null && pwd )"

for f in ${indir}/*; do
    fname=$(basename $f ".fasta")
    fasta=$(basename $f)
    out=${outdir}/${fname}
    fastadir="${out}/fasta"
    chunkdir="${out}/chunks"
    tabuldir="${out}/tabular"
    kofamdir="${out}/kofam"

    TFLAGS="-clean -sformat pearson -frame 6"

    echo "*****Setting up $(basename $f)..."

    mkdir -p ${out}
    mkdir -p ${fastadir}
    mkdir -p ${chunkdir}
    mkdir -p ${tabuldir}
    mkdir -p ${kofamdir}

    echo "Copying file..."
    echo "cp $f ${fastadir}/"
    cp $f ${fastadir}
    echo "Done."

    echo "Translating fasta file..."
    transeq ${TFLAGS} -sequence ${fastadir}/${fasta} -outseq ${fastadir}/${fname}.faa
    echo "Done."

    echo "Chunking file..."
    ${DIR}/split_fasta.py ${fastadir}/${fname}.faa ${chunkdir}/
    echo "Done."

    echo "*****Done."
    echo ""
done
