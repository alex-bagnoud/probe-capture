#!/bin/bash
#
#SBATCH --job-name=prob-ko
#SBATCH --cpus-per-task=4
#SBATCH --mem=1GB
#SBATCH --output=/proj/hilts/Probe-Capture/logs/log.out
#SBATCH --err=/proj/hilts/Probe-Capture/logs/log.err
#SBATCH --mail-type=ALL

module load kofamscan

WRKDIR=/proj/hilts/Probe-Capture
MODELS=$WRKDIR/kofamscan/db/all_profiles.hmm

for faa in $WRKDIR/fasta/*.faa; do
    BNAME=$(basename ${faa} ".faa")
    OUTDIR=${WRKDIR}/output/${BNAME}
    TBL=${OUTDIR}/tabular/all_profiles
    KOUT=${OUTDIR}/kofam/${BNAME}.KO.txt
    echo "Running profile search..."
    echo "hmmsearch --cpu 4 -T0 --tblout=${TBL} ${MODELS} ${faa}"
    hmmsearch --cpu 4 -T0 --tblout=${TBL} ${MODELS} ${faa}
    echo "Done."

    echo "Running kofam scan annotation script..."
    echo "exec_annotation -o ${KOUT} -r -p ${MODELS} -f mapper --cpu 4 --tmp-dir ${OUTDIR} ${faa}"
    exec_annotation -o ${KOUT} -r -p ${MODELS} -f mapper --cpu 4 --tmp-dir ${OUTDIR} ${faa}
    echo "Done."
done
