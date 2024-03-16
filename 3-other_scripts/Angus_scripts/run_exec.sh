#!/bin/bash
#
#SBATCH --job-name=prob-ko
#SBATCH --cpus-per-task=4
#SBATCH --mem=1GB
#SBATCH --output=/proj/hilts/Probe-Capture/logs/log_%A-%a.out
#SBATCH --err=/proj/hilts/Probe-Capture/logs/log_%A-%a.err
#SBATCH --mail-type=ALL

module load kofamscan

CPU=6

WRKDIR=/proj/hilts/Probe-Capture
MODELS=$WRKDIR/kofamscan/db/all_profiles.hmm
FAAPREFIX=$1
FAA=${WRKDIR}/genomes/${FAAPREFIX}/${FAAPREFIX}.faa
TBL=${WRKDIR}/genomes/${FAAPREFIX}/tabular/all_profiles
DTBL=${WRKDIR}/genomes/${FAAPREFIX}/all_profiles.domains.tsv
OUT=${WRKDIR}/genomes/${FAAPREFIX}/${FAAPREFIX}.KO.txt
TMPDIR=${WRKDIR}/genomes/tmp/
KOLIST=${WRKDIR}/kofamscan/db/ko_list
KFLAGS="-r -f detail --cpu ${CPU} --tmp-dir ${TMPDIR}"


# Run the HMM
echo "Running HMMs..."
hmmsearch --cpu ${CPU} --tblout=${TBL} --domtblout=${DTBL} ${MODELS} ${FAA} >/dev/null
echo "Done."

echo "Creating directory sturcture..."
mkdir -p ${TMPDIR}
mkdir -p ${TMPDIR}/tabular
mkdir -p ${WRKDIR}/genomes/${FAAPREFIX}/results
# Copy the table to the tmp dir so that the naming is correct
cp ${TBL} ${TMPDIR}/tabular/all_profiles

echo "Running exec_annotation..."
echo "/usr/bin/time exec_annotation ${KFLAGS} -o ${OUT} -p ${MODELS} -k ${KOLIST} ${FAA}"
/usr/bin/time exec_annotation ${KFLAGS} -o ${OUT} -p ${MODELS} -k ${KOLIST} ${FAA}

# Clean up
echo "Cleaning up..."
rm ${TMPDIR}/tabular/all_profiles

echo "Done."
