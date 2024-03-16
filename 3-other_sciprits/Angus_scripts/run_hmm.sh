#!/bin/bash
#
#SBATCH --job-name=prob-ko
#SBATCH --cpus-per-task=3
#SBATCH --mem=1GB
#SBATCH --output=/proj/hilts/Probe-Capture/logs/log_%A-%a.out
#SBATCH --err=/proj/hilts/Probe-Capture/logs/log_%A-%a.err
#SBATCH --mail-type=ALL

module load kofamscan

WRKDIR=/proj/hilts/Probe-Capture
MODELS=$WRKDIR/kofamscan/db/all_profiles.hmm
FAAPREFIX=$1
FAA=${WRKDIR}/samples/${FAAPREFIX}/chunks/${FAAPREFIX}_${SLURM_ARRAY_TASK_ID}.faa
TBL=${WRKDIR}/samples/${FAAPREFIX}/tabular/all_prokaryotes_${SLURM_ARRAY_TASK_ID}
DTBL=${WRKDIR}/samples/${FAAPREFIX}/domains/all_prokaryotes_${SLURM_ARRAY_TASK_ID}

hmmsearch --cpu 4 -T0 --tblout=${TBL} --domtblout=${DTBL} ${MODELS} ${FAA}

