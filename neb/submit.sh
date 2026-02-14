#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=06:00:00
#SBATCH --gres=lscratch:100
#SBATCH --signal=B:USR1@60

set -uo pipefail

SUBMIT_DIR=$(pwd)
SCRATCH_DIR="/lscratch/${USER}/orca_${SLURM_JOB_ID}"

mkdir -p "$SCRATCH_DIR"
rsync -a "$SUBMIT_DIR/" "$SCRATCH_DIR/"
cd "$SCRATCH_DIR"

cleanup() {
    echo "Signal/Error caught. Syncing progress to $SUBMIT_DIR..."
    rsync -a --exclude='*tmp*' "$SCRATCH_DIR/" "$SUBMIT_DIR/"
    rm -rf "$SCRATCH_DIR"
}
trap cleanup EXIT USR1

module load ORCA/6.1
ORCA_EXE=$(which orca)

echo "Starting NEB at $(date)"
$ORCA_EXE neb.inp > "$SUBMIT_DIR/neb.log"

echo "Starting OPT at $(date)"
$ORCA_EXE opt.inp > "$SUBMIT_DIR/opt.log"