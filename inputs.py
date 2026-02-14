GOAT_INP = """%PAL NPROCS 8 END
%MAXCORE 4000
%base "goat"
! XTB GOAT

*xyzfile 0 2 guess.xyz
"""

OPT_INP = """%PAL NPROCS 8 END
%MAXCORE 4000
%base "opt"
! WB97X-3C OPT

*xyzfile 0 2 goat.xyz
"""

FREQ_INP = """%PAL NPROCS 8 END
%MAXCORE 4000
%base "freq"
! M062X def2-TZVPP def2-TZVPP/c DEFGRID3 TightSCF SlowConv Opt NumFreq

%geom
 MaxIter 500
end

*xyzfile 0 2 opt.xyz
"""

CALC_INP = """%PAL NPROCS 8 END
%MAXCORE 4000
%base "calc"
! CCSD(T)-F12/RI cc-pVDZ-F12 cc-pVDZ-F12-CABS cc-pVTZ/c

*xyzfile 0 2 freq.xyz
"""

STATIONARY_SH = """#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=12:00:00
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

echo "Starting GOAT at $(date)"
$ORCA_EXE goat.inp > goat.log
cp "goat.log" "$SUBMIT_DIR/goat.log"

echo "Starting OPT at $(date)"
$ORCA_EXE opt.inp > opt.log
cp "opt.log" "$SUBMIT_DIR/opt.log"

echo "Starting FREQ at $(date)"
$ORCA_EXE freq.inp > freq.log
cp "freq.log" "$SUBMIT_DIR/freq.log"

echo "Starting CALC at $(date)"
$ORCA_EXE calc.inp > calc.log

if grep -q "ORCA TERMINATED NORMALLY" calc.log; then
    rsync -a "$SCRATCH_DIR/" "$SUBMIT_DIR/"

    echo "Pruning extra files in $SUBMIT_DIR..."
    cd "$SUBMIT_DIR"
    find . -type f ! -name 'freq.xyz' ! -name 'freq.hess' ! -name 'slurm-*.out' ! -name '*.log' -delete
        
else
    echo "ORCA error detected in calc.log."
    rsync -a "$SCRATCH_DIR/" "$SUBMIT_DIR/"
    find . -type f ! -name '*.xyz' ! -name '*.hess' ! -name 'slurm-*.out' ! -name '*.log' ! -name '*.inp' ! -name '*.sh' -delete

fi

# Disable the trap so the cleanup function doesn't rsync a second time on exit
trap - EXIT USR1
rm -rf "$SCRATCH_DIR"
"""

SCAN_INP = """%PAL NPROCS 8 END
%MaxCore 4000
%base "scan"
! XTB  ScanTS
%geom
 [SCAN] end
 TS_Mode {B [ATOMS]} end
 TS_Active_Atoms {[ATOMS]} end
end

* xyzfile 0 2 guess.xyz
"""

TRANS_FREQ_INP = """%PAL NPROCS 8 END
%MAXCORE 4000
%base "freq"
! M062X def2-TZVPP def2-TZVPP/c DEFGRID3 TightSCF SlowConv TightOpt OptTS Freq

%geom
 Recalc_Hess 5
 TS_Mode {B [ATOMS]} end
 TS_Active_Atoms {[ATOMS]} end
 MaxIter 500
end

*xyzfile 0 2 scan.xyz
"""

TRANSITION_SH = """#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=12:00:00
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

echo "Starting SCAN at $(date)"
$ORCA_EXE scan.inp > scan.log
cp "scan.log" "$SUBMIT_DIR/scan.log"

echo "Starting FREQ at $(date)"
$ORCA_EXE freq.inp > "$SUBMIT_DIR/freq.log"
cp "freq.log" "$SUBMIT_DIR/freq.log"

echo "Starting CALC at $(date)"
$ORCA_EXE calc.inp > calc.log

if grep -q "ORCA TERMINATED NORMALLY" calc.log; then
    rsync -a "$SCRATCH_DIR/" "$SUBMIT_DIR/"

    echo "Pruning extra files in $SUBMIT_DIR..."
    cd "$SUBMIT_DIR"
    find . -type f ! -name 'freq.xyz' ! -name 'freq.hess' ! -name '*trj*' ! -name 'slurm-*.out' ! -name '*.log' -delete
        
else
    echo "ORCA error detected in calc.log."
    rsync -a "$SCRATCH_DIR/" "$SUBMIT_DIR/"
    find . -type f ! -name '*.xyz' ! -name '*.hess' ! -name 'slurm-*.out' ! -name '*.log' ! -name '*.inp' ! -name '*.sh' -delete

fi

# Disable the trap so the cleanup function doesn't rsync a second time on exit
trap - EXIT USR1
rm -rf "$SCRATCH_DIR"
"""

"""

refined_file=$(ls scan.*.refined.xyz 2>/dev/null)
image_num=$(echo "$refined_file" | grep -oP "\d+")

prev=$(printf "%03d" $((10#$image_num - 1)))
cp scan.${prev}.xyz scan_prev.xyz
cp scan.${prev}.gbw scan_prev.gbw

curr=$(printf "%03d" $((10#$image_num)))
cp scan.${curr}.xyz scan_curr.xyz

next=$(printf "%03d" $((10#$image_num + 1)))
cp scan.${next}.xyz scan_next.xyz
cp scan.${next}.gbw scan_next.gbw

"""