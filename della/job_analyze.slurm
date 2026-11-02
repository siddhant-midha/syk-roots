#!/bin/bash
#SBATCH --job-name=syk-roots-an
#SBATCH --account=dabanin
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:30:00
#SBATCH --output=logs/%x_%A_%a.out
# Zeros of the annealed partition function for one (n, k): line SLURM_ARRAY_TASK_ID+1 of $PAIRS = "n k".
# Reads out/n{n}_k{k}/*.npz, writes results/n{n}_k{k}.{json,npz}.
set -euo pipefail
cd /scratch/gpfs/SARANGG/sm7456/jobs/syk-roots
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
read -r N K < <(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$PAIRS")
mkdir -p results
$HOME/miniconda3/envs/ldpc-dev/bin/python -u code/analyze.py out/n${N}_k${K}/*.npz --out results/n${N}_k${K} ${EXTRA:-}
