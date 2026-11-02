#!/bin/bash
#SBATCH --job-name=syk-roots
#SBATCH --account=dabanin
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:10:00
#SBATCH --output=logs/%x_%A_%a.out
# Annealed SYK partition function by ED (repo siddhant-midha/syk-roots, code/ = copy of the repo).
# Task = line SLURM_ARRAY_TASK_ID+1 of $TASKS: "n k seed time_limit_seconds". Each task samples disorder until its
# time limit and writes out/n{n}_k{k}/s{seed}.npz (checkpointed every 10 min). Submit e.g.
#   sbatch --array=0-383 --export=ALL,TASKS=tasks_small.txt job.slurm
#   sbatch --array=0-15 --cpus-per-task=4 --mem=8G --time=03:00:00 --export=ALL,TASKS=tasks_n20.txt job.slurm
set -euo pipefail
cd /scratch/gpfs/SARANGG/sm7456/jobs/syk-roots
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
read -r N K SEED TL < <(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$TASKS")
mkdir -p out/n${N}_k${K}
/usr/bin/time -v $HOME/miniconda3/envs/ldpc-dev/bin/python -u code/sample.py --n $N --k $K --seed $SEED \
    --samples 1000000000000 --time-limit $TL --out out/n${N}_k${K}/s${SEED}.npz 2>&1 | grep -vE "^\s+(Page|Volun|Invol|Swaps|File|Socket|Signals|Average|Exit|Minor|Major|Percent|System|User|Command)"
