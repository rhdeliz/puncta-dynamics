#!/bin/bash -l

#SBATCH -o ./job.out.%j
#SBATCH -e ./job.err.%j
#SBATCH -D ./
#SBATCH -J PPort2
#SBATCH --mail-type=END
#SBATCH --mail-user=slurm@deliz.me
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=72
#SBATCH --time=6:00:00

module purge
module load jdk/8.265
module load gcc/10 impi/2021.2
conda activate r_env

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=threads

cd /raven/u/deliz/dynamics_pipeline/supplemental_r_scripts/PhasePortraitAnalysis

Rscript --vanilla UserInput.R /raven/u/deliz/new_pipeline/Output4 /raven/u/deliz/dynamics_pipeline/supplemental_r_scripts/PhasePortraitAnalysis /raven/u/deliz/new_pipeline/Output4/PhasePortraitAnalysis 5 2 F 0
