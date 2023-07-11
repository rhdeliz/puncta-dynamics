#!/bin/bash -l

#SBATCH -o ./job.out.%j
#SBATCH -e ./job.err.%j
#SBATCH -D ./
#SBATCH -J PhasePortrait
#SBATCH --mail-type=ALL
#SBATCH --mail-user=slurm@deliz.me
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=72
#SBATCH --time=5:00:00

module purge
module load jdk/8.265 gcc/10 impi/2021.2 fftw-mpi R/4.0.2
echo 'modules loaded'

conda activate dynamics_pipeline
echo 'conda activated'

export OMP_NUM_THREDS=144

cd /raven/u/deliz/dynamics_pipeline/supplemental_r_scripts/PhasePortraitAnalysis
Rscript --vanilla --verbose UserInput.R '/raven/u/deliz/new_pipeline /raven/u/deliz/dynamics_pipeline/supplemental_r_scripts/PhasePortraitAnalysis /raven/u/deliz/new_pipeline/phase_portraits 5 2 199 399 0.01 0.95 0'

sleep 10
