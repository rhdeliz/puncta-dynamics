#!/bin/bash -l

#SBATCH -o ./job.out.%j
#SBATCH -e ./job.err.%j
#SBATCH -D ./
#SBATCH -J batch_3
#SBATCH --mail-type=ALL
#SBATCH --mail-user=slurm@deliz.me
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=72
#SBATCH --time=3:00:00

# Load all needed packages
module purge
module load jdk/8.265 gcc/10 impi/2021.2 fftw-mpi R/4.0.2
echo 'modules loaded'
conda activate dynamics_pipeline
echo 'conda activated'

# Specify parameters
## Path of parameters table
path=$'/raven/u/deliz/new_pipeline/pending_processing/batch_3/Input/parameter_tables'

## Scripts folder
cd /raven/u/deliz/dynamics_pipeline

## Cores for parallel processing in R
export OMP_NUM_THREDS=144

# Run scripts
## Python scripts
python mission_control.py $path 12

## Run R Scripts
Rscript --vanilla --verbose r_scripts/extract_intensity.R $path
Rscript --vanilla --verbose r_scripts/colocalization.R $path
Rscript --vanilla --verbose r_scripts/compile_tables.R $path
# Rscript --vanilla --verbose r_scripts/compress_everything.R $path

sleep 10
