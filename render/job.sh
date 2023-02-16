#!/bin/bash -l
#SBATCH -J rose
#SBATCH -o /path/to/rose/states/rose.out
#SBATCH -e /path/to/rose/states/rose.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rose-slurm@unkaktus.art
#SBATCH --ntasks=21
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=7
#SBATCH --time=48:00:00

source /home/SPACK2023/share/spack/setup-env.sh
module load apptainer-1.0.3-gcc-12.2.0-aojy6ca

cd /path/to/rose
srun -n 21 ./render.sh states/state.pvsm