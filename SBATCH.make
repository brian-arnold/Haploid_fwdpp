#!/bin/bash
#SBATCH -J make
#SBATCH -o out
#SBATCH -e err
#SBATCH -p test
#SBATCH -n 1
#SBATCH -t 50
#SBATCH --mem=3000

source new-modules.sh
module load gcc/7.1.0-fasrc01

make
