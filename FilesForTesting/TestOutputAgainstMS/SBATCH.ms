#!/bin/bash
#SBATCH -J perl
#SBATCH -o out
#SBATCH -e err
#SBATCH -p test
#SBATCH -n 1
#SBATCH -t 400
#SBATCH --mem=3000

./msdir/ms 100 100 -t 1000 -r 0 100000 -c 1000 300 > ms_output.txt
