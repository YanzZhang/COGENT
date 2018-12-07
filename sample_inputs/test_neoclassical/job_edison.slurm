#!/bin/bash -l
#SBATCH -q debug
#SBATCH -N 2
#SBATCH -t 00:20:00
#SBATCH -L SCRATCH

srun -n64 /global/homes/y/yaz148/test_collision/COGENT/exec/cogent.Linux.64.CC.ftn.OPT.MPI.ex neoclassical.in >output.txt
