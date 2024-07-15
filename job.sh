#!/bin/bash
#SBATCH -J Task3
#SBATCH -o ./%x.%j.%N.out
#SBATCH -D ./build
#SBATCH --get-user-env
#SBATCH --clusters=cm2
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --mail-user=kai_arenja@hotmail.de
#SBATCH --mail-type=end
#SBATCH --time=9:00:00

OMP_NUM_THREADS=56 ./MolSim ../input/eingabe-Rayleigh-Taylor-3D.xml 
