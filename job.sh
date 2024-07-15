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
#SBATCH --time=30:00:00

./MolSim ../input/eingabe-Rayleigh-Taylor-3D.xml 
