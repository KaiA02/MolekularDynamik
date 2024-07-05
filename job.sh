#!/bin/bash
#SBATCH -J disk
#SBATCH -o ./%x.%j.%N.out
#SBATCH -D ./build
#SBATCH --get-user-env
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
#SBATCH --mem=150MB
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --mail-user=kai_arenja@hotmail.de
#SBATCH --mail-type=end
#SBATCH --time=00:10:00

./MolSim ../input/eingabe-disk.xml 
