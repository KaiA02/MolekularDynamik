#!/bin/bash
#SBATCH -J "name"
#SBATCH -o ./%x.%j.%N.out
#SBATCH -D ./MolekularDynamik/build
#SBATCH --get-user-env
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
#SBATCH --mem=200MB
#SBATCH --get-user-env
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --mail-user=kai_arenja@hotmail.de
#SBATCH --mail-type=end
#SBATCH --time=00:30:00

./MolSim ../input/"name".xml 
