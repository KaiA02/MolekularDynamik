#!/bin/bash
#SBATCH --job-name=my_molsim_job         # Name des Jobs
#SBATCH --output=my_molsim_job_out.txt    # Datei für Standardausgabe
#SBATCH --error=my_molsim_job_err.txt     # Datei für Fehlerausgabe
#SBATCH --time=04:00:00                  # Maximale Laufzeit (Stunden:Minuten:Sekunden)
#SBATCH --clusters=serial                # Cluster wählen (hier: serial)
#SBATCH --partition=serial_std           # Partition wählen (hier: serial_std)
#SBATCH --mem=4gb                        # Speicher pro Job
#SBATCH --cpus-per-task=1                # Anzahl der CPUs pro Job
#SBATCH --mail-type=END                  # E-Mail-Benachrichtigung bei Jobende
#SBATCH --mail-user=ge47jow@tum.de # Ihre E-Mail-Adresse für Benachrichtigungen

# Laden der benötigten Module
module load gcc
module load cmake

# Navigieren zum Verzeichnis, in dem das Programm kompiliert wird
cd MolekularDynamik/
mkdir build
cd build

# Kompilieren des Programms mit gprof
cmake ..
make
gcc -pg -o MolSim ../input/eingabe-equilibrium.xml

# Ausführen des Programms mit einem XML-Eingabedateiargument
./MolSim ../input/eingabe-equilibrium.xml

# Optional: Speichern der gmon.out-Datei für das Profiling
cp gmon.out /MolekularDynamik/build/gmon.out
