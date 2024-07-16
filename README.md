# MolSim
=========

"GroupB"

## Program Execution:

1. **mkdir build && cd $_**
2. **cmake ..**
3. **make**
4. **ctest**
5. **./MolSim XML-input-file** (e.g **./MolSim ../input/eingabe-sonne.xml**)
6. **make doc**
7. **gprof MolSim gmon.out > analysis.txt** (for profiling. Uncomment this in cmakeList.txt: set(CMAKE_CXX_FLAGS "
   ${CMAKE_CXX_FLAGS} -pg") )

## XML input-file consists of:

1. Start_time as Double
2. End_Time as Double
3. Time_Step as Double
4. SmoothLJ as boolean (true for smoothLJ and false for normal LJ)
4. Output type: "xyz" or "vtk"
5. Input type: "cube" or "sonne"
6. ... see xsd file
7. (optional) performance measurement: true or false
8. (optional) log level: "off", "error", "warn", "info", "debug", "trace"

## Quick start guide:

1. adapt inputFile to your need
2. follow Steps described in Program Execution
3. open Paraview
4. open generated outputFiles in Paraview and add a Glyph to the Data
5. run Simulation

## Thermostat:

- First set the **thermostatON** option in xml on **YES**
- In the beginning the system is set to the initial temperature.
- **n_thermostat** means after how many steps the thermostat is applied.
- Reaching and holding the target temperature is the aim, therefore we have the **delta_temp** to say what the absolute
  change of temperature allowed is (for the steps, declared at **n_thermostat**)
- If no **delta_temp** is given, then we see declare the **delta_temp** as 0 in xml file, but see it as infinity
- If **n_thermostat** = 0, then the thermostat is applied directly.

## Membrane:

- **n1**, **n2**, **n3** are count of Particles in Membrane in x, y and z direction
- **distance** is distance between Particles in Membrane
- Particles in Membrane are always of **type 0** so they can be identified more easy.

## using Checkpoints (saving the state after simulation)

- after every simulation a new file Output.txt is created in the Input directory
- in case you want to build your new Simulation on top of another simulation, you can also input this output.txt file at
  the end of your console input.
- Step by step:
    - first call first simulation: **./MolSim ../input/eingabe-equilibrium.xml**
      (here you generate the input/output.txt file)
    - then call second simulation on top of first: **./MolSim ../input/eingabe-sphere.xml ../input/output.txt**

## Linux cluster

1. Login on "ssh -CX di57yey@lxlogin1.lrz.de"
2. Steps

- projekt clonen:  
  rm -rf MolekularDynamik  
  #main branch (falls spezifische branch gewünscht ist "-b branchName" anhängen)  
  git clone https://github.com/KaiA02/MolekularDynamik

- in projekt module laden und cmake/make ausführen  
  cd MolekularDynamik  
  module load xerces-c  
  module load cmake/3.21.4  
  module load gcc/11  
  mkdir build && cd build  
  cmake ..  
  make -j 5

- script ausführen  
  salloc --partition=cm2_inter -t 10   (10min is time (max 2h))   
  ./MolSim ../input/"name".xml

- oder so script ausführen  
  cd ..  
  sbatch job.sh


- status prüfen  
  squeue --clusters=cm2_tiny --start

- output ordner auf pc laden (im normalen terminal machen)  
  #zielquelle individuell einstellen     
  #-r für ganzen ordner  
  scp -r di57yey@lxlogin1.lrz.de:MolekularDynamik/output "zielquelle" (z.B. /Users/kaiarenja/Desktop/)

# create new c++ classes from xsd file:

xsdcxx cxx-tree --std c++11 /mnt/c/Users/joshu/CLionProjects/MolekularDynamik/src/inputReader/simulation.xsd
xsdcxx cxx-tree --std c++11 /home/kaiarenja/CLionProjects/MolekularDynamik/src/inputReader/simulation.xsd  
xsdcxx cxx-tree --std c++11 /mnt/c/Users/jh/CLionProjects/Molekulardynamik/Assignment5/src/inputReader/simulation.xsd

## explaination for realocateParticles in LCParticleContainer


## Parallel Strategies

# Strategy 1: **OpenMP for Loops**
  This strategy uses OpenMP to parallelize loops with #pragma omp parallel for. It is ideal for evenly distributed workloads, allowing multiple particles to be processed simultaneously, which enhances performance for large datasets.
  Key Features:
  -Simple implementation.
  -Efficient for evenly distributed tasks.
# Strategy 2: **OpenMP Task Parallelism**
  This strategy employs OpenMP's task-based parallelism using #pragma omp task. It provides better load balancing and resource utilization, particularly for workloads that vary significantly.
  Key Features:
  -Handles uneven workloads efficiently.
  -Uses dynamic task scheduling.
  Set the ParallelStrategy variable in the Calculations class to switch between these strategies based on your simulation needs.
# HowToDo:
  to use the Parallel Strategy you specify your Choice in the input in the ParallelStrategy Variable. 0 => no Strategy, 1 => Strategy1, 2 => Strategy2. After this you call for example: OMP_NUM_THREADS=14 ./MolSim ../input/eingabe-Rayleigh-Taylor-3D.xml




