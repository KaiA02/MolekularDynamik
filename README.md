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

## XML input-file consists of:

1. Start_time as Double
2. End_Time as Double
3. Time_Step as Double
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


## using Checkpoints (saving the state after simulation)
- after every simulation a new file Output.txt is created in the Input directory
- in case you want to build your new Simulation on top of another simulation, you can also input this output.txt file at
the end of your console input.
- Step by step: 
  - first call first simulation: **./MolSim ../input/eingabe-equilibrium.xml** 
            (here you generate the input/output.txt file)
  - then call second simulation on top of first: **./MolSim ../input/eingabe-sphere.xml ../input/output.txt**


# create new c++ classes from xsd file:
xsdcxx cxx-tree --std c++11 /mnt/c/Users/joshu/CLionProjects/MolekularDynamik/src/input/simulation.xsd  
xsdcxx cxx-tree --std c++11 /home/kaiarenja/CLionProjects/MolekularDynamik/src/inputReader/simulation.xsd

## explaination for realocateParticles in LCParticleContainer
input: is an array of for integers that represent all 6 boundarys. each Integer can be 1(_outflow_ boundary) or 2(_refelcting_ boundary) or 3(_periodic_ boundary)


