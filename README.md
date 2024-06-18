# MolSim
=========

"GroupB"

## Program Execution:

1. mkdir build && cd $_
2. cmake ..
3. make
4. ctest
5. ./MolSim XML-input-file (e.g ./MolSim ../input/eingabe-sonne.xml)
6. make doc

## XML input-file consists of:

1. Start_time as Double
2. End_Time as Double
3. Time_Step as Double
4. Output type: "xyz" or "vtk"
5. Input type: "cube" or "sonne"
6. (optional) performance measurement: true or false
7. (optional) log level: "off", "error", "warn", "info", "debug", "trace"


## Quick start guide:
1. adapt inputFile to your need
2. follow Steps described in Program Execution
3. open Paraview
4. open generated outputFiles in Paraview and add a Glyph to the Data
5. run Simulation

## Thermostat:
- First set the "thermostatON" option in xml on "YES"
- In the beginning the system is set to the initial temperature.
- n_thermostat means after how many steps the thermostat is applied.
- Reaching and holding the target temperature is the aim, therefore we have the delta_temp to say what the absolute change of temperature allowed is (for the steps, declared at n_thermostat)
- If no delta_temp is given, then we see declare the delta_temp as 0 in xml file, but see it as infinity
- If n_thermostat = 0, then the thermostat is applied directly.


## create new c++ classes from xsd file:
xsdcxx cxx-tree --std c++11 /mnt/c/Users/joshu/CLionProjects/MolekularDynamik/src/input/simulation.xsd d
xsdcxx cxx-tree --std c++11 /home/kaiarenja/CLionProjects/MolekularDynamik/src/inputReader/simulation.xsd

## explaination for realocateParticles in LCParticleContainer
input: is an array of for integers that represent all 6 boundarys. each Integer can be 1(outflow boundary) or 2(refelcting boundary) or 3(periodic boundary)

## reflecting Boundary
still working on the implementation.
in LCParticleContainer have a look at getBoundaryParticles(), handleBoundaryAction(), getInfluencingBoundarys(Particle* p)
Idea is to call handleBoundaryAction() after calculating LJF Forces on all particles. in the Methode handleBoundaryAction() all Boundary Particles are collected by getBoundaryParticles(). For each Particle the Borders, that will affect the Particle are collected. In the Future we have to generate a Halo Particle (look at handleBoundaryAction() to see specification for haloParticle) and calculate LJF between these two Particles and add force to Boundary Particle. But only generate Halo Particle, when distance to wall ist smaller than ((root 6 of 2) * omega) / 2.


