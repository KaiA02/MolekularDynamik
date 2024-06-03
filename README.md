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


## create new c++ classes from xsd file:
xsdcxx cxx-tree --std c++11 /mnt/c/Users/joshu/CLionProjects/MolekularDynamik/src/input/simulation.xsd






