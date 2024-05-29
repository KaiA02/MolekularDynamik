MolSim
===

"GroupB"

Program Execution:

mkdir build && cd $_  
cmake ..  
make  
ctest  
./MolSim XML-input-file  
make doc

XML input-file consists of

1. Start_time as Double
2. End_Time as Double
3. Time_Step as Double
4. Output type: "xyz" or "vtk"
5. Input type: "cube" or "sonne"
6. (optional) performance measurement: true or false
7. (optional) log level: "off", "error", "warn", "info", "debug", "trace"

./MolSim ../input/eingabe.xml

create new c++ classes from xsd file:
xsdcxx cxx-tree --std c++11 /mnt/c/Users/joshu/CLionProjects/MolekularDynamik/src/input/simulation.xsd

./MolSim ../input/eingabe-sonne.txt 0 1000 0.014 1 1 0  
./MolSim ../input/eingabe-cube.txt 0 5 0.0002 1 2 0  
./MolSim ../input/eingabe-disk.txt 0 10 0.00005 1 3 0
