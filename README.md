MolSim
===

"GroupB"

Program Execution: 

mkdir build && cd $_  
cmake ..  
make  
ctest  
./MolSim input  
make doc  
  
Input consist of  
1. Path to Inputfile  
2. Start_time as Double
3. End_Time as Double  
4. Time_Step as Double  
5. "1" for Output in vtk, "2" for Output in xyz
6. "1" for Input like Assignment 1, "2" for Input like Cuboids


./MolSim ../input/eingabe-cube.txt 0 100 0.0002 1 2
