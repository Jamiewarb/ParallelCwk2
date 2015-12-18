# Parallel Computing - Coursework 2 (Open MPI)
The second coursework for Final Year Bath Comp Sci, Parallel Computing unit

This zip contains 3 folders:  
Source, Values, Testing

Source - the actual compilable .c code  
Values - the .txt file with the numbers that the program will use  
Testing - the testing document and the values utilised by it (code run times etc)  


Compile sequential.c using gcc -Wall -pthread -o sequential sequential.c -lrt  
Compile parallel.c using mpicc -Wall -o ompi-parallel ompi-parallel.c -lrt

Both programs can be run with these possible flags:  
-debug : The level of debug output: 0, 1, 2 (default: 0)  
-d : length of the square array (default: 10)  
-p : how precise the relaxation needs to be before the program ends (default: 0.0001)  
-g : (1 or 0) 0 to use values in from file specified in program, 1 to generate them randomly (default: 0)  
-f : string, path of the textfile to use (default: ../Values/values.txt)  

For example: 
* ./ompi-parallel -debug 2 -d 500 -p 0.01 -g 0 -f values.txt
* 
Excluding a flag uses the default value

By default, the program reads from a text file, as generating numbers makes testing more obtuse.  
Capability to run with generated numbers is, however, provided through -g 1.