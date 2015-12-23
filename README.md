# Parallel Computing - Coursework 2 (Open MPI)
The second coursework for Final Year Bath Comp Sci, Parallel Computing unit

####Contents
This zip contains 3 folders:  
Source, Values and Testing

**Source** - the actual compilable .c code  
**Values** - the .txt file with the numbers that the program will use  
**Testing** - the testing document, the values utilised by it and the test output logs  

####Compiling and Running
Compile sequential.c using `gcc -Wall -o sequential sequential.c -lrt`  
Compile ompi-parallel.c using `mpicc -Wall -o ompi-parallel ompi-parallel.c`

Both programs can be run with these possible flags:  
-debug : The level of debug output: 0, 1, 2, 3 (*default: 0*)  
-d : integer length of the square array (*default: 10*)  
-p : how precise the relaxation needs to be before the program ends, as a double (*default: 0.0001*)  
-g : (1 or 0) 0 to use values from file specified, 1 to generate them randomly (*default: 0*)  
-f : string, path of the textfile to use (*default: ../Values/values.txt*)  

For example: 
* `./sequential -debug 2 -d 500 -p 0.01 -f values.txt`
* `mpirun -np 16 ./ompi-parallel -debug 1 -d 10000 -p 0.001`  

Excluding a flag will use the default value.

####Value Sets
By default, the program reads from the text file at Values/values.txt, as generating numbers makes testing more obtuse.  
Capability to run with generated numbers is, however, provided through -g 1.  

The program will read doubles from a file until enough have been read.  
Reaching EOF will rewind the file and continue reading from the start.