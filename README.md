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
-debug : The level of debug output: 0, 1, 2, 3 (*default: 1*)  
-d : integer length of the square array (*default: 1000*)  
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

####How does ompi-parallel code work?
* Initially, the master thread creates a matrix from an input file, looping through the file more than once if necessary.  
* Next, the master sets a start and end row for itself and each other process, then sends out the specific chunks of
the matrix to each process, including an extra row top and bottom, used in the relaxation process.  
* The master and all the slaves then relax their own portion of the matrix.  
* Each process tells the master if it went out of precision or not. If any process went out of precision, the Master tells everyone
that they need to continue.  
* When they continue, each process sends its top and bottom editable rows to its neighbouring processes. It also receives a top and bottom
un-editable row from those neighbouring processes. This is all the thread needs for its next iteration.  
* Once an iteration occurs where all processes stay within precision, the processes finally send their chunk of the matrix back to
the master thread, who combines it all for the finished result.  
