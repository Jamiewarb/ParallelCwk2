/* This program relaxes all non-boundary values in an array
 * using MPI parallelism.
 * The root process acts as a master and sends a portion of the
 * array to each child process. Master and child processes then
 * all relax the portion of the array assigned
 * to them, and the child processes send their partial relaxations to 
 * the master, who constructs the final array.
 * If any number changes more than the value of precision, the process
 * is repeated, until it does not happen, at which point the program is finished.
 * During the repeated steps, only the start and end row for each section are sent
 * to the slaves, as those were changed by neighbouring slaves. The section
 * itself is retained by the slave, to minimise communication.
 **/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <mpi.h>

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_RESET   "\x1b[0m"

#define send_data_tag 2001
#define return_data_tag 2002

double fRand(double, double);
void printArray(double*, int, int);

struct ProcessPiece {
    int startRow, endRow, elements;
};

int main(int argc, char *argv[]) {
	
	MPI_Status status;
	MPI_Request request;
	int an_id, my_id, root_process, num_procs, start_row, end_row, 
		num_elements_to_send, num_elements_to_receive;

	/* Values hard coded - ensure to update or use command line parameters detailed in README
	 * debug - The level of debug output: 0, 1, 2, 3
	 * dimension - how big the square array is
	 * precision - how precise the relaxation needs to be before the program ends
	 *
	 * generateNumbers - 0 to use values in setValues[][], 1 to generate them randomly
	 * textFile[] - text file to read numbers from. Needs to be set and filled in if generateNumbers == 0
	 */

	int debug = 2; /* Debug output: 0 no detail - 1 some detail - 2 most detail - 3 step by step detail */
 
	int dimension = 1000;
	double precision = 0.1;

	int generateNumbers = 0;
	// textFile needs to be set and filled in if generateNumbers == 0
	char textFile[] = "../Values/values.txt";
	
	
	// Whether to use ANSI colour codes in output (turn off when using Balena's output logging)
	int useAnsi = 1;

	/* End editable values */

	


	// We need to use these in analysing the time the program takes
	double diff;

	/* Parse command line input */
	int a;
	for (a = 1; a < argc; a++) { /* argv[0] is program name */
		if (strcmp(argv[a], "-d") == 0 || strcmp(argv[a], "-dimension") == 0) {
			if (a + 1 <= argc - 1) { /* Make sure we have more arguments */
				if (atoi(argv[a+1]) > 0) {
					a++;
					dimension = atoi(argv[a]);
				} else {
					fprintf(stderr, "LOG WARNING - Invalid argument for -d. Positive integer required. Using %d dimension as default.\n", dimension);
				}
			}
		} else if (strcmp(argv[a], "-p") == 0 || strcmp(argv[a], "-precision") == 0) {
			if (a + 1 <= argc - 1) { /* Make sure we have more arguments */
				if (atof(argv[a+1]) > 0.0) {
					a++;
					precision = atof(argv[a]);
				} else {
					fprintf(stderr, "LOG WARNING - Invalid argument for -p. Positive double required. Using %f precision as default.\n", precision);
				}
			}
		} else if (strcmp(argv[a], "-g") == 0 || strcmp(argv[a], "-generate") == 0) {
			if (a + 1 <= argc - 1) { /* Make sure we have more arguments */
				if (atoi(argv[a+1]) >= 0) {
					a++;
					generateNumbers = atoi(argv[a]);
				} else {
					fprintf(stderr, "LOG WARNING - Invalid argument for -g. Integer >= 0 required. Using %d dimension as default.\n", dimension);
				}
			}
		} else if (strcmp(argv[a], "-f") == 0 || strcmp(argv[a], "-filepath") == 0) {
			if (a + 1 <= argc - 1) { /* Make sure we have more arguments */
				a++;
				strncpy(textFile, argv[a], sizeof(textFile));
			}
		} else if (strcmp(argv[a], "-debug") == 0) {
			if (a + 1 <= argc - 1) { /* Make sure we have more arguments */
				if (atoi(argv[a+1]) >= 0) {
					a++;
					debug = atoi(argv[a]);
				} else {
					fprintf(stderr, "LOG WARNING - Invalid argument for -debug. Integer >= 0 required. Using %d debug as default.\n", debug);
				}
			}
		} else {
			/* Non optional arguments here, but we have none of those */
		}
	}

	// Check our precision isn't too extreme, and our dimension is of a large enough size
	if (precision < 0.0000000001) precision = 0.0000000001;
	if (dimension < 3) dimension = 3;

	// From here on, all processes execute this code
	MPI_Init(NULL, NULL);

	double startMPI = MPI_Wtime();

	root_process = 0;

	// Who am I, and how many of us are there?
	MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    if (my_id == root_process) {

    	// Executed by the master process

		// We need to load the text file to read the numbers from
		FILE *valueFile;

		if (!generateNumbers) {
			valueFile = fopen(textFile, "r");
			if (valueFile == NULL) {
				fprintf(stdout, "LOG ERROR - Failed to open file: %s. Please check the path is valid and the file exists. Exiting program\n", textFile);
				return 1;
			}
		}

		/* Set up two arrays, one to store current results & one to store changes */
		double *values = malloc(dimension * dimension * sizeof(double));
		double *newValues = malloc(dimension * dimension * sizeof(double));
		if (!values || !newValues) {
			fprintf(stderr, "LOG ERROR - Failed to malloc. Exiting program.\n");
			return 1;
		}
		
		srand((unsigned)time(NULL));

		int i, j;
		// We need i*j numbers, in both arrays
		for (i = 0; i < dimension; i++) {
			for (j = 0; j < dimension; j++) {
				// If we're going to generate the numbers, or use predetermined ones
				if (generateNumbers)
					values[i*dimension+j] = fRand(1, 2);
				else {
					// Read in the next double. If it ends up being EOF, rewind and read again to same place.
					if (fscanf(valueFile, "%lf", &values[i*dimension+j]) <= 0) {
						rewind(valueFile);
						fscanf(valueFile, "%lf", &values[i*dimension+j]);
					}
				}
				// And copy them to the new array as well
				newValues[i*dimension+j] = values[i*dimension+j];
			}
		}
		
		if (debug >= 1) {
			fprintf(stdout, "LOG FINE - Using %d cores.\n", num_procs);
			fprintf(stdout, "LOG FINE - Using array of dimension %d.\n", dimension);
			fprintf(stdout, "LOG FINE - Working to precision of %.10lf.\n", precision);
		}

		if (debug >= 2) {
		/* Display initial array for debugging */
			fprintf(stdout, "LOG FINEST - Working with array:\n");
			printArray(values, dimension, useAnsi);
		}

		int count = 0; // Count how many times we try to relax the square array
		int withinPrecision = 0; // 1 when pass entirely completed within precision, i.e. finished

		// Set some initial values that only need setting once, and can be reused
		int avg_rows = (floor((dimension - 2) / num_procs)) + 2; // We add two as we need the previous and next rows for relaxing
		int extra_rows = (dimension - 2) % num_procs; // Get number of extra, remainder rows
		int initial_end_row = avg_rows - 1; // Set end row of master thread, decrement for index (row 1 -> index 0)
		int e_rows, my_rows;

		struct ProcessPiece pieces[num_procs];

		while (!withinPrecision) {
			count++;
			withinPrecision = 1;
			if (count == 1) {
				// First time distributing the array, so we need to send the whole portion

				e_rows = extra_rows; // Create a temporary variable that we can decrement
				end_row = initial_end_row; // Set the initial end_row to that of the master thread's end row

				for (an_id = 1; an_id < num_procs; an_id++) {
					// Evenly distribute rows, give first threads extra rows until none left
					my_rows = avg_rows;
					if (e_rows > 0) {
						my_rows++;
						e_rows--;
					}

					start_row = end_row - 1; // Previous end row is our first editable row + we need the un-editable one before for relaxing
					end_row   = start_row + (my_rows - 1); // We subtract 1 to get index from row number. End row is final, un-editable row.

					num_elements_to_send = my_rows*dimension; // Multiply by dimension as we need to send the columns per row too

					pieces[an_id].startRow = start_row;
					pieces[an_id].endRow = end_row;
					pieces[an_id].elements = num_elements_to_send;

					MPI_Isend( &num_elements_to_send, 1 , MPI_INT,
										an_id, send_data_tag, MPI_COMM_WORLD, &request );

					MPI_Isend( &values[start_row*dimension], num_elements_to_send, MPI_DOUBLE,
										an_id, send_data_tag, MPI_COMM_WORLD, &request );
				}
			} else {
				// We've already sent the bulk of the array before, so just send outer top + bottom rows that changed
				for (an_id = 1; an_id < num_procs; an_id++) {
					start_row = pieces[an_id].startRow;
					end_row = pieces[an_id].endRow;

					// MPI_Isend top row
					MPI_Isend( &values[start_row*dimension], dimension, MPI_DOUBLE,
											an_id, send_data_tag, MPI_COMM_WORLD, &request );
					// MPI_Isend bottom row
					MPI_Isend( &values[end_row*dimension], dimension, MPI_DOUBLE,
											an_id, send_data_tag, MPI_COMM_WORLD, &request );
				}
			}

			
			// Then we need to process our part of the array
			for (i = 1; i < avg_rows - 1; i++) { // Skip top and bottom rows
				for (j = 1; j < dimension - 1; j++) { // Skip left and right columns
					// Store relaxed number into new array
					newValues[i*dimension+j] = (values[(i-1)*dimension+j] + values[(i+1)*dimension+j] 
								  + values[i*dimension+(j-1)] + values[i*dimension+(j+1)]) / 4.0;
					/* If the numbers changed more than precision, we need to do it again */
					if (withinPrecision == 1 && fabs(values[i*dimension+j] - newValues[i*dimension+j]) > precision) {
						withinPrecision = 0;
					}
				}
			}

			// Now we need to receive the slave's portions
			e_rows = extra_rows;
			end_row = initial_end_row;
			int returnedWP = 1;
			int rows;

			// Make sure they've actually received it by now
			MPI_Wait(&request, &status);
		
			// Redo the same logic we did when handing out the rows to find out how many each process had
			for (an_id = 1; an_id < num_procs; an_id++) {

				// Get the rows and elements of this process that we saved earlier
				start_row = pieces[an_id].startRow;
				end_row = pieces[an_id].endRow;

				num_elements_to_receive = pieces[an_id].elements;

				// Store the new relaxed portioned array here 
				// Russel, can this be improved instead of malloc & free each time in the loop? num_elements_to_receive may change...
				double *tempVals = malloc(num_elements_to_receive * sizeof(double)); 
				if (!tempVals) {
					fprintf(stderr, "LOG ERROR - Failed to malloc. Exiting program.\n");
					return 1;
				}

				MPI_Recv( tempVals, num_elements_to_receive, MPI_DOUBLE, 
		               				an_id, return_data_tag, MPI_COMM_WORLD, &status);

				MPI_Recv( &returnedWP, 1, MPI_INT, 
									an_id, return_data_tag, MPI_COMM_WORLD, &status);

				// Check if we went out of precision
				if (withinPrecision && !returnedWP) {
					withinPrecision = 0;
				}

				rows = num_elements_to_receive / dimension;

				// Finally merge this portion of the array in to newValues
				for (i = 1; i < rows - 1; i++) { // Skip top and bottom rows
					for (j = 1; j < dimension - 1; j++) { // Skip left and right columns
						newValues[(start_row+i)*dimension+j] = tempVals[i*dimension+j];
					}
				}
				free(tempVals);
			}
			// We need to tell all slaves to finish or continue
			for (an_id = 1; an_id < num_procs; an_id++) {
				MPI_Isend( &withinPrecision, 1 , MPI_INT,
								an_id, send_data_tag, MPI_COMM_WORLD, &request);
			}


			MPI_Wait(&request, &status);

			// Finally swap the array pointers
			double *tempValues = values;
			values = newValues;
			newValues = tempValues;

			if (debug >= 3) {
				fprintf(stdout, "LOG GRANULAR - Array at step %d:\n", count);
				printArray(values, dimension, useAnsi);
			}
		}

		if (debug >= 1) 
			fprintf(stdout, "\nLOG FINE - Program complete. Relaxation count: %d.\n", count);
		if (debug >= 2) {
			fprintf(stdout, "LOG FINEST - Final array:\n");
			printArray(values, dimension, useAnsi);
		}
		
		// Your work here is done
		free(values);
		free(newValues);

		fprintf(stdout, "Program complete.\n");

		double endMPI = MPI_Wtime();
		diff = endMPI-startMPI;
		
		printf("elapsed time = %lf seconds\n\n", diff);

		MPI_Wait(&request, &status);


	} else {
		// Executed by all slaves
		int withinPrecision = 0;
		int firstTime = 1;
		int elements_to_receive = 0;

		MPI_Recv(&elements_to_receive, 1, MPI_INT, 
		               				root_process, send_data_tag, MPI_COMM_WORLD, &status);
		// Declare two arrays, one to receive current values & one to store changes
		double *values = malloc(elements_to_receive * sizeof(double));
		double *newValues = malloc(elements_to_receive * sizeof(double));

		int i, j;
		int rows = elements_to_receive / dimension;

		if (!values || !newValues) {
			fprintf(stderr, "LOG ERROR - Failed to malloc. Exiting program.\n");
			return 1;
		}

		while (!withinPrecision) {
			if (firstTime) {
				firstTime = 0;

		        MPI_Recv(values, elements_to_receive, MPI_DOUBLE, 
		               				root_process, send_data_tag, MPI_COMM_WORLD, &status);
		        for (i = 0; i < rows - 0; i++) {
					for (j = 0; j < dimension - 0; j++) {
						newValues[i*dimension+j] = values[i*dimension+j];
					}
				}
				
			} else {
				// We've already got the main bulk of the array so get extra rows
				MPI_Recv(values, dimension, MPI_DOUBLE, 
		               				root_process, send_data_tag, MPI_COMM_WORLD, &status);

				MPI_Recv(&values[(rows-1)*dimension], dimension, MPI_DOUBLE, 
		               				root_process, send_data_tag, MPI_COMM_WORLD, &status);

			}

			withinPrecision = 1; // Will be set to 0 if a number changes more than precision

			for (i = 1; i < rows - 1; i++) { // Skip top and bottom rows
				for (j = 1; j < dimension - 1; j++) { // Skip left and right columns
					// Store relaxed number into new array
					newValues[i*dimension+j] = (values[(i-1)*dimension+j] + values[(i+1)*dimension+j] 
								  + values[i*dimension+(j-1)] + values[i*dimension+(j+1)]) / 4.0;
					/* If the numbers changed more than precision, we need to do it again */
					if (withinPrecision == 1 && fabs(values[i*dimension+j] - newValues[i*dimension+j]) > precision) {
						withinPrecision = 0;
					}
				}
			}

			// Send back our relaxed array
			MPI_Send( newValues, elements_to_receive, MPI_DOUBLE, 
								root_process, return_data_tag, MPI_COMM_WORLD);
			
			// Send back our value for withinPrecision
			MPI_Send( &withinPrecision, 1, MPI_INT, 
								root_process, return_data_tag, MPI_COMM_WORLD);
		
			// Wait for master to tell us if we need to continue or stop
			MPI_Recv(&withinPrecision, 1, MPI_INT, 
		               				root_process, send_data_tag, MPI_COMM_WORLD, &status);

			double *tempValues = values;
			values = newValues;
			newValues = tempValues;

		}
		free(values);
		free(newValues);
	}
	// Finalize the MPI environment.
	MPI_Finalize();
	return 0;
}

void printArray(double *values, int dimension, int useAnsi) {
	int i, j;
	for (i = 0; i < dimension; i++) {
		for (j = 0; j < dimension; j++) {
			if (i == 0 || i == dimension-1 || j == 0 || j == dimension -1) {
				if (useAnsi) fprintf(stdout, ANSI_COLOR_RED);
			}
			fprintf(stdout, "%f ", values[i*dimension+j]);
			if (useAnsi) fprintf(stdout, ANSI_COLOR_RESET);
		}
		fprintf(stdout, "\n");
	}
}

double fRand(double fMin, double fMax) {
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}
