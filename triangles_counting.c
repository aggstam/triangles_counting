// --------------------------------------------------------------------------------
//
// This program implements Scalable Triangle Counting algorithm
// introduced by S.Acerm, A.Yasar, S.Rajamanickam, M.Wolf and U.Catalyurek.
// In this version, algorithm calculation are performed in serial.
// Graph is read from an input file created by RandomGraph generator by S.Pettie 
// and V.Ramachandran.
//
// Author: Angelos Stamatiou, February 2020
//
// --------------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

FILE *fin; // Input file.
int nodes_count; // Graph nodes count.
int **matrix; // Graph nodes matrix.

// This function calculates triangles in a Graph, by performing (M*M).*M calculation, for
// the elements below the Graph matrix diagonal.
void calculate_triangles() {
	int i,j,k;
	
	int triangles_count = 0;
	int **mult = (int **)malloc(sizeof(int *) * nodes_count + sizeof(int) * nodes_count * nodes_count);
	int *ptr = (int *)(mult + nodes_count);
	for (i=0; i<nodes_count; i++) {
		mult[i] = (ptr + nodes_count * i);
		for (j=0; j<i; j++) {
			mult[i][j] = 0;	
			for (k=0; k<i; k++) {
				mult[i][j] += matrix[i][k] * matrix[k][j];
			}
			triangles_count += mult[i][j] * matrix[i][j];
		}
	}
	
	printf("Graph contains: %d triangles\n", triangles_count);
	free(mult);	
}

// This function initializes the Graph matrix, by reading the input file.
void initialize_matrix() {
	int i,j,fscanf_result;
	double w;

	matrix = (int **)malloc(sizeof(int *) * nodes_count + sizeof(int) * nodes_count * nodes_count);
	int *ptr = (int *)(matrix + nodes_count);
	for(i=0; i<nodes_count; i++) {
		matrix[i] = (ptr + nodes_count * i);		
		for (j=0; j<nodes_count; j++) {
			matrix[i][j] = 0;			
		}
	}			
	fscanf_result = fscanf(fin, "%d", &i);
	while (fscanf_result != 1 || i != -1) {
		fscanf_result = fscanf(fin, "%d %lf \n", &j, &w);
		if (i != j && matrix[i][j] == 0) {
			if (i < j) {
				matrix[j][i] = 1;
			} else if (i > j) {
				matrix[i][j] = 1;
			}
		}				
		fscanf_result = fscanf(fin, "%d", &i);
	}
}

// Auxiliary function that displays a message in case of wrong input parameters.
// Inputs:
//		char *compiled_name: Programms compiled name.
void syntax_message(char *compiled_name) {
	printf("Correct syntax:\n");
	printf("%s <input-file> \n", compiled_name);
	printf("where: \n");
	printf("<input-file> is the file containing a generated graph by RandomGraph that the algorithm will use.\n");
}

// This function checks run-time parameters validity and
// retrieves D-step value, input and output file names.
// Inputs:
//		char **argv: The run-time parameters.
// Output:
// 		1 --> Parameters read succussfully.
//		0 --> Something went wrong.
int read_parameters(char **argv) {
	char *input_filename = argv[1];
	if (input_filename == NULL) {
		printf("Input file parameter missing.\n");
		syntax_message(argv[0]);
		return 0;
	}
			
	fin = fopen(input_filename, "r");
	if (fin == NULL) {
		printf("Cannot open input file %s.\n", input_filename);
		return 0;		
	}

	printf("Counting Triangles of Graph retrieved from input file: %s\n", input_filename);
	return 1;
}

int main(int argc, char **argv) {	
	// Run-time parameters check.
	if (!read_parameters(argv)) {
		printf("Program terminates.\n");
		return -1;	
	}

	int fscanf_result = fscanf(fin, "%d \n", &nodes_count); // Retrieve Graph nodes count.
	if (fscanf_result != 1 || nodes_count > 0) {
		printf("Nodes count: %d\n", nodes_count);
		printf("Algorithm started, please wait...\n");
		initialize_matrix();
		clock_t t1 = clock();
		calculate_triangles(); // Perform algorith calculations.			
		clock_t t2 = clock();
		printf("Algorithm finished!\n");
		printf("Time spend: %f secs\n", ((float)t2 -t1) / CLOCKS_PER_SEC);	
	} else {
		printf("File is empty.\n");		
	}
	fclose(fin);
	printf("Program terminates.\n");
	return 0;
}
