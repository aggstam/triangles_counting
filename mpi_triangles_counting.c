// --------------------------------------------------------------------------------
//
// This program implements Scalable Triangle Counting algorithm
// introduced by S.Acerm, A.Yasar, S.Rajamanickam, M.Wolf and U.Catalyurek.
// In this version, MPI standar is used for parallelizing the algorithm.
// Graph Matrix is conformally partioned to blocks, based on MPI processes.
// Each process is assigned with a block and calculate a local triangles count.
// In order to calculate the local triangles count, required blocks are retrieved
// from other processes, based on the algorithm.
// Graph is read from an input file created by RandomGraph generator by S.Pettie 
// and V.Ramachandran.
//
// Author: Angelos Stamatiou, February 2020
//
// --------------------------------------------------------------------------------

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

FILE *fin;                          // Input file.
int nodes_count;                    // Graph nodes count.
int blocks_count;                   // Blocks to process count.
int **matrix;                       // Graph nodes matrix.
int **block;                        // Assigned block.
int **multiplication_blocks;        // Received blocks for matrix multiplication.
int **el_multiplication_blocks;     // Received blocks for multiplication per element.
int **qmap;                         // Mapping of assigned blocks to processes.
int q_i,q_j;                        // Position of assigned processes.
int rank;                           // MPI Process rank.
int size;                           // MPI Process count.
int q;                              // Blocks table size.
int n;                              // Block size.
int mult_blocks_receive_counter;    // Blocks received count for matrix multiplication.
int el_mult_blocks_receive_counter; // Blocks received count for multiplication per element.
int triangles_count;                // Triangles found by each process.

// This function calculates triangles in a Graph, by performing (M*M).*M calculation, for
// the elements below the Graph matrix diagonal.
void calculate_triangles()
{
    int triangles_count = 0;
    int **mult = (int**)malloc(sizeof(int*) * nodes_count + sizeof(int) * nodes_count * nodes_count);
    if (mult == NULL) {
      printf("Error: malloc for mult failed.\n");
      exit(1);
    }
    int *ptr = (int*)(mult + nodes_count);
    for (int i = 0; i < nodes_count; i++) {
        mult[i] = (ptr + nodes_count * i);
        for (int j = 0; j < i; j++) {
            mult[i][j] = 0;    
            for (int k = 0; k < i; k++) {
                mult[i][j] += matrix[i][k] * matrix[k][j];    
            }
            triangles_count += mult[i][j] * matrix[i][j];
        }
    }
    
    printf("Graph contains: %d triangles\n", triangles_count);
    free(mult);    
}

// This function initializes the Graph matrix, by reading the input file.
void initialize_matrix()
{
    int i,j,fscanf_result;
    double w;
    matrix = (int**)malloc(sizeof(int*) * nodes_count + sizeof(int) * nodes_count * nodes_count);
    if (matrix == NULL) {
      printf("Error: malloc for matrix failed.\n");
      exit(1);
    }
    int *ptr = (int*)(matrix + nodes_count);
    for(i = 0; i < nodes_count; i++) {
        matrix[i] = (ptr + nodes_count * i);        
        for (j = 0; j < nodes_count; j++) {
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
//      char *compiled_name: Programms compiled name.
void syntax_message(char *compiled_name)
{
    printf("Correct syntax:\n");
    printf("%s <input-file> \n", compiled_name);
    printf("where: \n");
    printf("<input-file> is the file containing a generated graph by RandomGraph that the algorithm will use.\n");
}

// This function checks run-time parameters validity and
// retrieves D-step value, input and output file names.
// Inputs:
//      char **argv: The run-time parameters.
// Output:
//      1 --> Parameters read succussfully.
//      0 --> Something went wrong.
int read_parameters(char **argv)
{
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

// This function reduces local triangles counts to master process P0.
// If all process have been assigned with blocks, MPI_Reduce is used,
// else each assigned slave process sends its local triangles count to the master process.
void reduce_triangles_counts()
{
    int total_triangles=0;
    if (blocks_count==size) {
        MPI_Reduce(&triangles_count, &total_triangles, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);        
    } else {
        if (rank==0) {
            MPI_Status status;
            int process_triangles,i;
            for (i = 1; i < blocks_count; i++) {
                MPI_Recv(&process_triangles, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
                total_triangles += process_triangles;
            }
            total_triangles += triangles_count;
        } else {
            MPI_Send(&triangles_count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);        
        }        
    }
    if (rank == 0) {
        printf("Graph contains: %d triangles\n", total_triangles);    
    }    
}

// In this function, master process calculates its local triangles count, by performing 
// the (M*M).*M calculation, based on the algorithm.
void calculate_local_triangles_master()
{
    int i,j,k;
    for (i = 1; i < q; i++) {
        MPI_Send(&block[0][0], n*n, MPI_INT, qmap[i][0], 0, MPI_COMM_WORLD);
    }

    triangles_count = 0;
    int** mult = (int**)malloc(sizeof(int*) * nodes_count + sizeof(int) * n * n);
    if (mult == NULL) {
      printf("Error: malloc for mult failed.\n");
      exit(1);
    }
    int* ptrm = (int*)(mult + nodes_count);
    for (i=0; i<n; i++) {
        mult[i] = (ptrm + n * i);
        for (j = 0; j < n; j++) {
            mult[i][j] = 0;    
            for (k = 0; k < n; k++) {
                mult[i][j] += block[i][k] * block[k][j];    
            }
            triangles_count += mult[i][j] * block[i][j];
        }
    }
    
    free(mult);    
    free(matrix);
    free(block);    
    free(qmap);
}

// In this function, each slave process calculates is local triangles count, by performing
// the (M*M).*M calculation, based on the algorithm.
void calculate_local_triangles_slave()
{
    int i,j,k,b;
    triangles_count = 0;
    int **mult = (int**)malloc(sizeof(int *) * n + sizeof(int) * n * n);
    if (mult == NULL) {
      printf("Error: malloc for mult failed.\n");
      exit(1);
    }
    int *mult_ptr = (int*)(mult + n);

    for (i = 0; i < n; i++) {
        mult[i] = (mult_ptr + n * i);
    }

    for (b = 0; b < mult_blocks_receive_counter; b++) {
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                mult[i][j]=0;
                for (k = 0; k < n; k++) {
                    mult[i][j] += block[i][k] * multiplication_blocks[k][(b*n)+j];        
                }
                triangles_count += mult[i][j] * el_multiplication_blocks[i][(b*n)+j];
            }
        }
    }
    
    // Computing remaining Block(assigned) for processes in the Q map diagonal, as multiplication_blocks 
    // and el_multiplication_blocks matrices are not the same size.
    if (q_i == q_j) {
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                mult[i][j] = 0;    
                for (k = 0; k < n; k++) {
                    mult[i][j] += block[i][k] * block[k][j];    
                }
                triangles_count += mult[i][j] * block[i][j];
            }
        }    
    }
    
    free(mult);
    free(block);
    free(multiplication_blocks);
    free(el_multiplication_blocks);
    free(qmap);    
}

// In this function, each slave process retrieves required Blocks from other processes,
// in order to perform its calculations. Each process finds which processes have its required
// Blocks, based on the Q Map. Each process requires blocks from the row of its column
// in the Q Map, for the matrix multiplication, and the previous blocks of its row in the Q Map,
// for the multiplication per element.
void receive_blocks_from_other_processes()
{
    int i,j,k;
    MPI_Status status;

    // Allocate memory for retrieving required Blocks.
    multiplication_blocks = (int**)malloc(sizeof(int*) * n + sizeof(int) * n * ((q_j+1)*n));
    if (multiplication_blocks == NULL) {
      printf("Error: malloc for multiplication_blocks failed.\n");
      exit(1);
    }
    int *multiplication_blocks_ptr = (int*)(multiplication_blocks + n);
    el_multiplication_blocks = (int**)malloc(sizeof(int*) * n + sizeof(int) * n * ((q_j+1)*n));
    if (el_multiplication_blocks == NULL) {
      printf("Error: malloc for el_multiplication_blocks failed.\n");
      exit(1);
    }
    int *el_multiplication_blocks_ptr = (int*)(el_multiplication_blocks + n);
    int **buffer = (int**)malloc(sizeof(int*) * n + sizeof(int) * n * n);
    if (buffer == NULL) {
      printf("Error: malloc for buffer failed.\n");
      exit(1);
    }
    int *buffer_ptr = (int*)(buffer + n);
    for (i = 0; i < n; i++) {
        multiplication_blocks[i] = (multiplication_blocks_ptr + ((q_j+1)*n) * i);
        el_multiplication_blocks[i] = (el_multiplication_blocks_ptr + ((q_j+1)*n) * i);
        buffer[i] = (buffer_ptr + n * i);    
    }

    // Retrieve blocks from the row of process column in Q Map.
    // Processes in column 0 require only the Block of master process.
    mult_blocks_receive_counter=0;
    if (q_j == 0) {
        mult_blocks_receive_counter++;        
        MPI_Recv(*buffer, n*n, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                multiplication_blocks[i][j] = buffer[i][j];
            }
        }
    } else {
        for (j = 0; j < q_j + 1; j++) {
            if (qmap[j][q_i] != rank) {
                mult_blocks_receive_counter++;    
                MPI_Recv(*buffer, n*n, MPI_INT, qmap[q_j][j], 0, MPI_COMM_WORLD, &status);
                for (i = 0; i < n; i++) {
                    for (k = 0; k < n; k++) {
                        multiplication_blocks[i][(j*n)+k] = buffer[i][k];
                    }
                }
            }        
            
        }
    }
    
    // Send assigned Blocks to the processes in the column of process row in Q Map. 
    for (i = q_i; i < q; i++) {
        if (qmap[i][q_i] != rank) {
            MPI_Send(&block[0][0], n*n, MPI_INT, qmap[i][q_i], 0, MPI_COMM_WORLD);        
        }            
    }

    // Retrieve blocks from the previous processes in the row of process in Q Map.
    el_mult_blocks_receive_counter=0;
    for (j = 0; j < q_j; j++) {
        el_mult_blocks_receive_counter++;
        MPI_Recv(*buffer, n*n, MPI_INT, qmap[q_i][j], 1, MPI_COMM_WORLD, &status);
        for (i = 0; i < n; i++) {
            for (k = 0; k < n; k++) {
                el_multiplication_blocks[i][(j*n)+k] = buffer[i][k];
            }
        }
    }

    // Copy assigned Block in last position of el_multiplication_blocks.
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            el_multiplication_blocks[i][(el_mult_blocks_receive_counter*n)+j] = block[i][j];
        }
    }

    // Send assigned Block to the next processes in the row of process in Q Map.
    for (j = q_j + 1; j < q_i + 1; j++) {
        MPI_Send(&block[0][0], n*n, MPI_INT, qmap[q_i][j], 1, MPI_COMM_WORLD);
    }
    free(buffer);
}

// This function retrieves Q mapping from master process.
void retrieve_q_mapping()
{
    int i,j;
    qmap = (int**)malloc(sizeof(int*) * q + sizeof(int) * q * q);
    if (qmap == NULL) {
      printf("Error: malloc for qmap failed.\n");
      exit(1);
    }
    int *qptr = (int*)(qmap + q);
    for (i = 0; i < q; i++) {
        qmap[i] = (qptr + q * i);
    }
    // Retrieve Q mapping from P0.
    MPI_Bcast(*qmap, q*q, MPI_INT, 0, MPI_COMM_WORLD);
    for (i = 0; i < q; i++) {
        for (j = 0; j < i+1; j++) {
            if (qmap[i][j] == rank) {
                q_i=i;
                q_j=j;
                i=q;
                j=q+1;
            }            
        }            
    }
}

// This function retrieves Process assigned Block.
void receive_block()
{
    int i,j,k;
    MPI_Status status;
    // Retrieve Block size in order to allocate memory for retrieving the Block.
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    block = (int**)malloc(sizeof(int*) * n + sizeof(int) * n * n);
    if (block == NULL) {
      printf("Error: malloc for block failed.\n");
      exit(1);
    }
    int* ptr = (int*)(block + n);
    for (i = 0; i < n; i++) {
        block[i] = (ptr + n * i);            
    }
    // Retrieve assigned Block by P0.
    MPI_Recv(*block, n*n, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);    
}

// This function, used by master process, calculates Q mapping and send each process its assigned Block.
void scatter_blocks()
{
    int i,j,p,q_row,q_column,p_row,p_column,counter;
    
    qmap = (int**)malloc(sizeof(int*) * q + sizeof(int) * q * q);
    if (qmap == NULL) {
      printf("Error: malloc for qmap failed.\n");
      exit(1);
    }
    int *qptr = (int*)(qmap + q);
    for (i = 0; i < q; i++) {
        qmap[i] = (qptr + q * i);
    }

    // Broadcast Block size to rest processes.
    n = nodes_count/q;
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);    
    block = (int**)malloc(sizeof(int*) * n + sizeof(int) * n * n);
    if (block == NULL) {
      printf("Error: malloc for block failed.\n");
      exit(1);
    }
    int *ptr = (int*)(block + n);
    for (i=0; i<n; i++) {
        block[i] = (ptr + n * i);
    }

    // Calculate Q mapping, and send each process its assigned Block.
    // Master process will position each process Id sequentially in a top to bottom manner.
    // E.G a 4x4 Q Map will look like:
    //             [0, 0, 0, 0]
    //            [1, 4, 0, 0]
    //            [2, 5, 7, 0]
    //            [3, 6, 8, 9]
    q_row = 1;
    q_column = 0;
    counter = 0;
    for (p = 1; p < size; p++) {
        if (counter < blocks_count-1) {
            p_row=0;
            p_column=0;
            for (i = (q_row*n); i < (q_row*n) + n; i++) {
                for (j = (q_column*n); j < (q_column*n) + n; j++) {
                    block[p_row][p_column] = matrix[i][j];
                    p_column++;
                }
                p_column=0;
                p_row++;        
            }
            MPI_Send(&block[0][0], n*n, MPI_INT, p, 0, MPI_COMM_WORLD);
            qmap[q_row][q_column] = p;
            if (q_row<q-1) {                    
                q_row++;
            } else {
                q_column++;
                q_row = q_column;                    
            }
            counter++;
        }                    
    }

    // Master process keeps Matrix first block.
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            block[i][j] = matrix[i][j];
        }        
    }

    // Broadcast Q Mapping to rest processes.
    MPI_Bcast(&qmap[0][0], q*q, MPI_INT, 0, MPI_COMM_WORLD);
}

int main(int argc, char **argv)
{
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    
    // Each process calculates Q value and blocks count, 
    // in order to know if it will be assigned with a block.    
    q = (-1+sqrt(1+8*size))/2; // P=Q(Q+1)/2=>2size=q^2+q=>q^2+q-2size=0
    blocks_count=((q*(q-1))/2)+q; // Blocks of Q Table that will be computed.
    if (rank == 0) {
        // Run-time parameters check.
        if (!read_parameters(argv)) {
            printf("Program terminates.\n");
            MPI_Abort(MPI_COMM_WORLD, -1);    
        }
        int fscanf_result = fscanf(fin, "%d \n", &nodes_count); // Retrieve Graph nodes count.
        if (fscanf_result != 1 || nodes_count > 0) {
            printf("Nodes count: %d\n", nodes_count);
            // Graph matrix is required to be conformally partitioned to blocks.
            if (nodes_count%q==0) {
                printf("Algorithm started, please wait...\n");
                initialize_matrix();
                clock_t t1 = clock();
                // When q==1, only one process is required.    
                if (q > 1) {
                    scatter_blocks(); // Calculate Q mapping and send each process its assigned Block.
                    calculate_local_triangles_master(); // Calculate local triangles count.
                    reduce_triangles_counts(); // Retrieve local triangles count from other processes.
                } else {
                    calculate_triangles(); // Perform serial algorithm.
                }                                                
                clock_t t2 = clock();
                printf("Algorithm finished!\n");
                printf("Time spend: %f secs\n", ((float)t2 -t1) / CLOCKS_PER_SEC);
                fclose(fin);        
            } else {
                printf("Graph cannot be conformally partitioned.\n");
                printf("Choose a different processes size.\n");
                fclose(fin);
                MPI_Abort(MPI_COMM_WORLD, -1);
            }
        } else {
            printf("File is empty.\n");
            fclose(fin);
            MPI_Abort(MPI_COMM_WORLD, -1);
        }    
    } else if (q > 1 && rank < blocks_count) {
        receive_block(); // Retrieve assigned block.
        retrieve_q_mapping(); // Retrieve Q mapping.
        receive_blocks_from_other_processes(); // Retrieve required blocks from other processes.
        calculate_local_triangles_slave(mult_blocks_receive_counter); // Calculate local triangles count.
        reduce_triangles_counts(); // Send local triangles count to process P0.
    }
    MPI_Finalize();
}
