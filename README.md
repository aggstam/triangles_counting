# triangles_counting
This program implements Scalable Triangle Counting algorithm introduced by S.Acerm, A.Yasar, S.Rajamanickam, M.Wolf and U.Catalyurek [1].
<br>
Graph is read from an input file created by RandomGraph generator by S.Pettie and V.Ramachandran [2].
<br>
Two implementations are included, one executing the algorithm in serial, and one using the MPI Standard.
<br>
MPI implementation requires *openmpi* package to be installed.
<br>
Graph Matrix is conformally partioned to blocks, based on MPI processes.
<br>
Each process is assigned with a block and calculate a local triangles count.
<br>
In order to calculate the local triangles count, required blocks are retrieved from other processes, based on the algorithm.

## Usage
Both version can be invocted via the Makefile, or by directly compiling and executing.

### Make usage
#### Normal code
```shell
$ make
```
To include a different input file:
```shell
$ make FILE={file_path}
```

#### MPI code
```shell
$ make mpi
```
To configure different how many processes to use:
```shell
$ make mpi PROCESSES={processes}
```
To include a different input file:
```shell
$ make mpi FILE={file_path}
```

### Direct usage
#### Normal code
Compilation:
```shell
$ gcc -o triangles_counting triangles_counting.c
```
Execution:
```shell
$ ./triangles_counting {input_file}
```

#### MPI code
Compilation:
```shell
$ mpicc -lm -o mpi_triangles_counting mpi_triangles_counting.c
```
Execution:
```shell
$ mpiexec -np {processes} ./mpi_triangles_counting {input_file}
```

## Execution examples
### Normal code
```shell
$ make
Executing normal code...
gcc -o triangles_counting triangles_counting.c
./triangles_counting grph_triangles
Counting Triangles of Graph retrieved from input file: grph_triangles
Nodes count: 1000
Algorithm started, please wait...
Graph contains: 14 triangles
Algorithm finished!
Time spend: 1.403673 secs
Program terminates.
```

### MPI code
```shell
$ make mpi
Executing MPI code...
mpicc -lm -o mpi_triangles_counting mpi_triangles_counting.c
mpiexec -np 4 ./mpi_triangles_counting grph_triangles
Counting Triangles of Graph retrieved from input file: grph_triangles
Nodes count: 1000
Algorithm started, please wait...
Graph contains: 14 triangles
Algorithm finished!
Time spend: 1.029489 secs
```

## References
[1] S. Acer, A. Yaşar, S. Rajamanickam, M. Wolf and Ü. V. Catalyürek, "Scalable Triangle Counting on Distributed-Memory Systems," 2019 IEEE High Performance Extreme Computing Conference (HPEC), 2019, pp. 1-5, doi: 10.1109/HPEC.2019.8916302.
<br>
[2] http://www.dis.uniroma1.it/challenge9/code/Randgraph.tar.gz
