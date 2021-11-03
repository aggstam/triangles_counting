This program implements Scalable Triangle Counting algorithm introduced by S.Acerm, A.Yasar, S.Rajamanickam, M.Wolf and U.Catalyurek [1].
<br>
Graph is read from an input file created by RandomGraph generator by S.Pettie and V.Ramachandran [2].
<br>
Two implementations are included, one executing the algorithm in serial, and one using the MPI Standard.
<br>
Graph Matrix is conformally partioned to blocks, based on MPI processes.
<br>
Each process is assigned with a block and calculate a local triangles count.
<br>
In order to calculate the local triangles count, required blocks are retrieved from other processes, based on the algorithm.
<br>
<br>
References:
<br>
[1] S. Acer, A. Yaşar, S. Rajamanickam, M. Wolf and Ü. V. Catalyürek, "Scalable Triangle Counting on Distributed-Memory Systems," 2019 IEEE High Performance Extreme Computing Conference (HPEC), 2019, pp. 1-5, doi: 10.1109/HPEC.2019.8916302.
<br>
[2] http://www.dis.uniroma1.it/challenge9/code/Randgraph.tar.gz
