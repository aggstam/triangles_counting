FILE = grph_triangles
PROCESSES = 4

all:
	$(info Executing normal code...)
	gcc -o triangles_counting triangles_counting.c
	./triangles_counting $(FILE)

mpi:
	$(info Executing MPI code...)
	mpicc -lm -o mpi_triangles_counting mpi_triangles_counting.c
	mpiexec -np $(PROCESSES) ./mpi_triangles_counting $(FILE)

clean:
	rm -f triangles_counting mpi_triangles_counting

.PHONY: all mpi clean
