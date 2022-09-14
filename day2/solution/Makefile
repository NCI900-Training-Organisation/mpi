CC:=mpicc
CFLAGS:=-g -Wall -O3 -lm

win: laplace_mpi_win.c mesh.c solver.c
	${CC} ${CFLAGS} -o laplace_mpi_win laplace_mpi_win.c mesh.c solver.c

collective: laplace_mpi_win.c mesh.c solver.c
	${CC} ${CFLAGS} -o laplace_mpi_collective laplace_mpi_collective.c mesh.c solver.c

io: laplace_mpi_io.c mesh.c solver.c
	${CC} ${CFLAGS} -o laplace_mpi_io laplace_mpi_io.c mesh.c solver.c
	
clean:
	rm -f *.o laplace_mpi_win