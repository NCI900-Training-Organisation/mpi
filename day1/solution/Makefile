CC:=mpicc
CFLAGS:=-g -Wall -O3 -lm

blocking: laplace_mpi_blocking.c mesh.c solver.c
	${CC} ${CFLAGS} -o laplace_mpi_blocking laplace_mpi_blocking.c mesh.c solver.c

nonblocking: laplace_mpi_nonblocking.c mesh.c solver.c
	${CC} ${CFLAGS} -o laplace_mpi_nonblocking laplace_mpi_nonblocking.c mesh.c solver.c

persistent: laplace_mpi_persistent.c mesh.c solver.c
	${CC} ${CFLAGS} -o laplace_mpi_persistent laplace_mpi_persistent.c mesh.c solver.c

clean:
	rm -f *.o laplace_mpi_blocking laplace_mpi_nonblocking laplace_mpi_persistent
