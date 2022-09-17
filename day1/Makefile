CC:=mpicc
CFLAGS:=-g -Wall -O3 
CLIBS:= -lmpiP -lm -lbfd -liberty -lunwind


blocking: laplace_mpi_blocking.c mesh.c solver.c
	${CC} ${CFLAGS} -o laplace_mpi_blocking laplace_mpi_blocking.c mesh.c solver.c -lm

nonblocking: laplace_mpi_nonblocking.c mesh.c solver.c
	${CC} ${CFLAGS} -o laplace_mpi_nonblocking laplace_mpi_nonblocking.c mesh.c solver.c -lm

persistent: laplace_mpi_persistent.c mesh.c solver.c
	${CC} ${CFLAGS} -o laplace_mpi_persistent laplace_mpi_persistent.c mesh.c solver.c -lm

blockingP: laplace_mpi_blocking.c mesh.c solver.c 
	${CC} ${CFLAGS} -o laplace_mpiP_blocking laplace_mpi_blocking.c mesh.c solver.c $(CLIBS)
    
nonblockingP: laplace_mpi_nonblocking.c mesh.c solver.c 
	${CC} ${CFLAGS} -o laplace_mpiP_nonblocking laplace_mpi_nonblocking.c mesh.c solver.c $(CLIBS)
    
persistentP: laplace_mpi_persistent.c mesh.c solver.c 
	${CC} ${CFLAGS} -o laplace_mpiP_persistent laplace_mpi_persistent.c mesh.c solver.c $(CLIBS)
    
clean:
	rm -f *.o laplace_mpi_blocking laplace_mpi_nonblocking laplace_mpi_persistent
