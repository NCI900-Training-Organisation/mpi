CC:=mpicc
CFLAGS:=-g -Wall -O3 
CLIBS:= -lmpiP -lm -lbfd -liberty -lunwind

win: laplace_mpi_win.c mesh.c solver.c
	${CC} ${CFLAGS} -o laplace_mpi_win laplace_mpi_win.c mesh.c solver.c -lm

collective: laplace_mpi_win.c mesh.c solver.c
	${CC} ${CFLAGS} -o laplace_mpi_collective laplace_mpi_collective.c mesh.c solver.c -lm

io: laplace_mpi_io.c mesh.c solver.c
	${CC} ${CFLAGS} -o laplace_mpi_io laplace_mpi_io.c mesh.c solver.c -lm

winP: laplace_mpi_win.c mesh.c solver.c 
	${CC} ${CFLAGS} -o laplace_mpiP_win laplace_mpi_win.c mesh.c solver.c $(CLIBS)
clean:
	rm -f *.o laplace_mpi_win