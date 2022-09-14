/* =================================================================
mesh.c

Solve a model 2D Poisson equaton with Dirichlet boundary condition.

-Delta u = 2pi^2 * sin(pi x)sin(pi y) in [0,1]^2
       u = sin(pi x) sin(y) on boundary

The problem is discretised over a uniform mesh by finite difference 
method and the resulting linear system is solved by choices of Jacobi
or Gauss-Seidel.


Compile:  mpicc -g -Wall -O3 -lm -o fd_laplace-mpi_block fd_laplace-mpi_block.c 

Usage:  mpirun -np 4 ./fdd_laplace-mpi size tolerance method

Produced for NCI Training. 

Frederick Fung 2022
4527FD1D
====================================================================*/


double bnd_fc(int x, int y, double space);

double rhs_fc(int x, int y, double space);

void init_mesh(const int mesh_size,
               double submesh[][mesh_size], 
               double submesh_new[][mesh_size],
               double subrhs[][mesh_size],
               const int rank,
               const int cells,
               const int int_rows,
               const double space,
               int *ptr_rows);


