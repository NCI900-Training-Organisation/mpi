/* =================================================================
solver.c

Solve a model 2D Poisson equaton with Dirichlet boundary condition.

-Delta u = 2pi^2 * sin(pi x)sin(pi y) in [0,1]^2
       u = sin(pi x) sin(y) on boundary

The problem is discretised over a uniform mesh by finite difference 
method and the resulting linear system is solved by choices of Jacobi
or Gauss-Seidel.


Compile:  mpicc -g -Wall -O3 -lm -o fd_laplace-mpi_block fd_laplace-mpi_block.c 

Usage:  mpirun -np 4 ./laplace_mpi_blocking size tolerance method

Produced for NCI Training. 

Frederick Fung 2022
4527FD1D
====================================================================*/


double local_L2_residual(const int *ptr_to_rows, 
                         int mesh_size, 
                         double space, 
                         const double *ptr_submesh, 
                         const double *ptr_subrhs);

void Jacobi(int *ptr_to_rows, 
            int mesh_size, 
            double *ptr_submesh,  
            double *ptr_submesh_new, 
            const double *ptr_rhs, 
            double space);

void Jacobi_int(int *ptr_to_rows, 
                int mesh_size, 
                double *ptr_submesh,  
                double *ptr_submesh_new, 
                const double *ptr_rhs, 
                double space);

void Jacobi_top(int *ptr_to_rows, 
                int mesh_size, 
                double *ptr_submesh,  
                double *ptr_submesh_new, 
                const double *ptr_rhs, 
                double space);


void Jacobi_bottom(int *ptr_to_rows, 
                   int mesh_size, 
                   double *ptr_submesh,  
                   double *ptr_submesh_new, 
                   const double *ptr_rhs, 
                   double space);