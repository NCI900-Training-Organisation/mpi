/* =================================================================
solver.c

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


#include<stdio.h>
#include <stdlib.h>
#include<string.h>
#include<math.h>
#include<mpi.h>
#include "solver.h"

#define MPIIO
//#define MPI_DEBUG


double local_L2_residual(const int *ptr_to_rows, int mesh_size, double space, const double *ptr_submesh, const double *ptr_subrhs)
{
    
    /* discrete L2 norm of the residual by a given approximation */
    double residual = 0.0f;

      for (int i = 1; i< *ptr_to_rows -1; i++){
        for (int j = 1; j< mesh_size-1; j++){
             
             double diff =  ( 4* *(ptr_submesh +  (i * mesh_size) +j ) - *(ptr_submesh + ((i-1)*mesh_size) +j) - *(ptr_submesh +((i+1)*mesh_size)+j)\
                             - *(ptr_submesh + (i *mesh_size) + j-1 ) - *(ptr_submesh + (i*mesh_size) +j +1)) / (space *space )- *(ptr_subrhs + i*mesh_size +j);
             diff = space *space * diff;
             residual += pow(diff, 2);
              
        }
    }
     
    return residual;
}


void Jacobi(int *ptr_to_rows, int mesh_size, double *ptr_submesh,  double *ptr_submesh_new, const double *ptr_rhs, double space)
{
    /* jacobi on the interior points */
    for (int i = 1; i< *ptr_to_rows-1 ; i++){
       for ( int j = 1; j< mesh_size-1; j++){      
        /* Update new approximation */
        *(ptr_submesh_new + i *mesh_size +j) =  space * space * (*(ptr_rhs +i * mesh_size +j)) *0.25 + ( *(ptr_submesh + ((i-1)* mesh_size) +j) \
                                                + *(ptr_submesh + ((i+1)*mesh_size) +j) +  *(ptr_submesh + (i *mesh_size) +j -1) + *(ptr_submesh + (i*mesh_size) +j+1) ) * 0.25;
      }     
    }

    /* swap current and new */
   for (int i = 1; i< *ptr_to_rows-1 ; i++){
       for ( int j = 1; j< mesh_size-1; j++){
           
         *(ptr_submesh + i * mesh_size +j ) = *(ptr_submesh_new + i * mesh_size +j);
      
      }     
    }
}

void Jacobi_int(int *ptr_to_rows, int mesh_size, double *ptr_submesh,  double *ptr_submesh_new, const double *ptr_rhs, double space)
{

        /* jacobi on the interior points */
    for (int i = 2; i< *ptr_to_rows-2 ; i++){
       for ( int j = 1; j< mesh_size-1; j++){
           
        /* Update new approximation */
        *(ptr_submesh_new + i *mesh_size +j) =  space * space * (*(ptr_rhs +i * mesh_size +j)) *0.25 + ( *(ptr_submesh + ((i-1)* mesh_size) +j) \
                                                + *(ptr_submesh + ((i+1)*mesh_size) +j) +  *(ptr_submesh + (i *mesh_size) +j -1) + *(ptr_submesh + (i*mesh_size) +j+1) ) * 0.25;

      }     
}
}

void Jacobi_top(int *ptr_to_rows, int mesh_size, double *ptr_submesh,  double *ptr_submesh_new, const double *ptr_rhs, double space)
{
    
    int i = *ptr_to_rows - 2;
    for (int j = 1; j< mesh_size-1; j++ )
         *(ptr_submesh_new + i *mesh_size +j) =  space * space * (*(ptr_rhs +i * mesh_size +j)) *0.25 + ( *(ptr_submesh + ((i-1)* mesh_size) +j) \
                                                + *(ptr_submesh + ((i+1)*mesh_size) +j) +  *(ptr_submesh + (i *mesh_size) +j -1) + *(ptr_submesh + (i*mesh_size) +j+1) ) * 0.25;                                                      
}


void Jacobi_bottom(int *ptr_to_rows, int mesh_size, double *ptr_submesh,  double *ptr_submesh_new, const double *ptr_rhs, double space)
{
    
    int i = 1;

    for (int j = 1; j< mesh_size-1; j++ )
               
         *(ptr_submesh_new + i *mesh_size +j) =  space * space * (*(ptr_rhs +i * mesh_size +j)) *0.25 + ( *(ptr_submesh + ((i-1)* mesh_size) +j) \
                                                + *(ptr_submesh + ((i+1)*mesh_size) +j) +  *(ptr_submesh + (i *mesh_size) +j -1) + *(ptr_submesh + (i*mesh_size) +j+1) ) * 0.25;

                                                       
}