/* =================================================================
solver.c

This file stores routines for numerical solvers. Currently it has 
jacobi method and a routine that calculates residual based on 5-point
stencils.

Prepared for NCI Training. 

Frederick Fung 2022
4527FD1D
====================================================================*/
#include<stdio.h>
#include <stdlib.h>
#include<string.h>
#include<math.h>
#include<mpi.h>
#include "solver.h"


double local_L2_residual(const int *ptr_to_rows, int mesh_size, double space, const double *restrict ptr_submesh, const double *restrict ptr_subrhs)
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


void Jacobi(int *ptr_to_rows, int mesh_size, double *restrict ptr_submesh,  double *restrict ptr_submesh_new,  const double *restrict ptr_rhs, double space)
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



void Jacobi_int(int *ptr_to_rows, int mesh_size, double *restrict ptr_submesh,  double *restrict ptr_submesh_new,  const double *restrict ptr_rhs, double space)
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

void Jacobi_top(int *ptr_to_rows, int mesh_size,  double *restrict ptr_submesh, double *restrict ptr_submesh_new,  const double *restrict ptr_rhs, double space)
{
    
    int i = *ptr_to_rows - 2;
    for (int j = 1; j< mesh_size-1; j++ )
         *(ptr_submesh_new + i *mesh_size +j) =  space * space * (*(ptr_rhs +i * mesh_size +j)) *0.25 + ( *(ptr_submesh + ((i-1)* mesh_size) +j) \
                                                + *(ptr_submesh + ((i+1)*mesh_size) +j) +  *(ptr_submesh + (i *mesh_size) +j -1) + *(ptr_submesh + (i*mesh_size) +j+1) ) * 0.25;                                                      
}


void Jacobi_bottom(int *ptr_to_rows, int mesh_size,  double *restrict ptr_submesh,  double *restrict ptr_submesh_new,  const double *restrict ptr_rhs, double space)
{
    
    int i = 1;

    for (int j = 1; j< mesh_size-1; j++ )
               
         *(ptr_submesh_new + i *mesh_size +j) =  space * space * (*(ptr_rhs +i * mesh_size +j)) *0.25 + ( *(ptr_submesh + ((i-1)* mesh_size) +j) \
                                                + *(ptr_submesh + ((i+1)*mesh_size) +j) +  *(ptr_submesh + (i *mesh_size) +j -1) + *(ptr_submesh + (i*mesh_size) +j+1) ) * 0.25;

                                                       
}