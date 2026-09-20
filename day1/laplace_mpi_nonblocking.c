/* =================================================================
laplace_mpi_nonblocking.c

Solve a model 2D Poisson equaton with Dirichlet boundary condition.

-Delta u = 2pi^2 * sin(pi x)sin(pi y) in [0,1]^2
       u = sin(pi x) sin(pi y) on boundary

The problem is discretised over a uniform mesh by finite difference 
method and the resulting linear system is solved by Jacobi iteration.


Compile:  mpicc -g -Wall -O3 -o laplace_mpi_nonblocking laplace-mpi_nonblocking.c -lm
        mesh.c solver.c

Usage:  mpirun -np 4 ./laplace_mpi_nonblocking mesh_size max_iter Jacobi

Prepared for NCI Training. 

Frederick Fung 2022
4527FD1D

Please leave comments at frederick.fung@anu.edu.au
====================================================================*/
#include<stdio.h>
#include <stdlib.h>
#include<string.h>
#include<math.h>
#include<mpi.h>
#include <errno.h>
#include <limits.h>
#include "mesh.h"
#include "solver.h"

//#define MPIIO
//#define MPI_DEBUG

int main(int argc, char *argv[]){

int rank, cells; 
MPI_Init(&argc, &argv);

MPI_Comm world = MPI_COMM_WORLD;
MPI_Comm_rank(world, &rank);
MPI_Comm_size(world, &cells);

/* Rank 0 validates input; all ranks take the same success or failure path. */
int args[2] = {0, 0}; /* mesh_size, max_iter */
int valid = 1;
if (rank == 0) {
    valid = argc == 4 && strcmp(argv[3], "Jacobi") == 0;
    if (valid) {
        for (int i = 0; i < 2; i++) {
            char *end;
            errno = 0;
            long value = strtol(argv[i + 1], &end, 10);
            if (errno != 0 || end == argv[i + 1] || *end != '\0' ||
                value < 0 || value > INT_MAX) {
                valid = 0;
                break;
            }
            args[i] = (int)value;
        }
        /* The split Jacobi update requires at least two owned rows per rank. */
        valid = valid && args[0] >= 4 && (args[0] - 2) / cells >= 2;
    }
    if (!valid) {
        fprintf(stderr, "Usage: %s [mesh_size] [max_iter] Jacobi\n"
                "Use integer values, max_iter >= 0, and at least two interior rows per rank.\n",
                argv[0]);
    }
}
MPI_Bcast(&valid, 1, MPI_INT, 0, world);
if (!valid) {
    MPI_Finalize();
    return EXIT_FAILURE;
}
MPI_Bcast(args, 2, MPI_INT, 0, world);
int mesh_size = args[0], max_iter = args[1];
double space = 1.0 / (mesh_size - 1);
if (rank == 0) printf("Jacobi METHOD IS IN USE\n");

/* number of interior rows in each process */
int int_rows = (mesh_size -2) / cells ;

/* calc the remaining rows */
int extra_rows = (mesh_size -2 ) - int_rows * cells;

/* total number of rows per cell, adding top and bottom ghost rows */
int rows = int_rows + 2;

/* add top extra rows */
int rows_top = rows + extra_rows; 
int *ptr_rows = NULL;

/* assign extra rows to the last cell */ 
if (rank == (cells - 1)) ptr_rows = &rows_top;
else ptr_rows = &rows;

/* alloc mem for meshes held in each cell */
double (*submesh)[mesh_size] = malloc(sizeof *submesh * *ptr_rows);

/* jacobi method requires to store a copy for updates */ 
double (*submesh_new)[mesh_size] = malloc(sizeof *submesh_new * *ptr_rows);

/* alloc mem for rhs held in each cell */
double (*subrhs)[mesh_size] = malloc(sizeof *subrhs * *ptr_rows);
if (submesh == NULL || submesh_new == NULL || subrhs == NULL) {
    fprintf(stderr, "Mesh allocation failed on rank %d\n", rank);
    MPI_Abort(world, EXIT_FAILURE);
}


/* setup mesh config */
init_mesh(mesh_size, submesh, submesh_new, subrhs, rank, cells, int_rows, space, ptr_rows);

int highertag=1, lowertag=2;

#TODO: binds both top and bottom bnd comm requests to one single arrary
MPI_Status top_bnd_status[2], bottom_bnd_status[2];
MPI_Request top_bnd_requests[2],  bottom_bnd_requests[2];
int top_flag, bottom_flag;

/* Assign topology to the ranks */
int upper = rank +1;
if (upper >= cells) upper = MPI_PROC_NULL;
int lower = rank -1;
if (lower < 0) lower = MPI_PROC_NULL;

int iter = 0;
while (iter< max_iter)
{
    iter+=1;
   
    /* communicate to the higher rank process */
    MPI_Irecv(submesh[*ptr_rows -1], mesh_size, MPI_DOUBLE, upper, highertag, world, &top_bnd_requests[0]);
    MPI_Isend(submesh[1], mesh_size, MPI_DOUBLE, lower, highertag, world, &bottom_bnd_requests[1]);

    /* communicate to the lower rank process */
    MPI_Irecv(submesh[0], mesh_size, MPI_DOUBLE, lower, lowertag, world, &bottom_bnd_requests[0]);
    MPI_Isend(submesh[*ptr_rows-2], mesh_size, MPI_DOUBLE, upper, lowertag, world, &top_bnd_requests[1]);

    /* Jacobi Interior */
    Jacobi_int(ptr_rows, mesh_size, &submesh[0][0], &submesh_new[0][0], &subrhs[0][0], space);

    #TODO: wait on all requests altogether
    /* Test on either the top or bottom layer */
    if ( (MPI_Testall(2, top_bnd_requests, &top_flag, top_bnd_status) > 0) || (MPI_Testall(2, bottom_bnd_requests, &bottom_flag, bottom_bnd_status) > 0))
    {
       MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /* if the top layer is ready */
    if (top_flag){
        /* perform jacobi on the top bnd */
        Jacobi_top(ptr_rows, mesh_size, &submesh[0][0], &submesh_new[0][0], &subrhs[0][0], space);

    /* if the the bottom layer is ready */
    if (bottom_flag){
    /* perform jacobi on the bottom bnd */
    Jacobi_bottom(ptr_rows, mesh_size, &submesh[0][0], &submesh_new[0][0], &subrhs[0][0], space);
    }
    
    /* if the bottom layer is yet ready */
    else{
        /* wait on the bottom layer */
        MPI_Waitall(2, bottom_bnd_requests, bottom_bnd_status);
        Jacobi_bottom(ptr_rows, mesh_size, &submesh[0][0], &submesh_new[0][0], &subrhs[0][0], space);        
        }
    }
    /* if the top layer is yet ready but the buttom is */
    else if (bottom_flag){
        /* perform jacobi bottom ready */
        Jacobi_bottom(ptr_rows, mesh_size, &submesh[0][0], &submesh_new[0][0], &subrhs[0][0], space);
        /* wait on the top layer */
        MPI_Waitall(2, top_bnd_requests, top_bnd_status);
        Jacobi_top(ptr_rows, mesh_size, &submesh[0][0], &submesh_new[0][0], &subrhs[0][0], space);
        }
    /* if neither of the top and bottom is ready, then wait on both */
    else {
        MPI_Waitall(2, bottom_bnd_requests, bottom_bnd_status);
        Jacobi_bottom(ptr_rows, mesh_size, &submesh[0][0], &submesh_new[0][0], &subrhs[0][0], space);
        MPI_Waitall(2, top_bnd_requests, top_bnd_status);
        Jacobi_top(ptr_rows, mesh_size, &submesh[0][0], &submesh_new[0][0], &subrhs[0][0], space);
    }

    /* swap current with new approx */
    for (int i = 1; i< *ptr_rows-1 ; i++){
        for ( int j = 1; j< mesh_size-1; j++){
                *(&submesh[0][0]+ i * mesh_size +j ) = *(&submesh_new[0][0]+ i * mesh_size +j);
        }     
    }
}

/* sync after solving the problem on each cell */
MPI_Send(submesh[1], mesh_size, MPI_DOUBLE, lower, lowertag, world);
MPI_Recv(submesh[*ptr_rows -1], mesh_size, MPI_DOUBLE, upper, lowertag, world, MPI_STATUS_IGNORE);
MPI_Send(submesh[*ptr_rows-2], mesh_size, MPI_DOUBLE, upper, highertag, world);
MPI_Recv(submesh[0], mesh_size, MPI_DOUBLE, lower, highertag, world, MPI_STATUS_IGNORE);

/* calc residual */
double residual, tot_res;
residual  = local_L2_residual(ptr_rows, mesh_size, space, &submesh[0][0], &subrhs[0][0]);
    
/* collecting residuals and returns to rank 0 */
MPI_Reduce(&residual, &tot_res, 1, MPI_DOUBLE, MPI_SUM, 0, world);   
if (rank == 0){
    tot_res = sqrt(tot_res);        
    printf("Final Residual %f after %d iterations.\n",  tot_res, max_iter); 
}

MPI_Finalize();
free(submesh);
free(submesh_new);
free(subrhs);

}