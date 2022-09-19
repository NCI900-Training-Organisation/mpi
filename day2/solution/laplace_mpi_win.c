/* =================================================================
laplace_mpi_win.c

Solve a model 2D Poisson equaton with Dirichlet boundary condition.

-Delta u = 2pi^2 * sin(pi x)sin(pi y) in [0,1]^2
       u = sin(pi x) sin(y) on boundary

The problem is discretised over a uniform mesh by finite difference 
method and the resulting linear system is solved by choices of Jacobi
or Gauss-Seidel.

Compile:  mpicc -g -Wall -O3 -o laplace_mpi_win laplace_mpi_win.c mesh.c solver.c -lm

Usage:  mpirun -np 4 ./laplace_mpi_win size max_iter method

Produced for NCI Training. 

Frederick Fung 2022
4527FD1D

Please leave comments at frederick.fung@anu.edu.au
====================================================================*/
#include<stdio.h>
#include <stdlib.h>
#include<string.h>
#include<math.h>
#include<mpi.h>
#include "mesh.h"
#include "solver.h"

#define MPIIO
//#define MPI_DEBUG

int main(int argc, char *argv[]){

int rank, cells; 
MPI_Init(&argc, &argv);

MPI_Comm world = MPI_COMM_WORLD;
MPI_Comm_rank(world, &rank);

MPI_Comm_size(world, &cells);

int mesh_size, max_iter; 

double space;

char *method;

/* build MPI_Datatype arg */
int blocklength[3] = {1, 1, 1};
MPI_Datatype type_list[3] = {MPI_INT, MPI_INT, MPI_DOUBLE};
MPI_Aint displacements[3];
MPI_Aint A_mesh_size, A_max_iter, A_space;
MPI_Datatype arg;

MPI_Get_address(&mesh_size, &A_mesh_size);
MPI_Get_address(&max_iter, &A_max_iter);
MPI_Get_address(&space, &A_space);

displacements[2] = A_space - A_max_iter;
displacements[1] = A_max_iter - A_mesh_size;
displacements[0] = 0;

MPI_Type_create_struct(3, blocklength, displacements, type_list, &arg);
MPI_Type_commit(&arg);



/* parse arguments on process 0 */
if (rank == 0){
    if (argc == 4){

    mesh_size = atof(argv[1]);

    max_iter = atof(argv[2]);

    method = argv[3];
    
    /* Broadcast args to the rest of processes */
    MPI_Bcast(&mesh_size, 1, arg, 0, world);
     
    if ((strcmp(method, "Jacobi") == 0)){
        printf("%s METHOD IS IN USE   \n ", method); }
    else {    
        printf( "Not a valid method\n");
        MPI_Finalize();
        exit(1);
         }
    }
    else {
        printf("Usage: %s [size] [max_iter] [method] \n", argv[0]);
        MPI_Finalize();
        exit(1);
    }
}
else {
    MPI_Bcast(&mesh_size, 1, arg, 0, world);
}

/* grid spacing */
space = (double) 1 / (mesh_size-1);

/* number of interior rows in each cell */
int int_rows = (mesh_size -2) / cells ;

/* calc the remaining rows */
int extra_rows = (mesh_size -2 ) - int_rows * cells;

/* total number of rows per cell, adding top and bottom ghost rows */
int rows = int_rows + 2;

/* add extra rows */
int rows_top = rows + extra_rows; /* one for ghost, one for bnd */

int *ptr_rows = NULL;

/* assign extra rows to the last cell */ 
if (rank == (cells - 1)) ptr_rows = &rows_top;

else ptr_rows = &rows;

if (rows <= 3) MPI_Abort(world, 1); /* add error handling */ 


/* alloc mem for meshes held in each cell */
double (*submesh)[mesh_size] = malloc(sizeof *submesh * *ptr_rows);

/* jacobi method requires to store updates */ 
double (*submesh_new)[mesh_size] = malloc(sizeof *submesh_new * *ptr_rows);

/* alloc mem for rhs held in each cell */
double (*subrhs)[mesh_size] = malloc(sizeof *subrhs * *ptr_rows);


init_mesh(mesh_size, submesh, submesh_new, subrhs, rank, cells, int_rows, space, ptr_rows);


/* window object for top and bottom bnd data */
MPI_Win bnd_win;

/* allocate mem for the win obj */
double *bnd;
MPI_Alloc_mem(2*sizeof(double )*mesh_size, MPI_INFO_NULL, &bnd);

/* create local memory window */
MPI_Win_create(bnd, 2*mesh_size *sizeof(double), sizeof(double),\
               MPI_INFO_NULL, world, &bnd_win);


/* Assign topology to the ranks */
int higher = rank +1;
if (higher >= cells) higher = MPI_PROC_NULL;
int lower = rank -1;
if (lower < 0) lower = MPI_PROC_NULL;

unsigned iter = 0;
while (iter<max_iter)
{
    iter+=1;

    /* as the target, sync with the origin process after local access completes */
    MPI_Win_fence(0, bnd_win);

    MPI_Put(submesh[1], mesh_size, MPI_DOUBLE, lower, 0, mesh_size, MPI_DOUBLE, bnd_win);
    MPI_Put(submesh[*ptr_rows - 2], mesh_size, MPI_DOUBLE, higher, mesh_size, mesh_size, MPI_DOUBLE, bnd_win);   

    /* as the origin, sync with the target to ensure that put completes */
    MPI_Win_fence(0, bnd_win);

    /* local window access */
    if (rank > 0 && rank < (cells -1)){
        for (int i=0; i< mesh_size; i++){
            submesh[0][i] = bnd[mesh_size+i]; /* bottom */
            submesh[*ptr_rows-1][i] = bnd[i]; /* top */            
        }
    }
    else if (rank ==0 ){ /* bottom rank only updates its top full row */
        for (int i = 0; i< mesh_size; i++ ){
            submesh[*ptr_rows-1][i] = bnd[i];
        }
    }
    else if (rank == (cells -1) ){ /* top rank only updates its bottom full row */ 
        for (int i = 0; i <mesh_size; i++){
            submesh[0][i] = bnd[mesh_size+i];
        }
    }

    Jacobi(ptr_rows, mesh_size, &submesh[0][0], &submesh_new[0][0], &subrhs[0][0], space);
}

MPI_Win_free(&bnd_win);

/* sync after solving the problem on each cell */
int highertag=1, lowertag=2;
MPI_Send(submesh[1], mesh_size, MPI_DOUBLE, lower, lowertag, world);
MPI_Recv(submesh[*ptr_rows -1], mesh_size, MPI_DOUBLE, higher, lowertag, world, MPI_STATUS_IGNORE);
MPI_Send(submesh[*ptr_rows-2], mesh_size, MPI_DOUBLE, higher, highertag, world);
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

