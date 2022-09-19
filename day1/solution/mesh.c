/* =================================================================
mesh.c

The mesh configuration for a unit square. It assumes the mesh is divided
vertically to the submeshes.

Produced for NCI Training. 

Frederick Fung 2022
4527FD1D
====================================================================*/


#include<stdio.h>
#include <stdlib.h>
#include<string.h>
#include<math.h>
#include "mesh.h"

#define MPIIO
//#define MPI_DEBUG


double bnd_fc(int x, int y, double space){
    /* Boundary function */

    double value = sin(M_PI * x * space) * sin(M_PI *y *space);
    //printf ("value %f\n", value);
    return value;
}

double rhs_fc(int x, int y, double space){
    /* rhs sin(pi*x)sin(pi*y) */
    double value =(2.0* M_PI*M_PI) * sin( M_PI * x * space) * sin ( M_PI * y * space); 
    return value;
}


void init_mesh(const int mesh_size,
               double submesh[][mesh_size], 
               double submesh_new[][mesh_size],
               double subrhs[][mesh_size],
               const int rank,
               const int cells,
               const int int_rows,
               const double space,
               int *ptr_rows)
{

/* relocate mem for the last rank (top slab) as it includes the extra_size work */
if (rank == (cells -1) ) {
    double (*tmp1)[mesh_size] = realloc(submesh, sizeof *submesh * *ptr_rows);
    double (*tmp2)[mesh_size] = realloc(submesh_new, sizeof *submesh_new * *ptr_rows);
    double (*tmp3)[mesh_size] = realloc(subrhs, sizeof *submesh * *ptr_rows);

    if (tmp1 == NULL || tmp2 == NULL || tmp3 == NULL){

        printf("Reallocation fails");
        exit(0);
    }

    else{

        submesh =tmp1;
        submesh_new = tmp2;
        subrhs = tmp3;
    }
    
    double x_coord, y_coord;

    /* initialise interior submesh values on the last rank */
    for (int i = 1; i< *ptr_rows-1; i++){
        for (int j = 0; j < mesh_size; j++){
           x_coord = j;
           y_coord = rank * int_rows + i;

           /* initialise side boundaris */
           if ( j == 0 || j == mesh_size-1 ){
               submesh[i][j] = bnd_fc(x_coord, y_coord, space); 
               subrhs[i][j] = 0.0;
           }
           
           /* iniitialise interior and rhs */
           else {
           submesh[i][j] = 0.0;
           subrhs[i][j] = rhs_fc(x_coord, y_coord, space);
           }
        }
    }
    for (int i = 0; i< mesh_size; i++) {   
        /* bottom ghost row */
        submesh[0][i] =0.0;
        /* top boundary row */
        submesh[*ptr_rows-1][i] =bnd_fc(i, mesh_size-1, space); 
    }
    
}

/* initialise shubmesh values on rank 0 */
else if (rank ==0 ) {
    double x_coord, y_coord;
    /* for interior rows */
    for (int i = 1; i< *ptr_rows-1; i++ ){
        for (int j =0; j< mesh_size ; j++){
            x_coord = j;
            y_coord = rank * int_rows +i ;
            
            /* initialise side boundaries */ 
            if ( j ==0 || j == mesh_size-1){
                submesh[i][j] = bnd_fc(x_coord, y_coord, space);
                subrhs[i][j] = 0.0;
            }

            /* initialise interior submesh values */
            else{
            submesh[i][j] = 0.0;
            subrhs[i][j] = rhs_fc(x_coord, y_coord, space);
            }            
        }
    }
    
    for (int i = 0; i<mesh_size; i++){     
        /* bottom boundary row */
        submesh[0][i] = bnd_fc(i, 0, space);
        /* top ghost tow */
        submesh[*ptr_rows-1][i] = 0.0;
    }
}

/* initialise middle cells */
else if (rank >0 && rank <cells -1){

    double x_coord, y_coord;

    /* for interior rows */
    for (int i = 1; i< *ptr_rows -1; i++){
        for (int j = 0; j< mesh_size; j++){

            x_coord =j;
            y_coord = rank * int_rows + i;

            /* initialise side boundaries */
            if (j == 0 || j == mesh_size-1 ){
                submesh[i][j] = bnd_fc(x_coord, y_coord, space);
                subrhs[i][j] = 0.0;
            }
            
            /* initialise interior submesh values */
            else{
            submesh[i][j] = 0.0;
            subrhs[i][j] = rhs_fc(x_coord, y_coord, space);
        }
    }
   }

   /* top and bottom ghost rows */
   for (int i=0; i< mesh_size; i++){
       submesh[0][i] = 0.0;
       submesh[*ptr_rows-1][i] = 0.0;
   }
}
}


