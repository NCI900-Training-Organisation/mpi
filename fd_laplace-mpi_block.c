/* =================================================================
fd_laplace-serial.c

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
    printf("mem alloc %d\n",*ptr_rows);

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
    printf("rank, size  %d, %d\n", rank, cells);
    
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


double local_l2_residual(const int *ptr_to_rows, int mesh_size, double space, const double *ptr_submesh, const double *ptr_subrhs){
    
    /* l2 norm of residual by a given approximation */
    double residual = 0.0f;

      for (int i = 1; i< *ptr_to_rows -1; i++){
        for (int j = 1; j< mesh_size-1; j++){

             
             double diff =  ( 4* *(ptr_submesh +  (i * mesh_size) +j ) - *(ptr_submesh + ((i-1)*mesh_size) +j) - *(ptr_submesh +((i+1)*mesh_size)+j)\
                             - *(ptr_submesh + (i *mesh_size) + j-1 ) - *(ptr_submesh + (i*mesh_size) +j +1)) / (space *space )- *(ptr_subrhs + i*mesh_size +j);
      
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

int main(int argc, char *argv[]){

int rank, cells; 

MPI_Init(&argc, &argv);

MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  printf("rank init %d\n", rank);

MPI_Comm_size(MPI_COMM_WORLD, &cells);
printf("size init %d\n", cells);

int mesh_size; 

double space;

double tolerance;

char *method;


/* parse arguments */
if (argc == 4){
mesh_size = atof(argv[1]);

tolerance = atof(argv[2]);

method = argv[3];

if ((strcmp(method, "Gauss-Seidel") == 0)){
   printf("%s METHOD IS IN USE   \n ", method);}
else if ((strcmp(method, "Jacobi") == 0)){
    printf("%s METHOD IS IN USE \n", method);
}
else {
    
    printf( "Not a valid method\n");
    exit(1);
}


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

printf("ptr size %d\n", *ptr_rows);

printf("size %d\n", rows);

if (rows <= 3) MPI_Abort(MPI_COMM_WORLD, 1); /* add error handling */ 


/* alloc mem for meshes held in each cell */
double (*submesh)[mesh_size] = malloc(sizeof *submesh * *ptr_rows);

/* jacobi method requires to store updates */ 
double (*submesh_new)[mesh_size] = malloc(sizeof *submesh_new * *ptr_rows);

/* alloc mem for rhs held in each cell */
double (*subrhs)[mesh_size] = malloc(sizeof *subrhs * *ptr_rows);


init_mesh(mesh_size, submesh, submesh_new, subrhs, rank, cells, int_rows, space, ptr_rows);


int uppertag=1, lowertag=2;

MPI_Status status;


/* Assign topology to the ranks */
int upper = rank +1;
if (upper >= cells) upper = MPI_PROC_NULL;
int lower = rank -1;
if (lower < 0) lower = MPI_PROC_NULL;



for (int iter = 0; iter < 100; iter++)
{

/* send bottom buffer to lower */
MPI_Send(submesh[1], mesh_size, MPI_DOUBLE, lower, lowertag, MPI_COMM_WORLD);
/* receive buffer to allocate to the top ghost layer */
MPI_Recv(submesh[*ptr_rows -1], mesh_size, MPI_DOUBLE, upper, lowertag, MPI_COMM_WORLD, &status);


#ifdef MPI_DEBUG
        printf("MPI process %d received value from rank %d, with tag %d and error code %d.\n", rank, status.MPI_SOURCE, status.MPI_TAG, status.MPI_ERROR);

#endif

/* send top buffer to upper */
MPI_Send(submesh[*ptr_rows-2], mesh_size, MPI_DOUBLE, upper, uppertag, MPI_COMM_WORLD);
/* receive buffer to allocate to the bottom ghost layer */
MPI_Recv(submesh[0], mesh_size, MPI_DOUBLE, lower, uppertag, MPI_COMM_WORLD, &status);

#ifdef MPI_DEBUG
        printf("MPI process %d received value from rank %d, with tag %d and error code %d.\n", rank, status.MPI_SOURCE, status.MPI_TAG, status.MPI_ERROR);

#endif

Jacobi(ptr_rows, mesh_size, &submesh[0][0], &submesh_new[0][0], &subrhs[0][0], space);



     
double residual, tot_res;
     residual  = local_l2_residual(ptr_rows, mesh_size, space, &submesh[0][0], &subrhs[0][0]);
    
   
MPI_Reduce(&residual, &tot_res, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   //printf("Reduce");
   if (rank == 0){

       tot_res = sqrt(tot_res);
        
       printf("Residual1  %f\n",  tot_res); }
}

/* sync after solving the problem on each cell */
MPI_Send(submesh[1], mesh_size, MPI_DOUBLE, lower, lowertag, MPI_COMM_WORLD);
MPI_Recv(submesh[*ptr_rows -1], mesh_size, MPI_DOUBLE, upper, lowertag, MPI_COMM_WORLD, &status);
MPI_Send(submesh[*ptr_rows-2], mesh_size, MPI_DOUBLE, upper, uppertag, MPI_COMM_WORLD);
MPI_Recv(submesh[0], mesh_size, MPI_DOUBLE, lower, uppertag, MPI_COMM_WORLD, &status);

/* calc residual */
double residual, tot_res;
residual  = local_l2_residual(ptr_rows, mesh_size, space, &submesh[0][0], &subrhs[0][0]);
    

/* collecting residuals and returns to rank 0 */
MPI_Reduce(&residual, &tot_res, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   
if (rank == 0){

    tot_res = sqrt(tot_res);
        
    printf("Residual1  %f\n",  tot_res); }

       


/* Output data */
char file_name[20];

#ifndef MPIIO
sprintf(file_name, "laplace-soln-rank-%d", rank);
FILE *fp = fopen(file_name, "w");
for (int i = 0; i<*ptr_to_rows; i++){
    for (int j = 0; j<mesh_size; j++){
        fwrite(&submesh[i][j], sizeof(double), 1, fp);
    }
}
fclose(fp);
#endif 


#ifdef MPIIO
sprintf(file_name, "laplace-soln-whole");

MPI_File fp;
MPI_Offset offset; 
MPI_Status IO_status;
MPI_File_open(MPI_COMM_WORLD, file_name, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
int count, w_count ;

if (rank ==(cells -1)){

    count = mesh_size * (*ptr_rows -1);
    offset=  (mesh_size * (int_rows +1 )) * sizeof(double) + (rank-1)*(mesh_size * int_rows) * sizeof(double);
    MPI_File_write_at(fp, offset, &submesh[1][0], count, MPI_DOUBLE, &IO_status );

    #ifdef MPI_DEBUG
    MPI_Get_count(&IO_status, MPI_DOUBLE, &w_count);
    printf("w_count for cell %d is %d from offset %lld\n", rank, w_count, offset);
    #endif

}

else if (rank == 0 ) {

    count =  mesh_size * (*ptr_rows -1);
    offset = 0;
    MPI_File_write_at(fp, offset, &submesh[0][0], count, MPI_DOUBLE, &IO_status);


    #ifdef MPI_DEBUG
    MPI_Get_count(&IO_status, MPI_DOUBLE, &w_count);
    printf("w_count for cell %d is %d from offset %lld\n", rank, w_count, offset);
    #endif MPI_DEBUG
}

else{

    count = (mesh_size) * int_rows; /* count number of full node */
    offset =  (mesh_size * (*ptr_rows -1 )) * sizeof(double) + (rank-1)*(mesh_size * (*ptr_rows -2 )) * sizeof(double);
    MPI_File_write_at(fp, offset, &submesh[1][0], count, MPI_DOUBLE, &IO_status );

    #ifdef MPI_DEBUG
    MPI_Get_count(&IO_status, MPI_DOUBLE, &w_count);
    printf("w_count for cell %d is %d from offset %lld\n", rank, w_count, offset);
    #endif
}

MPI_File_close(&fp);
#endif


MPI_Finalize();
free(submesh);
free(submesh_new);
free(subrhs);


}

else {
    printf("Usage: %s [size] [tolerance] [method] \n", argv[0]);
    exit(1);
}

}