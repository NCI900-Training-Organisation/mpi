/* =================================================================
fd_laplace-serial.c

Solve a model 2D Poisson equaton with Dirichlet boundary condition.

-Delta u = 2pi^2 * sin(pi x)sin(pi y) in [0,1]^2
       u = sin(pi x) sin(y) on boundary

The problem is discretised over a uniform mesh by finite difference 
method and the resulting linear system is solved by choices of Jacobi
or Gauss-Seidel.


Compile:  gcc -fopenmp -g -Wall -O3 -lm -o fd_laplace-serial  fd_laplace-serial.c

Usage:  ./fd_laplace-omp size tolerance method

Produced for NCI Training. 

Frederick Fung 2022
4527FD1D
====================================================================*/


#include<stdio.h>
#include <stdlib.h>
#include<string.h>
#include<math.h>

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


double l2_residual(int size, double space, double grid[size][size], double rhs[size][size]){
    
    /* l2 norm of residual by a given approximation */
    double residual = 0.0f;

    for (int i = 1; i < size -1 ; i++ ){
        for (int j = 1; j < size - 1; j++){
            
            /* residual = Ax^k -b */
            double diff = (4*grid[i][j]-grid[i-1][j]-grid[i+1][j]-grid[i][j+1]-grid[i][j-1])/ (space *space) -rhs[i][j];
           //  printf("i %d, j %d, IJ, %f, I-1J, %f, I+1J %f, IJ+1 %f, IJ-1 %f, stencil %f diff %f\n", i, j, grid[i][j], grid[i-1][j], grid[i+1][j], grid[i][j+1], grid[i][j-1],  4*grid[i][j]-grid[i-1][j]-grid[i+1][j]-grid[i][j+1]-grid[i][j-1],  diff);
            
            diff = space*space*diff;
            residual += pow(diff, 2);


        }
      }
    
    
    return sqrt(residual);

    
}


void Jacobi(double tolerance, int size, double grid[size][size], double rhs[size][size], double space){
    
    double (*grid_new)[size] = malloc(sizeof *rhs *size);

    //double residual = 1.0f;
    double *p_grid, *p_grid_tmp, *p_rhs, *swap;
    int iter = 0;
    int i, j;

    //printf("1 starting address of grid %p\n", (void*)&grid[0][0]);
    memcpy(&grid_new[0][0], &grid[0][0], size *size *sizeof(grid[0][0]));
    p_grid = &grid[0][0];
    p_grid_tmp = &grid_new[0][0];
    p_rhs = &rhs[0][0];
    
    /* iterates until tolerance is reached */
    do {
    for (i = 1; i< size -1; i++){    
        for ( j = 1; j< size -1; j++){
        
        /* Update new approximation */
        //*(p_grid_tmp + i *size +j) =  space * space * (*(p_rhs +i * size +j)) *0.25 + ( *(p_grid + ((i-1)* size) +j) + *(p_grid + ((i+1)*size) +j)
        //                                                                     /+  *(p_grid + (i *size) +j -1) + *(p_grid + (i*size) +j+1) ) * 0.25;
        grid_new[i][j] = space  *space * rhs[i][j] * 0.25 + ( grid[i-1][j]+ grid[i+1][j] + grid[i][j-1]+ grid[i][j+1])*0.25;
        //printf("Space %f Jacobi i %d, j %d, subij %f, IJ, %f, I-1J, %f, I+1J %f, IJ+1 %f, IJ-1 %f,  new %f\n", space , i, j, rhs[i][j], grid[i][j], grid[i-1][j], -grid[i+1][j], grid[i][j+1], grid[i][j-1],  \
                space  *space * rhs[i][j] * 0.25 + ( grid[i-1][j]+ grid[i+1][j] + grid[i][j-1]+ grid[i][j+1])*0.25);

         }     
    }

      //printf("2 p_grid  %p, p_grid_tmp %p\n", (void*)p_grid, (void*)p_grid_tmp);
      //printf("2 starting address of submesh %p, of submesh_tmp %p\n", (void*)&grid[0][0], (void*)&grid_tmp[0][0]);
       //  swap = p_grid_tmp;
       //  p_grid_tmp = p_grid;
        // p_grid = swap;
         
             /* swap current and new */
   for (int i = 1; i< size-1 ; i++){
       for ( int j = 1; j< size-1; j++){
           
        /* Update new approximation */
       // printf("update on rank %d %d-th row and %d-th col\n", rank, i, j);
         grid[i][j] = grid_new[i][j];
      
      }     
    }

    //printf("3 p_grid  %p, p_grid_tmp %p\n", (void*)p_grid, (void*)p_grid_tmp); 
    //printf("3 starting address of submesh %p, of submesh_tmp %p\n", (void*)&grid[0][0], (void*)&grid_tmp[0][0]);

         
         
        if (iter % 1 == 0){
         double residual2 = l2_residual(size, space, grid, rhs);
        
         printf("Residual2 after %d iteration: %f\n", iter, residual2);
        }
        iter+=1;
        
         

    }
   while (iter <1000 );
   //while (residual > tolerance);
}

void GS(double tolerance, int size, double grid[size][size], double rhs[size][size], double space){


    double residual = 1.0f;
    double *p_grid,  *p_rhs;
    int iter = 0;
    int i, j;


    
    p_grid = &grid[0][0];
    p_rhs = &rhs[0][0];
        
    /* iterates until tolerance is reached */
    do {
   
    iter +=1;
    
    for (i = 1; i< size -1; i++){
       for ( j = 1; j< size -1; j++){
           
        /* Update new approximation */
        *(p_grid + i *size +j) =  space * space * (*(p_rhs +i * size +j)) *0.25 + ( *(p_grid + ((i-1)* size) +j) + *(p_grid + ((i+1)*size) +j)
                                                                              +  *(p_grid + (i *size) +j -1) + *(p_grid + (i*size) +j+1) ) * 0.25;
      }     
    }
     

        /* Print out residuals after every 100 iterations */
        if (iter % 1 == 0){
        residual = l2_residual(size, space, grid, rhs);
        
        printf("Residual after %d iteration: %.10f\n", iter, residual);
        }


    }
   
   while (residual > tolerance);
}

int main(int argc, char *argv[]){


int size;

double space;

double tolerance;

char *method;


/* parse arguments */
if (argc == 4){
size = atof(argv[1]);

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


FILE *fp;
//clrscr();
fp=fopen("laplace-soln","w"); //output will be stored in this file

space = (double) 1 / (size-1);

double (*grid)[size] = malloc(sizeof *grid * size);

double (*rhs)[size] = malloc(sizeof *rhs * size);




/* Initialise meshes */


for (int i = 0; i< size; i++){
    for (int j = 0 ; j< size; j++){
        
        if (i == 0 || j == 0 || i == size-1 || j == size - 1 ){
            
            grid[i][j] = bnd_fc(i, j,  space);
            rhs[i][j] = 0.0;
         
        }
        else
       {
        grid[i][j] = 0.0;
        rhs[i][j] = rhs_fc(i,j, space);
       }
    }
    

      
    }


/* Smoothers */
if (strcmp(method, "Gauss-Seidel")==0)
GS(tolerance, size, grid, rhs, space);
else{
    double residual1 = l2_residual(size, space, grid, rhs);
    printf("Residual1 after %f\n", residual1);
Jacobi(tolerance, size, grid, rhs, space);}
//double residual1 = l2_residual(size, space, grid, rhs);
        
//printf("Residual1 after %f\n", residual1);

char data[10]; /* 8 bytes for double, 1 byte for space and 1 byte for null */
int re;
/* Output data  */
for (int i = 0; i< size; i++){
  for (int j = 0; j<size; j++)
               //fprintf(fp,"%.5f\t",grid[i][j]);
            //fprintf(fp,"\n");
            {
          //  re = sprintf(data, "%6f ", grid[i][j]); 
           // printf("re %d\n", re);

            fwrite(&grid[i][j], sizeof(grid[0][0]), 1, fp);
            }
            //fwrite("\n", sizeof(char), 1, fp);

}

fclose(fp);
exit(0);
}
else {
    printf("Usage: %s [size] [tolerance] [method] \n", argv[0]);
    exit(1);
}

}