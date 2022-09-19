/* =================================================================
solver.h

The header for solver.c

Preepared for NCI Training. 

Frederick Fung 2022
4527FD1D
====================================================================*/

double local_L2_residual(const int *ptr_to_rows, 
                         int mesh_size, 
                         double space, 
                         const double *restrict ptr_submesh, 
                         const double *restrict ptr_subrhs);

void Jacobi(int *ptr_to_rows, 
            int mesh_size, 
            double *restrict ptr_submesh,  
            double *restrict ptr_submesh_new, 
            const double *restrict ptr_rhs, 
            double space);

void Jacobi_int(int *ptr_to_rows, 
                int mesh_size, 
                double *restrict ptr_submesh,  
                double *restrict ptr_submesh_new, 
                const double *restrict ptr_rhs, 
                double space);

void Jacobi_top(int *ptr_to_rows, 
                int mesh_size, 
                double *restrict ptr_submesh,  
                double *restrict ptr_submesh_new, 
                const double *restrict ptr_rhs, 
                double space);


void Jacobi_bottom(int *ptr_to_rows, 
                   int mesh_size, 
                   double *restrict ptr_submesh,  
                   double *restrict ptr_submesh_new, 
                   const double *restrict ptr_rhs, 
                   double space);