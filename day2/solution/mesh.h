/* =================================================================
mesh.h

The header of mesh.c

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


