/* ==============================================================

The program simulates molecular dynamics governed by the Lennard-
Jones potential. 

The system is initialised in a 2D unit square with uniform distributed
point mass moledulars. The temperature T is set to be degree 1, and 
the velocities are sampled with Gaussian distribution with standard
deviation sqrt(T/2), mean zero.

The velocity verlet method is used to update the positions at each
time instance.
 

Compile gcc -g -Wall -O3 -lm -o md-serial md-serial.c 

Frederick Fung Apr 2022

================================================================ */



#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

#define N 10
#define rcut2 1.25992

typedef struct molecular *Mo;


/* Define molecular */

struct molecular{
    /* relative position in a local cartesian coord. */
    int id;
    double r_x;
    double r_y;
    double mass;
    double v_x;
    double v_y;
    double f_x;
    double f_y;
};


Mo create_molecular(int id, double r_x, double r_y, double mass,
                    double v_x, double v_y, double f_x, double f_y){

   
     struct molecular *m = malloc(sizeof(*m));
     
     m->id= id;
     m->r_x = r_x;
     m->r_y = r_y;
     m->mass = mass;
     m->v_x = v_x;
     m->v_y = v_y;
     m->f_x = f_x;
     m->f_y = f_y;
   printf("id %d\n", id);
     return m;


}

void destroy_molecular(Mo m){
    free(m);
};


/* Gaussian Random by Box Muller  */
double rand_normal(double mean, double stddev)
{
    static double n2 = 0.0;
    static int n2_cached = 0;
    if (!n2_cached)
    {
        double x, y, r;
        do
        {
            x = 2.0*rand()/RAND_MAX - 1;
            y = 2.0*rand()/RAND_MAX - 1;

            r = x*x + y*y;
        }
        while (r == 0.0 || r > 1.0);
        {
            double d = sqrt(-2.0*log(r)/r);
            double n1 = x*d;
            n2 = y*d;
            double result = n1*stddev + mean;
            n2_cached = 1;
            return result;
        }
    }
    else
    {
        n2_cached = 0;
        return n2*stddev + mean;
    }
}



void initialise_MD(Mo **MD){

   //printf("pointer address %p \n", (void* )&(MD));

   *MD = malloc(N*N * sizeof(MD));
   for (int i = 0; i<N; ++i)
    { 
    for (int j = 0; j<N; ++j ){
       
        double v_x = rand_normal(0.0, sqrt(0.5));
        double v_y = rand_normal(0.0, sqrt(0.5));
        
        /* uniform distribution of moleculars in a square domain */
        double r_x = i/ (double) N;
        double r_y = j/ (double) N;

        int id =  i*N +j;
        double f_x = 0.0;
        double f_y = 0.0;
        double mass = 1.0;
 
        (*MD)[id] = create_molecular(id, r_x, r_y, mass, v_x, v_y, f_x, f_y);
    
   
    }

    }

}

/* calculate force between each pair of moleculars */
void force_on_each(Mo *MD){
     
    for (int i = 0;  i< N*N ; i++)
        for(int j = 0; j < N*N; j++){
           
            
            double rdiff_x = fabs(MD[i]->r_x - MD[j]->r_x);
            double rdiff_y = fabs(MD[i]->r_y - MD[j]->r_y);
            double r2 = rdiff_x * rdiff_x + rdiff_y * rdiff_y;
           

            if (r2 < rcut2 && r2>0.0){
           
                double r2inv = 1/r2;
                double r6inv = r2inv * r2inv * r2inv;
           
                MD[i]->f_x += 24 * r6inv * r2inv * (2*r6inv + 1) * rdiff_x;
           
                MD[i]->f_y += 24 * r6inv * r2inv * (2*r6inv + 1) * rdiff_y;

            }
        }

}

/* use velocity verlet integration to update positions and velocities at each time instance */
void verlet_method(Mo *MD, int ts, double step){
     

   for (int i =0; i< ts ; i++){
       force_on_each(MD);
     for (int i = 0; i < N*N; i++){
         
        double r_old_x, r_old_y, r_new_x, r_new_y;

        if (i == 0){
            double r_old_x = MD[i]->r_x;
            double r_old_y = MD[i]->r_y;
        }
        /* under assumption mass = 1, a = f_x */ 
        r_new_x = 2.0 * MD[i]->r_x - r_old_x + step * step * MD[i]->f_x ;

        r_new_y = 2.0 * MD[i]->r_y - r_old_y + step * step * MD[i]->f_y ;

        MD[i]->v_x = (r_new_x - r_old_x) / (2.0 * step); 

        MD[i]->v_y = (r_new_y - r_old_y) / (2.0 * step);

        r_old_x = MD[i] -> r_x;

        r_old_y = MD[i] -> r_y;

        MD[i]->r_x = r_new_x;

        MD[i]->r_y = r_new_y;
     }

   }
}



int main (){

    srand(time(NULL));

    Mo *MD =NULL ;
    
    int ts = 10;

    double step = 0.01;

    initialise_MD(&MD);

    verlet_method(MD, ts, step);

    for (int i =0; i< N*N; i++){
        destroy_molecular(MD[i]);
    }

    *MD =NULL; 

}