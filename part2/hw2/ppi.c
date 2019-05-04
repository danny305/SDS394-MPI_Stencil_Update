#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "mpi.h"
#include <math.h>


#define f(x)   ( 4.0e0 / (1.0e0 + (x)*(x)) )

double mysecond();

int main(int argc, char *argv[]){

	MPI_Init(&argc, &argv);

	
  double PI25DT = 3.141592653589793238462643e0;

  double  pi, h, sum, sum_rank, x;
	double t0_int,t1_int, t0_red, t1_red,t_diff_int, t_diff_red;
	double min_int_t, max_int_t, avg_int_t, min_red_t, max_red_t, avg_red_t;
  int     n, i, ierr;
	int rank, nranks;
	MPI_Comm comm = MPI_COMM_WORLD;


	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nranks);
 
  char line[12];

	// Read in total number of intervals
	if (rank == 0){
		  scanf("%d \n",&n);
  		printf("\nScanned: %d \n", n);
	}


	//Send total number fo intervals to all tasks
	ierr = MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (rank ==0){
		printf("Sent broadcast from root! \n");
	}


	MPI_Barrier(MPI_COMM_WORLD);


  h    = 1.0e0/n;           /* Calculate the interval size */
  sum  = 0.0e0;
  
  t0_int = mysecond();
  for(i = rank + 1; i <= n; i += nranks)
    {
      x = h * ( (double)(i) - 0.5e0 );
      sum = sum + f(x);
    }
  t1_int = mysecond();

  sum_rank = h * sum;
  
	MPI_Barrier(MPI_COMM_WORLD);

	//Get pi from all tasks
	t0_red;	
	MPI_Reduce(&sum_rank, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);  
	t1_red;

	//Get time integration avg, min, max
	t_diff_int= t1_int - t0_int;
	MPI_Reduce(&t_diff_int, &avg_int_t, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&t_diff_int, &min_int_t, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&t_diff_int, &max_int_t, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	avg_int_t /= nranks;


	//Get reduction time avg, min, max
	t_diff_red = t1_red - t0_red;
	MPI_Reduce(&t_diff_red, &avg_red_t, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&t_diff_red, &min_red_t, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&t_diff_red, &max_red_t, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	avg_red_t /= nranks;

	if (rank == 0){
		printf("calc. pi:%20.16f  Error:%20.16f  %13.9f(sec)\n", pi, pi - PI25DT, t_diff_int );
		printf("AVG TIMES (number of ranks - %d): \n"\
					 "\tIntegration time: %13.9f(sec)  Reduction time: %13.9f(sec)\n",
					  nranks, avg_int_t, avg_red_t );
		printf("MIN TIMES (number of ranks - %d): \n"\
					 "\tIntegration time: %13.9f(sec)  Reduction time: %13.9f(sec)\n",
					  nranks, min_int_t, min_red_t );
		printf("MAX TIMES (number of ranks - %d): \n"\
					 "\tIntegration time: %13.9f(sec)  Reduction time: %13.9f(sec)\n",
					  nranks, max_int_t, max_red_t );
	}


	ierr = MPI_Finalize();
  return(0);
}

double mysecond()
{
   struct timeval tp;
   struct timezone tzp;
   int i;
   i = gettimeofday(&tp,&tzp);
   return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}


/* http://www-unix.mcs.anl.gov/mpi/www/ */
