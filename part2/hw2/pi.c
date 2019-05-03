#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#define f(x)   ( 4.0e0 / (1.0e0 + (x)*(x)) )

double mysecond();

int main(int argc, char *argv[])
{

  double PI25DT = 3.141592653589793238462643e0;

  double  pi, h, sum, x, t0,t1;
  int     n, i, ierr;
 
  char line[12];

  scanf("%d",&n);
  
  h    = 1.0e0/n;           /* Calculate the interval size */
  sum  = 0.0e0;
  
  t0 = mysecond();
  for(i = 1; i <= n; i++)
    {
      x = h * ( (double)(i) - 0.5e0 );
      sum = sum + f(x);
    }
  t1 = mysecond();

  pi = h * sum;
  
  printf(" calc. pi:%20.16f  Error:%20.16f  %13.9f(sec)\n", pi, pi - PI25DT, t1-t0 );

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
