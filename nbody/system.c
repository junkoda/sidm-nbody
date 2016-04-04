#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"


#ifdef WALLCLOCK
#include <mpi.h>
#endif

/* returns the number of cpu-ticks in seconds that
 * have elapsed. (or the wall-clock time)
 */
double second(void)
{
#ifdef WALLCLOCK
  return MPI_Wtime();
#else
  return ((double)clock())/CLOCKS_PER_SEC;
#endif  

  /* note: on AIX and presumably many other 32bit systems, 
   * clock() has only a resolution of 10ms=0.01sec 
   */ 
}

/* returns the time difference between two measurements 
 * obtained with second(). The routine takes care of the 
 * possible overflow of the tick counter on 32bit systems.
 */
double timediff(double t0, double t1)
{
  double dt;
  
  dt=t1-t0;

  if(dt<0)  /* overflow has occured (for systems with 32bit tick counter) */
    {
      dt=t1 + pow(2,32)/CLOCKS_PER_SEC - t0;
    }

  return dt;
}


/* returns the maximum of two double
 */
double dmax(double x,double y)
{
  if(x>y)
    return x;
  else
    return y;
}

/* returns the minimum of two double
 */
double dmin(double x,double y)
{
  if(x<y)
    return x;
  else
    return y;
}



/* returns the maximum of two integers
 */
int imax(int x,int y)
{
  if(x>y)
    return x;
  else
    return y;
}

/* returns the minimum of two integers
 */
int imin(int x,int y)
{
  if(x<y)
    return x;
  else
    return y;
}

