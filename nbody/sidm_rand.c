#ifdef SIDM

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#ifdef RANDOM_GSL
gsl_rng* rand_gen;
#endif
extern long   iseed;

#include "sidm_rand.h"




void init_rand(int seed, int ThisTask)
{
  /* Init rand() */
#ifdef RANDOM_GSL
  int i;
  const gsl_rng_type *rt;
  gsl_rng_env_setup();

  rt= gsl_rng_mt19937; //gsl_rng_default;
  rand_gen= gsl_rng_alloc(rt);
  //unsigned long seed= seed1 + seed2*ThisTask; //All.Seed1 + All.Seed2*ThisTask;
  gsl_rng_set(rand_gen, seed);
  for(i=0; i<1000000; i++)
    ran2(&iseed);

  if(ThisTask == 0)
    fprintf(stdout, "GSL mt19937 Random Number Initialized %d: %f\n",
	    ThisTask, ran2(&iseed));
  else
    ran2(&iseed);

  //  gsl_rng_free(r);  
#else
  iseed= -seed; // //All.Seed1 + All.Seed2*ThisTask);
  if(ThisTask == 0)
    fprintf(stdout, "Random Number Initialized %d: %f\n",ThisTask,ran2(&iseed));
  else
    ran2(&iseed);
#endif
}

#ifndef RANDOM_GSL

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-15
#define RNMX (1.0-EPS)

double ran2(long *idum)
{
  int j;
  long k;
  static long idum2= 123456789;
  static long iy= 0;
  static long iv[NTAB];
  float temp;

  if(*idum <= 0){
    if(-(*idum) < 1)
      *idum= 1;
    else
      *idum= -(*idum);
    idum2= (*idum);

    for(j=NTAB+7 ; j>=0 ; j--){
      k= (*idum)/IQ1;
      *idum= IA1*(*idum-k*IQ1) - k*IR1;
      if(*idum < 0)
    *idum += IM1;
      if(j < NTAB)
    iv[j]= *idum;
    }
    iy= iv[0];
  }
  k= (*idum)/IQ1;
  *idum= IA1*(*idum-k*IQ1)-k*IR1;
  if(*idum < 0)
    *idum += IM1;
  k= idum2/IQ2;
  idum2= IA2*(idum2-k*IQ2)-k*IR2;
  if(idum2 < 0)
    idum2 += IM2;
  j= iy/NDIV;
  iy= iv[j]-idum2;
  iv[j]= *idum;
  if(iy < 1)
    iy += IMM1;
  if((temp= AM*iy) > RNMX)
    return RNMX;
  
  return temp;
}
#endif

#ifndef INLINE
static void random_direction(double n[])
{
  double y1, y2, r2, sq1r;
  do{
    y1= 1.0-2.0*ran2(&iseed);
    y2= 1.0-2.0*ran2(&iseed);
    r2= y1*y1 + y2*y2;
  } while(r2 > 1.0);

  sq1r= sqrt(1.0-r2);
  n[0]= 2.0*y1*sq1r;
  n[1]= 2.0*y2*sq1r;
  n[2]= 1.0-2.0*r2;
}

#ifdef RANDOM_GSL
double ran2(long* idum)
{
  return gsl_rng_uniform(rand_gen);
}
#endif

#endif



#endif
