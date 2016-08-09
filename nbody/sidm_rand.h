#ifdef SIDM


#ifdef RANDOM_GSL
  #include <gsl/gsl_rng.h>
  gsl_rng* rand_gen;
#ifdef INLINE
 static inline double ran2(long* idum)
  {
    return gsl_rng_uniform(rand_gen);
  }
#else
  double ran2(long* idum);
#endif

#else
  double ran2(long* idum);
#endif

extern long iseed;
void init_rand(int seed, int ThisTask);

#ifdef INLINE
static inline void random_direction(double n[])
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
#else
void random_direction(double n[]);
#endif

#endif





