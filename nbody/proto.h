#ifndef INLINE_FUNC
#ifdef INLINE
#define INLINE_FUNC inline
#else
#define INLINE_FUNC
#endif
#endif

#include "forcetree.h"
#include "cooling.h"

void   advance(void);
void   allocate_commbuffers(void);
void   allocate_memory(void);
void   begrun(void);
void   check_omega(void);
void   close_outputfiles(void);
void   compute_accelerations(int mode);
void   compute_global_quantities_of_system(void);
void   compute_potential(void);
void   construct_timetree(void);
void   cooling_and_starformation(void);
void   delete_node(int i);
void   density(void);
void   determine_interior(void);
int    dissolvegas(void);
double dmax(double,double);
double dmin(double,double);
void   do_box_wrapping(void);
void   DomainDecomposition(void);
double drand48();
void   endrun(int);
void   energy_statistics(void);
void   ensure_neighbours(int mode);
void   every_timestep_stuff(void);
void   ewald_corr(double dx, double dy, double dz, double *fper); // INLINE
void   ewald_force(double x[3], double force[3]);
void   ewald_init(void);
double ewald_psi(double x[3]);
double ewald_pot_corr(double dx, double dy, double dz); // INLINE
int    find_ancestor(int i);
double find_next_outputtime(double time);
void   find_next_time(void);
int    find_next_time_walk(int node);
void   free_memory(void);
void   find_timesteps(int mode);   
double get_starformation_rate(int i);
void   gravity_tree(void);
void   hydro_force(void);
int    imax(int,int);
int    imin(int,int);
void   init(void);
void   insert_node(int i);
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE *stream);
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE *stream);
void   open_outputfiles(void);
#ifndef INLINE
double INLINE_FUNC periodic(double x);
#else
inline double periodic(double x) // jk Feb 26
{
  while(x > All.BoxHalf)
    x -=All.BoxSize;

  while(x < -All.BoxHalf)
    x+=All.BoxSize;

  return x;
}
#endif
double pot_integrand(double xx);
void   predict(double time);
void   predict_collisionless_only(double time);
void   predict_sph_particles(double time);
void   read_ic(char *fname);
void   read_ic_cluster(char *fname);
void   read_ic_cluster_gas(char *fname);
void   read_ic_cluster_wimp(char *fname);
int    read_outputlist(char *fname);
void   read_parameter_file(char *fname);
void   restart(int mod);
void   run(void);
void   savepositions(int num);
void   savepositions_ioformat1(int num);
double second(void);
void   set_softenings(void);
void   set_sph_kernel(void);
void   set_units(void);
void   setup_smoothinglengths(int desired_ngb);
void   setup_smoothinglengths_veldisp(int desired_ngb);
void   statistics(void);
double timediff(double t0,double t1);
void   veldisp(void);
void   veldisp_ensure_neighbours(int mode);


/* on some DEC Alphas, the correct prototype for pow() is missing,
   even when math.h is included ! */

double pow(double, double);

void   reflect(void);
