#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"



/* This routine computes the accelerations for all active particles. 
 * First, the gravitational forces are computed (this also reconstructs
 * the tree, if needed). Also note that the gas-particle tree will
 * in any case be updated in its geometrical properties.
 *
 * If gas particles are presented, the `interior' of the local domain
 * is determined. This region is guaranteed to contain only particles
 * local to the processor. This information will be used to reduce
 * communication in the hydro part.
 * The the density for active SPH particles is computed. If the number
 * of neighbours should be outside the allowed bounds, it will be
 * readjusted by the function ensure_neighbours(). 
 * Finally, the hydrodynamical forces are added.
 */
void compute_accelerations(int mode) 
{
  double tstart,tend;

  if(ThisTask==0)
    {
      printf("Start force computation...\n"); 
      fflush(stdout);
    }

  tstart=second();   /* measure the time for the full force computation */

  gravity_tree();    /* computes gravity accel. */

  tend=second();
  All.CPU_Gravity+= timediff(tstart,tend);


#ifdef VELDISP
  tstart=second();
  determine_interior();
  veldisp();
  veldisp_ensure_neighbours(mode); 
    // checks that number of neighbours is correct 
  tend=second();
  All.CPU_EnsureNgb+= timediff(tstart, tend);
#endif

  /* Wrong place for scattering. Sep 27, 2005 */
  // return Oct 9, 2005

#ifdef SIDM
  tstart=second();
  determine_interior();
  
  if(mode == 0){
    sidm();
    sidm_ensure_neighbours(mode);
  }
  
  tend=second();
  All.CPU_EnsureNgb+= timediff(tstart, tend);
#endif

  if(All.TotN_gas > 0) 
    {
      if(ThisTask==0)
	{
	  printf("Start density computation...\n"); 
	  fflush(stdout);
	}


      tstart=second();
#ifndef VELDISP
#ifndef SIDM
      determine_interior();
#endif
#endif
      density();   /* computes density, and pressure */      
    
      tend=second();
      All.CPU_Hydro+= timediff(tstart, tend);      
     
      
      tstart=second();
    
      ensure_neighbours(mode); /* checks that number of neighbours is correct */

      tend=second();
      All.CPU_Hydro+= timediff(tstart, tend); 
      All.CPU_EnsureNgb+= timediff(tstart, tend);

      if(ThisTask==0)
	{
	  printf("Start hydro-force computation...\n"); 
	  fflush(stdout);
	}

      tstart=second();

      hydro_force();      /* adds hydrodynamical accelerations 
			       and computes du/dt  */
#ifdef COOLING
      cooling_and_starformation();   /* do radiative cooling and star formation */
#endif

      tend=second();
      All.CPU_Hydro+= timediff(tstart, tend);
    }

  if(ThisTask==0)
    {
      printf("force computation done.\n"); 
      fflush(stdout);
    }

#ifdef DEB
  if(ThisTask==0)
    { 	
      FdDEB=fopen("accel.deb","w");
      fprintf(FdDEB," %d \n",All.NumCurrentTiStep);
      fclose(FdDEB);
    }
#endif
}























