#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"



/*  This function determines new timesteps for the active particles.
 *  If mode==1, or if more than 1/4 of the particles are active,
 *  the timetree is reconstructed from scratch at the end of the 
 *  routine, otherwise particles are one by one removed and reinserted
 *  into the timeline.
 */
void find_timesteps(int mode)    
{    
  double csnd=0, ac, v=0;
  double dt_courant=0, dt=0, dtold, a3inv;
#ifdef VELDISP
  double dt2;
#endif
#ifdef SIDM
  double C_max, dt_sidm;
  double h, hinv, hinv3;
  double C_Grho, dt_Grho; // For Grho timestep

  #if (CROSS_SECTION_TYPE == 2)
  double beta, v_dep;
  #endif
#endif

  int    i;
  int    count;
  double hubble_a=0 ,s_a;
  double t0, t1;
  
  t0=second();



#ifdef SIDM
  C_Grho = BALLINVERSE*(All.DesNumNgb+All.MaxNumNgbDeviation);
#endif

  if(All.ComovingIntegrationOn) 
    {
      hubble_a=All.Hubble*sqrt(All.Omega0/pow(All.Time,3) 
			       + (1-All.Omega0-All.OmegaLambda)/pow(All.Time, 2) + All.OmegaLambda);
      s_a= All.Hubble*sqrt(All.Omega0 + All.Time*(1-All.Omega0-All.OmegaLambda)+pow(All.Time,3)*All.OmegaLambda);
      a3inv= 1/(All.Time*All.Time*All.Time);
#ifdef SIDM

#if (CROSS_SECTION_TYPE == 0)
    C_max = SAFEFACTOR * 
            BALLINVERSE*(All.DesNumNgb+All.MaxNumNgbDeviation)*2*vmax*
            All.CrossSectionInternal/pow(All.Time,2)/s_a;
#elif (CROSS_SECTION_TYPE == 1)
    C_max = SAFEFACTOR * 
            BALLINVERSE*(All.DesNumNgb+All.MaxNumNgbDeviation)*
            All.CrossSectionInternal/pow(All.Time,2.5)/s_a;
#elif (CROSS_SECTION_TYPE == 2)
    vc= All.YukawaVelocity/sqrt(All.Time); /* In internal unit */
    if(2.0*vmax < vc/sqrt(3.0)) {
      beta= 2.0*vmax/vc;
      v_dep= 1.0/(1.0 + beta*beta);
      C_max = SAFEFACTOR * 
	      BALLINVERSE*(All.DesNumNgb+All.MaxNumNgbDeviation)*
	      2.0*vmax*v_dep*v_dep*
              All.CrossSectionInternal/pow(All.Time,2)/s_a;
    }
    else {
      C_max = SAFEFACTOR * 
              BALLINVERSE*(All.DesNumNgb+All.MaxNumNgbDeviation)*
              (3.0*sqrt(3.0)/16.0)*vc*
              All.CrossSectionInternal/pow(All.Time,2)/s_a;
    }
#elif (CROSS_SECTION_TYPE == 3)
    C_max = SAFEFACTOR * 
            BALLINVERSE*(All.DesNumNgb+All.MaxNumNgbDeviation)*
            All.CrossSectionInternal/pow(All.Time,2.0)*2*All.CrossSectionVelScale/s_a;
#elif (CROSS_SECTION_TYPE == 4)
    //vc= All.YukawaVelocity/sqrt(All.Time);
    C_max = SAFEFACTOR * 
            BALLINVERSE*(All.DesNumNgb+All.MaxNumNgbDeviation)*2*vmax*
            All.CrossSectionInternal/pow(All.Time,2)/s_a;
#endif

#endif
    }
  else
    {
      s_a= a3inv= 1;
#ifdef SIDM

#if (CROSS_SECTION_TYPE == 0)
    C_max = SAFEFACTOR * 
            BALLINVERSE*(All.DesNumNgb+All.MaxNumNgbDeviation)*2*vmax*
            All.CrossSectionInternal;
#elif (CROSS_SECTION_TYPE == 1)
    C_max = SAFEFACTOR * 
            BALLINVERSE*(All.DesNumNgb+All.MaxNumNgbDeviation)*
            All.CrossSectionInternal;
#elif (CROSS_SECTION_TYPE == 2)
    vc= All.YukawaVelocity;
    if(2.0*vmax < vc/sqrt(3.0)) {
      v_dep= 1.0/(1.0 + 2.0*vmax/vc);
      C_max = SAFEFACTOR * 
	      BALLINVERSE*(All.DesNumNgb+All.MaxNumNgbDeviation)*
	      2.0*vmax*v_dep*v_dep*
	      All.CrossSectionInternal;
    }
    else {
      C_max = SAFEFACTOR * 
              BALLINVERSE*(All.DesNumNgb+All.MaxNumNgbDeviation)*
              (3.0*sqrt(3.0)/16.0)*vc*
	      All.CrossSectionInternal;
    }
#elif (CROSS_SECTION_TYPE == 3)
    C_max = SAFEFACTOR * 
            BALLINVERSE*(All.DesNumNgb+All.MaxNumNgbDeviation)*2*All.CrossSectionVelScale*
            All.CrossSectionInternal;
#elif (CROSS_SECTION_TYPE == 4)
    //vc= All.YukawaVelocity;
    C_max = SAFEFACTOR * 
            BALLINVERSE*(All.DesNumNgb+All.MaxNumNgbDeviation)*2*vmax*
            All.CrossSectionInternal;

#endif
#endif
    }

  if(NumForceUpdate >= NumPart/4 && mode==0)
    mode=1;

  for(i=IndFirstUpdate,count=0; count<NumForceUpdate; count++, i=P[i].ForceFlag)
    {
      ac= sqrt(P[i].Accel[0]*P[i].Accel[0] + 
	       P[i].Accel[1]*P[i].Accel[1] + 
	       P[i].Accel[2]*P[i].Accel[2]);

      dtold= 2*(P[i].CurrentTime + P[i].MaxPredTime - 2*All.Time);

      if(P[i].Type==0)     /* gas particle */
	{
	  v= sqrt(P[i].VelPred[0]*P[i].VelPred[0] + 
		  P[i].VelPred[1]*P[i].VelPred[1] + 
		  P[i].VelPred[2]*P[i].VelPred[2]);
	  
	  SphP[i].Pressure = GAMMA_MINUS1*SphP[i].EgySpec *SphP[i].Density;
	  csnd  = sqrt(GAMMA*SphP[i].Pressure/SphP[i].Density);
	}

      switch(All.TypeOfTimestepCriterion)
	{
	case 0:
	  dt= sqrt(2*All.ErrTolIntAccuracy*All.SofteningTable[P[i].Type]/ac * s_a);
	  break;
	case 1:
	  dt= All.ErrTolVelScale/ac;
	  break;
#ifdef VELDISP
	case 2:	
	  if(P[i].Type==0)
	    {
	      if(All.ComovingIntegrationOn) 
		dt= All.ErrTolVelScale*(csnd/sqrt(All.Time))/ac;
	      else
		dt= All.ErrTolVelScale*csnd/ac;
	    }
	  else
	    dt= All.ErrTolIntAccuracy*P[i].VelDisp/ac;
	  break;
	case 3:
	  if(P[i].Type==0)
	    {
	      if(All.ComovingIntegrationOn) 
		dt= 3*All.ErrTolIntAccuracy*sqrt(1.5)*hubble_a*All.Time/sqrt(4*M_PI*All.G*SphP[i].Density*a3inv);
	      else
		dt= 3*All.ErrTolIntAccuracy*sqrt(1.5)/sqrt(4*M_PI*All.G*SphP[i].Density);
	    }
	  else
	    {
	      if(All.ComovingIntegrationOn) 
		dt= 3*All.ErrTolIntAccuracy*sqrt(1.5)*hubble_a*All.Time/sqrt(4*M_PI*All.G*P[i].DensVelDisp*a3inv);
	      else
		dt= 3*All.ErrTolIntAccuracy*sqrt(1.5)/sqrt(4*M_PI*All.G*P[i].DensVelDisp);
	    }
	  break;
	case 4:
	  if(All.ComovingIntegrationOn) 
	    {
	      if(P[i].Type==0)
		{
		  dt= 3*All.ErrTolIntAccuracy*sqrt(1.5)*hubble_a*All.Time/sqrt(4*M_PI*All.G*SphP[i].Density*a3inv);
		  dt2= All.ErrTolIntAccuracy*(csnd/sqrt(All.Time))/ac;
		}
	      else
		{
		  dt= 3*All.ErrTolIntAccuracy*sqrt(1.5)*hubble_a*All.Time/sqrt(4*M_PI*All.G*P[i].DensVelDisp*a3inv);
		  dt2= All.ErrTolIntAccuracy*P[i].VelDisp/ac;
		}
	      if(dt2<dt)
		dt=dt2;
	    }
	  else
	    {
	      if(P[i].Type==0)
		{
		  dt= 3*All.ErrTolIntAccuracy*sqrt(1.5)/sqrt(4*M_PI*All.G*SphP[i].Density);
		  dt2= All.ErrTolIntAccuracy*csnd/ac;
		}
	      else
		{
		  dt= 3*All.ErrTolIntAccuracy*sqrt(1.5)/sqrt(4*M_PI*All.G*P[i].DensVelDisp);
		  dt2= All.ErrTolIntAccuracy*P[i].VelDisp/ac;
		}
	      if(dt2<dt)
		dt=dt2;
	    }
	  break;
#endif 
	}
      
	
      if(P[i].Type==0)     /* gas particle */
	{
	  /* need to satisfy the Courant condition for gas particles */

	  if(All.ComovingIntegrationOn) /* comoving variables */
	    {
	      v*= sqrt(All.Time); /* physical velocity */
	      dt_courant = All.CourantFac * All.Time*hubble_a* All.Time*SphP[i].Hsml/
                 (All.Time*SphP[i].Hsml*fabs(sqrt(All.Time)*SphP[i].DivVel) + dmax(csnd,v)*(1+0.6*All.ArtBulkViscConst));
	    }
	  else
	    {
	      dt_courant = All.CourantFac  * SphP[i].Hsml/(SphP[i].Hsml*fabs(SphP[i].DivVel) + 
							     dmax(csnd,v)*(1+0.6*All.ArtBulkViscConst));
	    }

	  if(dt_courant<dt) 
	    dt=dt_courant;
	}
#ifdef SIDM
    else{
      h = P[i].HsmlVelDisp;
      hinv  = 1.0/h;
      hinv3 = hinv*hinv*hinv;
      dt_sidm = All.ProbabilityTol/(C_max*P[i].Mass*hinv3);
      
      if(dt_sidm < dt)
	dt= dt_sidm;

      // Additional G rho timestep criterion sep 23, 2005
      if(All.ComovingIntegrationOn) 
	dt_Grho= All.ErrTolDynamicalAccuracy*hubble_a*All.Time/sqrt(C_Grho*All.G*P[i].Mass*hinv3*a3inv);
      else
	dt_Grho= All.ErrTolDynamicalAccuracy/sqrt(C_Grho*All.G*P[i].Mass*hinv3);

      if(dt_Grho < dt)
	dt= dt_Grho;

    }
#endif

      if(dt>TIMESTEP_INCREASE_FACTOR*dtold)
	{
 	  if(mode!=2)
	    dt= TIMESTEP_INCREASE_FACTOR*dtold;
	}

      if(dt>= All.MaxSizeTimestep)
	{
	  dt= All.MaxSizeTimestep * (1.00 + 0.02*drand48());

	  /* NOTE: It's not good for the ordered binary tree if there
	   * are many particles with EXACTLY the same timestep.
	   * Therefore we randomize it here slightly around the value MaxSizeTimestep.
           * When all particles are constrained by All.MaxSizeTimestep (e.g. at the beginning)
           * they will then actually do all the SAME timestep, of size MaxSizeTimestep.
	   */
	}
	  
      if(dt < All.MinSizeTimestep) 
	{
	  if(P[i].Type==0)
	    {
	      printf("warning: Timestep wants to be below the limit `MinSizeTimestep'\n");
	      printf("Part-ID=%d  dt=%g dtc=%g ac=%g xyz=(%g|%g|%g)  sph: %g %g %g \n", 
		     P[i].ID, dt,  dt_courant, ac, P[i].PosPred[0], P[i].PosPred[1], P[i].PosPred[2],
		     SphP[i].Hsml, SphP[i].EgySpec, csnd); 
	    }
	  else
	    {
	      printf("warning: Timestep wants to be below the limit `MinSizeTimestep'\n");
	      printf("Part-ID=%d  dt=%g ac=%g xyz=(%g|%g|%g)\n", 
		     P[i].ID, dt, ac, P[i].PosPred[0], P[i].PosPred[1], P[i].PosPred[2]);
	    }
	  dt=All.MinSizeTimestep;
	  /* NOTE: It's not good for the ordered binary tree if there
	   * are many particles with EXACTLY the same timestep.
	   * Therefore we randomize it here slightly around the value MaxSizeTimestep.
	   */
	  dt*= 1.0 + 0.02*drand48(); 
	}

      if(mode==0)
        delete_node(i);

      P[i].MaxPredTime = P[i].CurrentTime + 0.5 * dt;

      if(mode==0)
        insert_node(i);
    }

  if(mode>0)
    construct_timetree();

  fflush(stdout);

  t1=second();
  All.CPU_TimeLine+= timediff(t0,t1);

}








