#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"



/* Calculate hydro force and rate of change of internal energy.
 * Also, we do the cooling here in an isochoric approximation,
 * with an implicit integration.
 */
void hydro_force(void)
{
  double  h_i, h_i2, hinv, hinv4;
  double  p_over_rho2_i, p_over_rho2_j, soundspeed_i, soundspeed_j;
  double  r, r2, u, mass_j, dt;
  double  dx,dy,dz;
  double  hfc, dwk_i, vdotr, vdotr2, visc, mu_ij, c_ij, rho_ij, h_ij, f1, f2;
  int     i,j,k,ii,n, numngb;
  double  h_j,dwk_j;
  int     jj;
  double  hubble_a=0, s_a_inverse, sqrt_of_a=0, prefac=0, hfc_egy, fac_vsic_fix=0, a3inv=0; 
  float   *r2list;
  int     *ngblist;
  double  tstart,tend;
  int     ntot,ntotleft,npleft,nthis;     
  int     nstart,nbuffer,ncount,nchunk;
  int     timelinecounter,startcounter;   
  int     *nrecv,*noffset;
  int     level,sendTask,recvTask;
  int     nexport, place;
  MPI_Status status;


  if(All.ComovingIntegrationOn)
    {
      /* Factors for comoving integration of hydro */

      sqrt_of_a= sqrt(All.Time);
      a3inv=  1/(All.Time*All.Time*All.Time);
      hubble_a= All.Hubble*sqrt(All.Omega0/(All.Time*All.Time*All.Time) 
			     + (1-All.Omega0-All.OmegaLambda)/(All.Time*All.Time) + All.OmegaLambda);
      s_a_inverse=1/(All.Hubble*sqrt(All.Omega0 + All.Time*(1-All.Omega0-All.OmegaLambda)+
				     (All.Time*All.Time*All.Time)*All.OmegaLambda));
      prefac= s_a_inverse/All.Time;
      fac_vsic_fix= hubble_a * All.Time*All.Time*All.Time;
    }

  /* `NumForceUpdate' gives the number of particles on this processor that want a force update */  

  MPI_Allreduce(&NumSphUpdate, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(NumSphUpdate==0)
    nthis= 1;
  else
    nthis = ((double)NumSphUpdate*(All.BunchSizeHydro-NTask))/ntot + 1; /* maximum size of chunck coming from this task */

  /* Note: it can happen that this processor does have NumForceUpdate==0, while ntot>0 */
    
  if(NumSphUpdate>0) 
    {
      nchunk= NumSphUpdate/nthis;                             /* number of chunks needed */  
      if((NumSphUpdate%nthis) > 0) 
	nchunk+=1;
    }
  else 
    nchunk=0;

  nrecv=   malloc(sizeof(int)*NTask);  /* list of particle numbers that constituate current bunch */
  noffset= malloc(sizeof(int)*NTask);  /* offsets of bunches in common list */

  ntotleft= ntot;              /* particles left for all tasks together */
  npleft= NumSphUpdate;        /* particles left for this task */

  nstart= IndFirstUpdate;      /* first particle for this task */
  startcounter=0;



  while(ntotleft>0)
    {
      if(nthis>npleft)
	nthis=npleft;

      for(i=nstart, ncount=0, nexport=0, timelinecounter=startcounter; 
	  timelinecounter<NumForceUpdate && ncount<nthis; 
	  i=P[i].ForceFlag, timelinecounter++)
	if(P[i].Type==0)
	  {
	    ncount++;

	    for(j=0; j<3; j++)
	      {
		if(P[i].PosPred[j] < (DomainMin[0][j]+ 4*SphP[i].Hsml))
		  break;
		if(P[i].PosPred[j] > (DomainMax[0][j]- 4*SphP[i].Hsml))
		  break;
	      }

	    if(j!=3)  /* particle lies NOT completely inside . needs to be sent to other processors */
	      {
		nexport++;
                P[i].Type |= 8;  /* flag the particle for export */
	      }
	  }

      MPI_Allgather(&nexport, 1, MPI_INT, nrecv, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Allreduce(&ncount, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      for(i=nbuffer=0; i<NTask; i++)  /* compute length of common list */
	nbuffer += nrecv[i];

      for(i=1,noffset[0]=0; i<NTask; i++) /* set-up offset-table */
	noffset[i]= noffset[i-1]+nrecv[i-1];

     /* fill in the particles at the right place */
      
      for(i=nstart, ncount=0, nexport=0, timelinecounter=startcounter; 
	  timelinecounter<NumForceUpdate && ncount<nthis; 
	  i=P[i].ForceFlag, timelinecounter++)
	{
	  if((P[i].Type&7)==0)
	    {
	      if(P[i].Type&8)
		{
		  place= noffset[ThisTask] + nexport; 
                  nexport++;
		}
	      else
		place= nbuffer + (ncount-nexport);

	      for(k=0;k<3;k++)
		{
		  HydroDataIn[place].Pos[k]= P[i].PosPred[k];
		  HydroDataIn[place].Vel[k]= P[i].VelPred[k];
		}
	      HydroDataIn[place].Hsml        = SphP[i].Hsml;
	      HydroDataIn[place].Mass        = P[i].Mass;
#ifdef SFR
	      HydroDataIn[place].FormedStellarMass = SphP[i].FormedStellarMass;
#ifdef CLOUDS
	      HydroDataIn[place].CloudMass =         SphP[i].CloudMass;
#endif
#endif
	      HydroDataIn[place].DensityPred = SphP[i].DensityPred;
	      HydroDataIn[place].Pressure    = SphP[i].Pressure;
	      HydroDataIn[place].Current     = P[i].CurrentTime;

	      /* calculation of F1 */
	      if(SphP[i].DensityPred > 0)
		{
		  p_over_rho2_i = SphP[i].Pressure/(SphP[i].DensityPred * SphP[i].DensityPred);
		  soundspeed_i  = sqrt(GAMMA*p_over_rho2_i*SphP[i].DensityPred);
		  f1 = fabs(SphP[i].DivVel)/
		    (fabs(SphP[i].DivVel)+SphP[i].CurlVel + 0.0001*soundspeed_i/SphP[i].Hsml);
		}
	      else
		f1=0;
	      HydroDataIn[place].F1 = f1;
  	      
	      ncount++;
	    }
	}
 

      /* now start big communication */

      tstart=second();      
      for(level=1;level<NTask;level++)
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ level; 

	  MPI_Sendrecv(&HydroDataIn[noffset[sendTask]],  nrecv[sendTask]*sizeof(struct hydrodata_in), MPI_BYTE, recvTask, TAG_ANY, 
	               &HydroDataIn[noffset[recvTask]],  nrecv[recvTask]*sizeof(struct hydrodata_in), MPI_BYTE, recvTask, TAG_ANY, MPI_COMM_WORLD, &status);
	}
      tend=second();
      All.CPU_CommSum+= timediff(tstart, tend);


      

      /* up to here, the procedure is quite similar to density(). Now force evaluation. */

      /* Main loop for hydro acceleration and adiabatic du/dt */

      for(i=0; i<(nbuffer+ncount-nexport); i++)  /* Loop over all particles in the buffer */
	{
	  HydroDataResult[i].DtEgySpec=0;
	  HydroDataResult[i].Acc[0]= HydroDataResult[i].Acc[1]= HydroDataResult[i].Acc[2]= 0;

	  if(HydroDataIn[i].DensityPred>0)  
	    {
	      p_over_rho2_i = HydroDataIn[i].Pressure/(HydroDataIn[i].DensityPred * HydroDataIn[i].DensityPred);
	      soundspeed_i  = sqrt(GAMMA*p_over_rho2_i*HydroDataIn[i].DensityPred);
	    }
	  else
	    p_over_rho2_i=soundspeed_i=0;

	  h_i = HydroDataIn[i].Hsml;
	  h_i2= h_i*h_i;


	  numngb= ngb_treefind_pairs(&HydroDataIn[i].Pos[0], HydroDataIn[i].Hsml, &ngblist, &r2list);
	  
	  for(n=0; n<numngb; n++)
	    {
	      j = ngblist[n]+1; 
	      r2= r2list[n];
	      h_j = SphP[j].Hsml;
	      
	      if(r2<h_i2 || r2<h_j*h_j)
		{
		  r = sqrt(r2);		      

		  if(r>0)
		    {		  
		      if(SphP[j].DensityPred>0)
			{
			  p_over_rho2_j = SphP[j].Pressure/(SphP[j].DensityPred*SphP[j].DensityPred);
			  soundspeed_j  = sqrt(GAMMA*p_over_rho2_j*SphP[j].DensityPred);
			}
		      else
			p_over_rho2_j=soundspeed_j=0;
		      
		      dx = HydroDataIn[i].Pos[0] - P[j].PosPred[0];
		      dy = HydroDataIn[i].Pos[1] - P[j].PosPred[1];
		      dz = HydroDataIn[i].Pos[2] - P[j].PosPred[2];
#ifdef PERIODIC
		      dx= periodic(dx);
		      dy= periodic(dy);
		      dz= periodic(dz);
#endif
		      vdotr = ( dx * (HydroDataIn[i].Vel[0] - P[j].VelPred[0])
			      + dy * (HydroDataIn[i].Vel[1] - P[j].VelPred[1])
			      + dz * (HydroDataIn[i].Vel[2] - P[j].VelPred[2]) );

		      if(All.ComovingIntegrationOn)
			{
			  vdotr2 = vdotr/sqrt_of_a + hubble_a*r2;
			}
		      else
			vdotr2 = vdotr; 
			  
		      
		      if(r2<h_i2)
			{
			  hinv = 1.0/h_i;
			  hinv4 = hinv*hinv*hinv*hinv;
			  
			  u = r*hinv;
			  ii = (int)(u*KERNEL_TABLE);
			  dwk_i= hinv4*( KernelDer[ii] + (KernelDer[ii+1]-KernelDer[ii])*(u-KernelRad[ii])*KERNEL_TABLE);
			}
		      else
			dwk_i=0;
		      
		      
		      if(r2<h_j*h_j)
			{
			  hinv = 1.0/h_j;
			  hinv4 = hinv*hinv*hinv*hinv;
			  
			  u = r*hinv;
			  jj = (int)(u*KERNEL_TABLE);
			  dwk_j= hinv4*( KernelDer[jj] + (KernelDer[jj+1]-KernelDer[jj])*(u-KernelRad[jj])*KERNEL_TABLE);
			}
		      else
			dwk_j=0;

		      if(vdotr2<0)   /* ... artificial viscosity */		      
			{
			  c_ij  = 0.5*(soundspeed_i + soundspeed_j);
			  h_ij  = 0.5*(h_i + h_j);
			  
			  if(All.ComovingIntegrationOn)
			    mu_ij  = All.Time * h_ij*vdotr2/(r2+0.01*h_ij*h_ij);
			  else
			    mu_ij  = h_ij*vdotr2/(r2+0.01*h_ij*h_ij);
			  
			  rho_ij = 0.5*(HydroDataIn[i].DensityPred + SphP[j].DensityPred);
			  
			  f1 = HydroDataIn[i].F1;
			  f2 = fabs(SphP[j].DivVel)/(fabs(SphP[j].DivVel) + SphP[j].CurlVel + 0.0001*soundspeed_j/SphP[j].Hsml);
			  
			  if(rho_ij>0)
			    visc = (-All.ArtBulkViscConst*mu_ij*c_ij + 2*All.ArtBulkViscConst*mu_ij*mu_ij)/rho_ij*(f1+f2)*0.5;
			  else
			    visc=0;
			  
			  /* .... end artificial viscosity evaluation */
			  /* now make sure that viscous acceleration is not too large */
			  
			  dt = 2*(All.Time - HydroDataIn[i].Current); 

			  if(dt>0 && (dwk_i + dwk_j)<0)
			    {
			      if(All.ComovingIntegrationOn) /* comoving variables */
				visc = dmin(visc, fac_vsic_fix* vdotr2/
					    (0.5*(HydroDataIn[i].Mass+P[j].Mass)*(dwk_i + dwk_j)*r*dt)); 
			      else
				visc = dmin(visc,  vdotr2/
				            (0.5*(HydroDataIn[i].Mass+P[j].Mass)*(dwk_i + dwk_j)*r*dt)); 
			    }
			}
		      else
			visc=0;

		      /* calculate final acceleration and rate of change of thermal energy */
		      
#ifdef SFR
#ifdef CLOUDS
		      mass_j = P[j].Mass - SphP[j].FormedStellarMass- SphP[j].CloudMass;
#else
		      mass_j = P[j].Mass - SphP[j].FormedStellarMass;
#endif
#else
		      mass_j=  P[j].Mass;
#endif		  
		      if(All.ComovingIntegrationOn) /* comoving variables */
			{
			  /* arithemtic mean for symmetrization */
			  /*
			  hfc  = prefac* 0.5*mass_j*(p_over_rho2_i + p_over_rho2_j + visc)*(dwk_i+dwk_j)/r;
			  */
			  /* or geometric mean. Both can be chosen */
			  hfc  = prefac* 0.5*mass_j*(2*sqrt(p_over_rho2_i*p_over_rho2_j) + visc)*(dwk_i+dwk_j)/r;

			  hfc_egy  = hfc*All.Time*sqrt_of_a;
			}
		      else
			{
			  /* arithemtic mean for symmetrization */
			  /*
			  hfc      = 0.5*mass_j*(p_over_rho2_i + p_over_rho2_j + visc)*(dwk_i+dwk_j)/r;
			  */
			  /* or geometric mean. Both can be chosen */
			  hfc      = 0.5*mass_j*(2*sqrt(p_over_rho2_i*p_over_rho2_j) + visc)*(dwk_i+dwk_j)/r;

			  hfc_egy  = hfc; 
			}

#ifdef SFR
#ifdef CLOUDS
		      hfc*= (HydroDataIn[i].Mass - HydroDataIn[i].FormedStellarMass - HydroDataIn[i].CloudMass)/
                             HydroDataIn[i].Mass;
#else
		      hfc*= (HydroDataIn[i].Mass - HydroDataIn[i].FormedStellarMass)/
 			     HydroDataIn[i].Mass;
#endif
#endif
		      
		      HydroDataResult[i].Acc[0] -= hfc*dx;
		      HydroDataResult[i].Acc[1] -= hfc*dy;
		      HydroDataResult[i].Acc[2] -= hfc*dz;
		      
		      HydroDataResult[i].DtEgySpec += 0.5 * hfc_egy * vdotr2;
		    }
		}
	    } 
	}
 


      
      /* communicate the results and sum them up */
      tstart=second();
      for(level=1; level<NTask; level++)
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ level; 

	  MPI_Sendrecv(&HydroDataResult[noffset[recvTask]], nrecv[recvTask]*sizeof(struct hydrodata_out), MPI_BYTE, recvTask, TAG_ANY, 
	               &HydroDataPartialResult[0],          nrecv[sendTask]*sizeof(struct hydrodata_out), MPI_BYTE, recvTask, TAG_ANY, MPI_COMM_WORLD, &status);

	  for(i=0; i<nrecv[ThisTask]; i++)
	    {
	      HydroDataResult[noffset[ThisTask] + i].DtEgySpec += HydroDataPartialResult[i].DtEgySpec;

	      for(k=0;k<3;k++)
		HydroDataResult[noffset[ThisTask] + i].Acc[k] += HydroDataPartialResult[i].Acc[k];
	    }
	}
      tend=second();
      All.CPU_CommSum+= timediff(tstart, tend);
	 




      /* transfer the result to the SPH-particles on this processor */

      for(i=nstart, ncount=0, nexport=0, timelinecounter=startcounter; 
	  timelinecounter<NumForceUpdate && ncount<nthis; 
	  i=P[i].ForceFlag, timelinecounter++)
	{
	  if((P[i].Type&7)==0)
	    {
	      if(P[i].Type&8)
		{
		  place= noffset[ThisTask] + nexport; 
                  nexport++;
		  P[i].Type&= 7;
		}
	      else
		place= nbuffer + (ncount-nexport);

	      /* add to the acceleration, do not reset it. (gravity is already computed) */
	      for(k=0;k<3;k++)   
		P[i].Accel[k] += HydroDataResult[place].Acc[k];
	      
	      SphP[i].DtEgySpec= HydroDataResult[place].DtEgySpec;

	      ncount++;
	    }
	}
      
      nstart=i;                  /* continue at this point in the timeline in next iteration */
      startcounter=timelinecounter;

      npleft-= ncount;
      ntotleft-= ntot;
    }

  free(nrecv);
  free(noffset);

 
#ifdef DEB
  if(ThisTask==0)
    { 	
      FdDEB=fopen("hydra.deb","w");
      fprintf(FdDEB," %d \n",All.NumCurrentTiStep);
      fclose(FdDEB);
    }
#endif
}








