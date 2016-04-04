#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"




/* This function computes the local density for each active SPH particle,
 * the number of neighbours in the current smoothing radius,
 * and the divergence and rotation of the velocity field.
 * The pressure is updated as well.
 * If a particle with its smoothing region is fully inside the local
 * domain, it is not exported to the other processors.
 */
void density(void)
{
  double  rotv[3];
  double  h=0, hinv=0, hinv3=0, hinv4=0;
  double  rho, divv, wk, dwk;
  double  dx, dy, dz, r, r2, u, mass_j;
  double  dt;
  int     i,j,k,ii,n, numngb;
  double  hubble_a, prefac=0;
  float   *r2list;
  int     *ngblist;
  double  tstart,tend;
  int     ntot,ntotleft,npleft,nthis;     
  int     nstart,nbuffer,ncount,nchunk;
  int     timelinecounter,startcounter;   
  int     *nrecv,*noffset;
  int     level,sendTask,recvTask;
  int     place, nexport;
  MPI_Status status;
  
  
  if(All.ComovingIntegrationOn)
    {
      /* Factors for comoving integration of hydro */
      hubble_a=All.Hubble*sqrt(All.Omega0/(All.Time*All.Time*All.Time) 
			       + (1-All.Omega0-All.OmegaLambda)/(All.Time*All.Time) + All.OmegaLambda);
      
      prefac= 1.0/(hubble_a*pow(All.Time, 1.5));
    }



  /* `NumSphUpdate' is the number of particles on this processor that want a density update */
  MPI_Allreduce(&NumSphUpdate, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  /* note: it can happen that this processor has NumForceUpdate==0, but that ntot>0 */
  
  /* compute maximum size of bunch that may come from this task */
  if(NumSphUpdate==0)
    nthis=1;
  else
    nthis = ((double)NumSphUpdate*(All.BunchSizeDensity-NTask))/ntot + 1; 


  if(NumSphUpdate>0) 
    {
      nchunk= NumSphUpdate/nthis;                             /* number of bunches needed */  
      if((NumSphUpdate%nthis)>0) nchunk+=1;
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
		if(P[i].PosPred[j] < (DomainMin[0][j]+SphP[i].Hsml))
		  break;
		if(P[i].PosPred[j] > (DomainMax[0][j]-SphP[i].Hsml))
		  break;
	      }

	    if(j!=3)  /* particle lies NOT completely inside . needs to be sent to other processors */
	      {
		nexport++;
                P[i].Type |= 8;
	      }
	  }

      MPI_Allgather(&nexport, 1, MPI_INT, nrecv, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Allreduce(&ncount, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      for(i=nbuffer=0; i<NTask; i++)  /* compute length of common list */
	nbuffer += nrecv[i];

      for(i=1, noffset[0]=0; i<NTask; i++) /* set-up offset-table */
	noffset[i]= noffset[i-1] + nrecv[i-1];

      /* fill in the own particles at the right place */
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

	      for(k=0; k<3; k++)
		{
		  DensDataIn[place].Pos[k]= P[i].PosPred[k];
		  DensDataIn[place].Vel[k]= P[i].VelPred[k];
		}
	      DensDataIn[place].Hsml = SphP[i].Hsml;
	      
	      ncount++;
	    }
	}
	 

      /* now start big communication */

      tstart=second();      
      for(level=1;level<NTask;level++)
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ level; 

	  MPI_Sendrecv(&DensDataIn[noffset[sendTask]], nrecv[sendTask]*sizeof(struct densdata_in), MPI_BYTE, recvTask, TAG_ANY, 
	               &DensDataIn[noffset[recvTask]], nrecv[recvTask]*sizeof(struct densdata_in), MPI_BYTE, recvTask, TAG_ANY, MPI_COMM_WORLD, &status);
	}
      tend=second();
      All.CPU_CommSum+= timediff(tstart, tend);



      /* ok, all preparations are done. Now density evaluation and velocity stuff */
      
      for(i=0; i<(nbuffer+ncount-nexport); i++)
	{
	  rho = divv = rotv[0] = rotv[1] = rotv[2] = 0;

	  numngb= ngb_treefind_variable(&DensDataIn[i].Pos[0], DensDataIn[i].Hsml, 0, &ngblist, &r2list);   

	  if(numngb>0)
	    {
	      h = DensDataIn[i].Hsml;
	      hinv  = 1.0/h;
	      hinv3 = hinv*hinv*hinv;
	      hinv4 = hinv3*hinv;
	    }

	  for(n=0; n<numngb; n++)
	    {
	      j  = ngblist[n]+1; 

	      dx = DensDataIn[i].Pos[0] - P[j].PosPred[0];
	      dy = DensDataIn[i].Pos[1] - P[j].PosPred[1];
	      dz = DensDataIn[i].Pos[2] - P[j].PosPred[2];
#ifdef PERIODIC
	      dx= periodic(dx);
	      dy= periodic(dy);
	      dz= periodic(dz);
#endif
	      r2 = dx*dx + dy*dy + dz*dz;
	      
	      r = sqrt(r2);
	      
	      if(r<h)
		{
		  u = r*hinv;
		  
		  ii = (int)(u*KERNEL_TABLE);
		  
		  wk =hinv3*( Kernel[ii]  + (Kernel[ii+1]-Kernel[ii])*(u-KernelRad[ii])*KERNEL_TABLE);
		  dwk=hinv4*( KernelDer[ii] + (KernelDer[ii+1]-KernelDer[ii])*(u-KernelRad[ii])*KERNEL_TABLE);
		  
#ifdef SFR
#ifdef CLOUDS
		  mass_j = P[j].Mass - SphP[j].FormedStellarMass- SphP[j].CloudMass;
#else
		  mass_j = P[j].Mass - SphP[j].FormedStellarMass;
#endif
#else
		  mass_j = P[j].Mass;
#endif
	
		  rho += mass_j * wk;
		  
		  if(r>0)
		    {
		      divv -= mass_j * dwk/r *
			( dx * (DensDataIn[i].Vel[0] - P[j].VelPred[0])
		        + dy * (DensDataIn[i].Vel[1] - P[j].VelPred[1])
		        + dz * (DensDataIn[i].Vel[2] - P[j].VelPred[2]) );
		  
		      rotv[0] += mass_j * dwk/r *
			  (  dz * (DensDataIn[i].Vel[1] - P[j].VelPred[1])
			   - dy * (DensDataIn[i].Vel[2] - P[j].VelPred[2]) );

		      rotv[1] += mass_j * dwk/r *
		           (  dx * (DensDataIn[i].Vel[2] - P[j].VelPred[2])
      	                    - dz * (DensDataIn[i].Vel[0] - P[j].VelPred[0]) );

		      rotv[2] += mass_j * dwk/r *
		           ( dy * (DensDataIn[i].Vel[0] - P[j].VelPred[0])
			   - dx * (DensDataIn[i].Vel[1] - P[j].VelPred[1]) );
		    }
		}
	      
	    }
	  
	  DensDataResult[i].Rho= rho;
	  DensDataResult[i].Div= divv;
	  DensDataResult[i].Ngb= numngb;
	  DensDataResult[i].Rot[0]= rotv[0];  
	  DensDataResult[i].Rot[1]= rotv[1];  
	  DensDataResult[i].Rot[2]= rotv[2];  
	}  



      /* now communicate contributions, and sum up */

      tstart=second();
      for(level=1; level<NTask; level++)
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ level; 
	
          MPI_Sendrecv(&DensDataResult[noffset[recvTask]], nrecv[recvTask]*sizeof(struct densdata_out), MPI_BYTE, recvTask, TAG_ANY,
	               &DensDataPartialResult[0], nrecv[sendTask]*sizeof(struct densdata_out), MPI_BYTE, recvTask, TAG_ANY, MPI_COMM_WORLD, &status);

	  for(i=0; i<nrecv[ThisTask]; i++)
	    {
	      DensDataResult[noffset[ThisTask] + i].Rho += DensDataPartialResult[i].Rho;
	      DensDataResult[noffset[ThisTask] + i].Ngb += DensDataPartialResult[i].Ngb;
	      DensDataResult[noffset[ThisTask] + i].Div += DensDataPartialResult[i].Div;
	      
	      for(k=0;k<3;k++)
		DensDataResult[noffset[ThisTask] + i].Rot[k] += DensDataPartialResult[i].Rot[k];
	    }
	}
      tend=second();
      All.CPU_CommSum+= timediff(tstart,tend);



      /* transfer the result to the particles */

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
	      
	      rotv[0]= DensDataResult[place].Rot[0];
	      rotv[1]= DensDataResult[place].Rot[1];
	      rotv[2]= DensDataResult[place].Rot[2];
	      divv=    DensDataResult[place].Div;
	      rho=     DensDataResult[place].Rho;

	      SphP[i].NumNgb= DensDataResult[place].Ngb;

	      SphP[i].CurlVel= sqrt(rotv[0]*rotv[0] + rotv[1]*rotv[1] + rotv[2]*rotv[2])/rho;
	      SphP[i].Density= SphP[i].DensityPred = rho;
	      SphP[i].DivVel= divv/rho;

	      if(All.ComovingIntegrationOn) /* comoving variables */
		{
		  SphP[i].DtDensity= - prefac*divv;
		  SphP[i].DtHsml= - SphP[i].Hsml*SphP[i].DtDensity/(3*SphP[i].Density);
		}
	      else
		{  /* normal space  */
		  SphP[i].DtDensity= -divv;
		  SphP[i].DtHsml= SphP[i].Hsml*SphP[i].DivVel/3;
		}
    
	      /* make sure that predicted values of density and smoothing length
		 cannot take on negative values */

	      dt = 2*(All.Time - P[i].CurrentTime);  /* this is the actual time-step */
	      if(dt>0)
		{
		  SphP[i].DtHsml+= SphP[i].Hsml/(2*dt)*(pow(((double)All.DesNumNgb)/SphP[i].NumNgb, 1.0/3)-1);
		  
		  SphP[i].DtDensity= dmax(-0.9*SphP[i].Density/dt , SphP[i].DtDensity);
		  SphP[i].DtHsml   = dmax(-0.9*SphP[i].Hsml/dt ,    SphP[i].DtHsml);
		}

	      SphP[i].Pressure = GAMMA_MINUS1*(SphP[i].EgySpecPred)*SphP[i].DensityPred;

	      ncount++;
	    }
	}



      nstart=i;      /* continue at this point in the timeline in next iteration */
      startcounter=timelinecounter;

      npleft-= ncount;
      ntotleft-= ntot;
    }

  free(nrecv);
  free(noffset);


#ifdef DEB
  if(ThisTask==0)
    { 	
      FdDEB=fopen("density.deb","w");
      fprintf(FdDEB," %d \n",All.NumCurrentTiStep);
      fclose(FdDEB);
    }
#endif
}



/*  this function wraps the distance x to the closest image
 *  for the given box size
 */
#ifndef INLINE
double INLINE_FUNC periodic(double x)
{
  while(x > All.BoxHalf)
    x -=All.BoxSize;

  while(x < -All.BoxHalf)
      x+=All.BoxSize;

  return x;
}
#endif




/* the function below detects particles that have a number of neighbours 
 * outside the allowed tolerance range. For these, particles the smoothing
 * length is adjusted accordingly. Note that the smoothing length is
 * not allowed to fall below the bound set by MinGasHsml
 */
void ensure_neighbours(int mode)
{
#define MAXITER 30

  int    i, ntot, last=0;
  float  *r2list;
  int    *ngblist, count, candidates;
  int    iter=0;
  double save;

#ifdef VELDISP
#define PPP P
#else
#ifdef SIDM
#define PPP P
#else
#define PPP SphP
#endif
#endif
 
  for(i=IndFirstUpdate, count=0, candidates=0; count<NumForceUpdate; i=P[i].ForceFlag, count++)   
    {  
      if(P[i].Type==0)
	{
	  if( SphP[i].NumNgb < (All.DesNumNgb-All.MaxNumNgbDeviation) || 
	     (SphP[i].NumNgb > (All.DesNumNgb+All.MaxNumNgbDeviation) && SphP[i].Hsml>(1.01*All.MinGasHsml)))
	    candidates++;
	}
    }

  /* Note: When the comoving softening is made time-dependent, the value of All.MinGasHsml will 
     drop slowly with time. So if we would check with SphP[i].Hsml==All.MinGasHsml, all the partciles
     that are bounded by the lower cut-off in the SPH resolution would have to be readjusted EVERY timestep,
     even so their smoothing length is already accurate. Now we only adjust Hsml to lower values if we are at least
     1% above the minimum. In this way, a varying minimimum resolution cut-off is followed efficiently.
  */

  MPI_Reduce(&candidates, &ntot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ntot, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  if(ntot>0)
    {
      if(ThisTask==0)
	{
	  printf("\n%d particles have too few/too many neighbours!\n", ntot);
	  printf("Now fixing that...\n"); 
	}
      
      /* we will here make use of  P[i].Accel[j] and store lower
       * and upper boundary of intervals for smoothing range bisection
       * Note that P[i].Accel[j] will be recomputed later on (force computation).
       */
      
      for(i=IndFirstUpdate, count=0; count<NumForceUpdate; i=P[i].ForceFlag, count++)   
	if(P[i].Type==0)
	  PPP[i].Left= PPP[i].Right= 0;
    
      do
	{
	  for(i=1, NumForceUpdate= NumSphUpdate= 0; i<=N_gas; i++) 
	    { 
	      if( SphP[i].NumNgb < (All.DesNumNgb-All.MaxNumNgbDeviation) || 
		 (SphP[i].NumNgb > (All.DesNumNgb+All.MaxNumNgbDeviation) && SphP[i].Hsml>(1.01*All.MinGasHsml)))
		{
		  if(PPP[i].Left>0 && PPP[i].Right>0)
		    if((PPP[i].Right-PPP[i].Left) < 1.0e-3 * PPP[i].Left)
		      continue;

		  if(NumForceUpdate==0)
		    IndFirstUpdate= i;
		  else
		    P[last].ForceFlag= i;
		  
		  NumForceUpdate++;
		  NumSphUpdate++;
		  last=i;
		  
		  if(SphP[i].NumNgb < (All.DesNumNgb-All.MaxNumNgbDeviation))
		    PPP[i].Left= dmax(SphP[i].Hsml, PPP[i].Left);
		  else
		    if(PPP[i].Right!=0)
		      {
			if(SphP[i].Hsml<PPP[i].Right)
			  PPP[i].Right= SphP[i].Hsml;
		      }
		    else
		      PPP[i].Right= SphP[i].Hsml;
		}
	    }
	  
	  MPI_Allreduce(&NumSphUpdate, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	  
      
	  if(ntot>0)
	    {
	      if(ThisTask==0)
		printf("ngb iteration %d.  still %d particles\n", iter, ntot);

	      for(i=IndFirstUpdate,count=0; count<NumForceUpdate; i=P[i].ForceFlag, count++)   
		{
		  if(iter >= 20)
		    {
		      printf("i=%d ID=%d Hsml=%g Left=%g Right=%g Ngbs=%d Right-Left=%g\n   pos=(%g|%g|%g)\n",
			     i, P[i].ID, SphP[i].Hsml, PPP[i].Left, PPP[i].Right, SphP[i].NumNgb, PPP[i].Right-PPP[i].Left,
			     P[i].PosPred[0], P[i].PosPred[1], P[i].PosPred[2]);
		    }
		  
		  if(iter == MAXITER)
		    {
		      printf("ThisTask=%d Mi=(%g|%g|%g) Ma=(%g|%g|%g)\n", ThisTask,
			     DomainMin[0][0], DomainMin[0][1], DomainMin[0][2],
			     DomainMax[0][0], DomainMax[0][1], DomainMax[0][2]);
		      printf("i=%d ID=%d coord=(%g|%g|%g)\n", i, P[i].ID, P[i].PosPred[0], P[i].PosPred[1], P[i].PosPred[2]);
		      printf("ngb_treefind= %g\n", sqrt(ngb_treefind(P[i].PosPred , All.DesNumNgb, 0, 0, &ngblist, &r2list)));
		    }

		  if(PPP[i].Left==0 || PPP[i].Right==0) 
		    {
		      if(PPP[i].Right==0 && SphP[i].NumNgb<15)
			{
			  SphP[i].Hsml= sqrt(ngb_treefind(P[i].PosPred , All.DesNumNgb, 0, 0, &ngblist, &r2list));  
			}
		      else
			{
			  SphP[i].Hsml=  SphP[i].Hsml*( 0.5 + 0.5*pow(SphP[i].NumNgb/((double)All.DesNumNgb), -1.0/3));
			}
		    }
		  else
		    {
		      SphP[i].Hsml=  0.5*(PPP[i].Left + PPP[i].Right);
		    }

		  if(SphP[i].Hsml < All.MinGasHsml)
		    SphP[i].Hsml= All.MinGasHsml;
		}
	      
	      density();

	      iter++;

	      if(iter > MAXITER)
		{
		  fprintf(stdout, "failed to converge in function ensure_neighbours()\n");
		  endrun(1155);
		}
	    }
	}
      while(ntot>0);
      

      if(mode==0)   /* restore timeline to active particles */
	{
	  save= All.TimeStep;
	  find_next_time();
	  All.TimeStep= save;
	}
      else  /* make all particles active again */
	{
	  for(i=1; i<=NumPart; i++) 
	    P[i].ForceFlag=i+1;
	  
	  P[NumPart].ForceFlag=1; IndFirstUpdate=1; NumForceUpdate=NumPart; NumSphUpdate=N_gas;
	}
    }
#undef MAXITER
}


/* This function finds a region of space that is guaranteed to contain 
 * only predicted gas coordinates from the local domain. With this
 * information, the communication in the hydro-dynamical computations
 * can be reduces. 
 *
 * Each processor first collects the extensions of all the gas-domains 
 * on all the other processors. This extension is stored in DomainMax[],
 * and DomainMin[], and is updated in the function ngb_update_nodes(),
 * or the tree construction.
 * Then the local region is cut with each of the other regions.
 * If there is a non-vanishing overlap, the local interior region is reduced 
 * by this overlap while trying to keep its volume as large as possible.
 * (we can't have in the largest volume possible, since we 
 * continue to describe the interior region with a simple rectangular domain.) 
 */
void determine_interior(void)
{
  int   i, j, j0, j1, j2, flag, type;
  float common_min[3], common_max[3], shift[3];
  float innerMax[3], innerMin[3];
  float V, Vmax;
#ifdef PERIODIC
  int ix, iy, iz;
#endif


  for(type=0; type<5; type++)
    {

#ifndef VELDISP
#ifndef SIDM
      if(type>0)
	break;
#endif
#endif

      if(Ntype[type]>0)
	{
	  MPI_Allgather(&DomainMin[type][0], 3, MPI_FLOAT, InteriorMin, 3, MPI_FLOAT, MPI_COMM_WORLD);
	  MPI_Allgather(&DomainMax[type][0], 3, MPI_FLOAT, InteriorMax, 3, MPI_FLOAT, MPI_COMM_WORLD);

	  for(i=0; i<NTask; i++)
	    {
	      if(i!=ThisTask)
		{
#ifdef PERIODIC
		  for(ix=-1, shift[0]=-All.BoxSize; ix<=1; ix++, shift[0]+=All.BoxSize)
		    for(iy=-1, shift[1]=-All.BoxSize; iy<=1; iy++, shift[1]+=All.BoxSize)
		      for(iz=-1, shift[2]=-All.BoxSize; iz<=1; iz++, shift[2]+=All.BoxSize)
#else
			shift[0]= shift[1]= shift[2]= 0;
#endif
		  {
		    for(j=0, flag=1; j<3; j++)
		      {
			/* find overlapping region */
			common_min[j]= dmax(DomainMin[type][j], InteriorMin[i*3+j] + shift[j]);
			common_max[j]= dmin(DomainMax[type][j], InteriorMax[i*3+j] + shift[j]);
			
			if(common_max[j] <= common_min[j]) /* no overlap in this dimension, so we are done here */
			  {
			    flag=0; 
			    break;
			  }
		      }
		    
		    if(flag) /* we need to reduce the local domain */
		      {
			/* find the reduction that leaves the largest volume */
			
			Vmax=-1;
			memcpy(innerMin, DomainMin[type], sizeof(float)*3);
			memcpy(innerMax, DomainMax[type], sizeof(float)*3);
			
			for(j0=0, j1=1, j2=2; j0<3; j0++, j1++, j2++)
			  { 
			    if(j1>2) 
			      j1=0;
			    
			    if(j2>2) 
			      j2=0;
			    
			    V= (common_min[j0]   - DomainMin[type][j0])*
			      (DomainMax[type][j1] - DomainMin[type][j1])*
			      (DomainMax[type][j2] - DomainMin[type][j2]);
			    
			    if(V>Vmax)
			      {
				Vmax= V;
				memcpy(innerMin, DomainMin[type], sizeof(float)*3);
				memcpy(innerMax, DomainMax[type], sizeof(float)*3);
				
				innerMax[j0]= common_min[j0];
			      }
			    
			    V= (DomainMax[type][j0] - common_max[j0])*
			      (DomainMax[type][j1] - DomainMin[type][j1])*
			      (DomainMax[type][j2] - DomainMin[type][j2]);
			    
			    if(V>Vmax)
			      {
				Vmax= V;
				memcpy(innerMin, DomainMin[type], sizeof(float)*3);
				memcpy(innerMax, DomainMax[type], sizeof(float)*3);
				
				innerMin[j0]= common_max[j0];
			      }
			  }
			
			memcpy(DomainMin[type], innerMin,  sizeof(float)*3);
			memcpy(DomainMax[type], innerMax,  sizeof(float)*3);
		      }
		  }
		}
	    }
	}
    }
}















