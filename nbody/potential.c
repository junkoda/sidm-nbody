#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>


#include "allvars.h"
#include "proto.h"



/* This function computes the gravitational potential for ALL the particles.
 * It expects that the particles are predicted to the current time.
 * The routine constructs a new force-tree. 
 */
void compute_potential(void)  
{
  int     i,k;
  int     level,sendTask,recvTask;
  int     np,ntot,ntotleft,npleft,nthis;     
  int     nstart,nbuffer,ncount,nchunk;
  int     *nrecv,*noffset;
  double  fac;
  double  t0,t1;
  MPI_Status status;
  double r2;

  t0=second();

  if(All.ComovingIntegrationOn)
    set_softenings();

  if(ThisTask==0)
    {
      printf("Start computation of potential for all particles...\n"); 
      fflush(stdout);  
    }


  np=NumPart;
  
  MPI_Allreduce(&np, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(NoCostFlag==0)
    force_costevaluate();

  force_treebuild();
  All.NumForcesSinceLastTreeConstruction= All.TreeUpdateFrequency*All.TotNumPart ; /* ensures that new tree will be constructed next time*/
  NoCostFlag=1;


 
  ntotleft= ntot;  /* particles left for all tasks together */
  npleft= np;      /* particles left for this task */
  nstart= 1;       /* particle index of start of current bunch for this task */

  nthis = ((double)np*All.BunchSizeForce)/ntot; /* maximum size of chunck coming from this task */

  nchunk= np/nthis;                             /* number of chunks needed */  
  if((np%nthis)>0) nchunk+=1;

  nrecv=malloc(sizeof(int)*NTask);  /* particle numbers of each task constituating current bunch */
  noffset=malloc(sizeof(int)*NTask);

  while(ntotleft>0)
    {
      if(nthis>npleft)
	nthis=npleft;

      for(i=nstart,ncount=0; i<=NumPart && ncount<nthis; i++)
	ncount++;
      
      MPI_Allgather(&ncount, 1, MPI_INT, nrecv, 1, MPI_INT, MPI_COMM_WORLD);

      for(i=nbuffer=0; i<NTask; i++)
	nbuffer += nrecv[i];

      for(i=1,noffset[0]=0; i<NTask; i++)
	noffset[i]= noffset[i-1]+nrecv[i-1];
	 
      for(i=nstart, ncount=0; i<=NumPart && ncount<nthis; i++, ncount++)
	{
	  for(k=0;k<3;k++)
	    GravDataIn[noffset[ThisTask] + ncount].Pos[k]= P[i].PosPred[k];
	    
	  GravDataIn[noffset[ThisTask] + ncount].Type= P[i].Type;
	  GravDataIn[noffset[ThisTask] + ncount].OldAcc= P[i].OldAcc;
	}


      /* now start big communication */

      for(level=1; level<NTask; level++)
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ level; 

	  MPI_Sendrecv(&GravDataIn[noffset[sendTask]], nrecv[sendTask]*sizeof(struct gravdata_in), MPI_BYTE, recvTask, TAG_ANY,
		       &GravDataIn[noffset[recvTask]], nrecv[recvTask]*sizeof(struct gravdata_in), MPI_BYTE, recvTask, TAG_ANY ,
		       MPI_COMM_WORLD,&status);
	}




      /* ok, all preparations are done, we can compute the partial potentials
       * for the particles in the buffer
       */
      for(i=0;i<nbuffer;i++)
	force_treeevaluate_potential(i);    
    
      
      /* communicate the results and sum them up 
       */
      for(level=1; level<NTask; level++)
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ level; 
	  
	  MPI_Sendrecv(&GravDataPotential[noffset[recvTask]], nrecv[recvTask], MPI_DOUBLE, recvTask, TAG_ANY,
		       &GravDataPartialPotential[0], nrecv[sendTask], MPI_DOUBLE, recvTask, TAG_ANY,
		       MPI_COMM_WORLD, &status);
	
	  for(i=0; i<nrecv[ThisTask]; i++)
	    GravDataPotential[noffset[ThisTask] + i] += GravDataPartialPotential[i];
	}


      /* transfer the result to the particles 
       */
      for(i=nstart,ncount=0; i<=NumPart && ncount<nrecv[ThisTask]; i++)
	{
	  P[i].Potential =  GravDataPotential[noffset[ThisTask]+ncount];
	  P[i].Potential += P[i].Mass/All.SofteningTable[P[i].Type];  /* removes self energy */
	  
	  ncount++;
	}
      nstart= i;
      npleft -= nrecv[ThisTask];

      ntotleft -= nbuffer;
    }

  free(nrecv);
  free(noffset);




  if(All.ComovingIntegrationOn)
    {
      fac=0.5*All.Omega0*All.Hubble*All.Hubble;
 
      for(i=1; i<=NumPart; i++)
	{
#ifdef PERIODIC
	  P[i].Potential = All.G*P[i].Potential;
#else
	  for(k=0, r2=0; k<3; k++)
	    r2 += P[i].PosPred[k]*P[i].PosPred[k];

	  P[i].Potential = All.G*P[i].Potential - fac*r2;
#endif
	}
    }
  else
    {
      fac= -0.5*All.OmegaLambda*All.Hubble*All.Hubble;

      for(i=1;i<=NumPart;i++)
	{
	  P[i].Potential *= All.G;

	  if(fac!=0)
	     {
	       for(k=0,r2=0;k<3;k++)
		 r2 += P[i].PosPred[k]*P[i].PosPred[k];
	       
	       P[i].Potential += fac*r2;
	     }
	}
    }

  
  if(ThisTask==0)
    {
      printf("potential done.\n"); fflush(stdout);
    }

  t1=second();

  All.CPU_Potential+= timediff(t0,t1);
}

























