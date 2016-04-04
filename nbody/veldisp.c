#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"


#ifdef VELDISP

/* This function computes the local velocity dispersion
 * of each particle (among particles of the same type)
 */
void veldisp(void)
{
  double  vsum[3], v2sum[3];
  int     i,ii,j,k,n, numngb;
  double  h, hinv, hinv3, r, u, wk, rho;	  
  float   *r2list;
  int     *ngblist;
  double  tstart,tend;
  int     ntot,ntotleft,npleft,nthis;     
  int     nstart,nbuffer,ncount,nchunk;
  int     timelinecounter,startcounter;   
  int     *nrecv,*noffset;
  int     level,sendTask,recvTask;
  int     place, nexport;
  int     num_collisionless;
  MPI_Status status;
  
  num_collisionless= NumForceUpdate-NumSphUpdate;

  MPI_Allreduce(&num_collisionless, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  /* compute maximum size of bunch that may come from this task */
  if(num_collisionless==0)
    nthis=1;
  else
    nthis = ((double)(num_collisionless)*(All.BunchSizeVelDisp-NTask))/ntot + 1; 


  if(num_collisionless>0) 
    {
      nchunk= num_collisionless/nthis;                             /* number of bunches needed */  
      if((num_collisionless%nthis)>0) nchunk+=1;
    }
  else
    nchunk=0;


  nrecv=   malloc(sizeof(int)*NTask);  /* list of particle numbers that constituate current bunch */
  noffset= malloc(sizeof(int)*NTask);  /* offsets of bunches in common list */

  ntotleft= ntot;              /* particles left for all tasks together */
  npleft= num_collisionless;   /* particles left for this task */

  nstart= IndFirstUpdate;      /* first particle for this task */
  startcounter=0;

  while(ntotleft>0)
    {
      if(nthis>npleft)
	nthis=npleft;

      for(i=nstart, ncount=0, nexport=0, timelinecounter=startcounter; 
	  timelinecounter<NumForceUpdate && ncount<nthis; 
	  i=P[i].ForceFlag, timelinecounter++)
	if(P[i].Type>0)
	  {
	    ncount++;
	    
	    for(j=0; j<3; j++)
	      {
		if(P[i].PosPred[j] < (DomainMin[P[i].Type][j] + P[i].HsmlVelDisp))
		  break;
		if(P[i].PosPred[j] > (DomainMax[P[i].Type][j] - P[i].HsmlVelDisp))
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
	if(P[i].Type>0)
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
		VelDispDataIn[place].Pos[k]= P[i].PosPred[k];
	      }
	    VelDispDataIn[place].Hsml  = P[i].HsmlVelDisp;
	    VelDispDataIn[place].Type =  P[i].Type&7;
	    
	    ncount++;
	  }

      /* now start big communication */

      tstart=second();      
      for(level=1;level<NTask;level++)
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ level; 
	  
	  MPI_Sendrecv(&VelDispDataIn[noffset[sendTask]], nrecv[sendTask]*sizeof(struct veldispdata_in), MPI_BYTE, recvTask, TAG_ANY, 
	               &VelDispDataIn[noffset[recvTask]], nrecv[recvTask]*sizeof(struct veldispdata_in), MPI_BYTE, recvTask, TAG_ANY, MPI_COMM_WORLD, &status);
	}
      tend=second();
      All.CPU_CommSum+= timediff(tstart, tend);

      /* ok, all preparations are done. Now density evaluation and velocity stuff */

      for(i=0; i<(nbuffer+ncount-nexport); i++)
	{
	  rho= vsum[0] = vsum[1] = vsum[2] = v2sum[0] = v2sum[1] = v2sum[2] = 0;

	  numngb= ngb_treefind_variable(&VelDispDataIn[i].Pos[0], 
					VelDispDataIn[i].Hsml, 
					VelDispDataIn[i].Type, 
					&ngblist, &r2list);   

	  h = VelDispDataIn[i].Hsml;
	  hinv  = 1.0/h;
	  hinv3 = hinv*hinv*hinv;
	  
	  for(n=0; n<numngb; n++)
	    {
	      j  = ngblist[n]+1; 
	      r  = r2list[n];
	      
	      if(r<h)
		{
		  u = r*hinv;
		  ii = (int)(u*KERNEL_TABLE);
		  wk =hinv3*( Kernel[ii]  + (Kernel[ii+1]-Kernel[ii])*(u-KernelRad[ii])*KERNEL_TABLE);
		  rho += P[j].Mass * wk;
		}	

	      for(k=0; k<3; k++)
		{
		  vsum[k] +=  P[j].VelPred[k];
		  v2sum[k] += P[j].VelPred[k]*P[j].VelPred[k];
		}
	    }
	  
	  VelDispDataResult[i].Ngb= numngb;
	  VelDispDataResult[i].Rho= rho;

	  for(k=0; k<3; k++)
	    {
	      VelDispDataResult[i].Vsum[k]= vsum[k];
	      VelDispDataResult[i].V2sum[k]= v2sum[k];
	    }
	}  

      /* now communicate contributions, and sum up */

      tstart=second();
      for(level=1; level<NTask; level++)
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ level; 
	
          MPI_Sendrecv(&VelDispDataResult[noffset[recvTask]], nrecv[recvTask]*sizeof(struct veldispdata_out), MPI_BYTE, recvTask, TAG_ANY,
	               &VelDispDataPartialResult[0], nrecv[sendTask]*sizeof(struct veldispdata_out), MPI_BYTE, recvTask, TAG_ANY, MPI_COMM_WORLD, &status);

	  for(i=0; i<nrecv[ThisTask]; i++)
	    {
	      VelDispDataResult[noffset[ThisTask] + i].Ngb += VelDispDataPartialResult[i].Ngb;
	      VelDispDataResult[noffset[ThisTask] + i].Rho += VelDispDataPartialResult[i].Rho;
	      
	      for(k=0;k<3;k++)
		{
		  VelDispDataResult[noffset[ThisTask] + i].Vsum[k] += VelDispDataPartialResult[i].Vsum[k];
		  VelDispDataResult[noffset[ThisTask] + i].V2sum[k] += VelDispDataPartialResult[i].V2sum[k];
		}
	    }
	}
      tend=second();
      All.CPU_CommSum+= timediff(tstart,tend);

      /* transfer the result to the particles */
      
      for(i=nstart, ncount=0, nexport=0, timelinecounter=startcounter; 
	  timelinecounter<NumForceUpdate && ncount<nthis; 
	  i=P[i].ForceFlag, timelinecounter++)
	if(P[i].Type>0)
	  
	  {
	    if(P[i].Type&8)
	      {
		place= noffset[ThisTask] + nexport; 
		nexport++;
		P[i].Type&= 7;
	      }
	    else
	      place= nbuffer + (ncount-nexport);
	    
	    P[i].NgbVelDisp = VelDispDataResult[place].Ngb;
	    P[i].DensVelDisp = VelDispDataResult[place].Rho;
	    
	    P[i].VelDisp = 0;
	    
	    if(P[i].NgbVelDisp>0)
	      {
		for(k=0; k<3; k++)
		  {
		    VelDispDataResult[place].V2sum[k]/= P[i].NgbVelDisp;
		    VelDispDataResult[place].Vsum[k]/= P[i].NgbVelDisp;
		    
		    P[i].VelDisp += VelDispDataResult[place].V2sum[k] - 
		      VelDispDataResult[place].Vsum[k]*VelDispDataResult[place].Vsum[k];
		  }
	      }
	    
	    if(P[i].VelDisp > 0)
	      P[i].VelDisp= sqrt(P[i].VelDisp);
	    
	    ncount++;
	  }
      
      nstart=i;      /* continue at this point in the timeline in next iteration */
      startcounter=timelinecounter;
      
      npleft-= ncount;
      ntotleft-= ntot;
    }

  free(nrecv);
  free(noffset);
}






/* the function below detects particles that have a number of neighbours 
 * outside the allowed tolerance range. For these, particles the smoothing
 * length is adjusted accordingly. Note that the smoothing length is
 */
void veldisp_ensure_neighbours(int mode)
{
#define MAXITER 30

  int    i, ntot, last=0;
  float  *r2list;
  int    *ngblist, count, candidates;
  int    iter=0;
  double save;
  
  for(i=IndFirstUpdate, count=0, candidates=0; count<NumForceUpdate; i=P[i].ForceFlag, count++)   
    {
      if(P[i].Type>0)
	{
	  if( P[i].NgbVelDisp < (All.DesNumNgb-All.MaxNumNgbDeviation) || 
	      (P[i].NgbVelDisp > (All.DesNumNgb+All.MaxNumNgbDeviation)))
	    candidates++;
	}
    }

  MPI_Reduce(&candidates, &ntot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ntot, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  if(ntot>0)
    {
      if(ThisTask==0)
	{
	  printf("\n%d particles have too few/too many neighbours in veldisp calculation!\n", ntot);
	  printf("Now fixing that...\n"); 
	}
      
      for(i=IndFirstUpdate, count=0; count<NumForceUpdate; i=P[i].ForceFlag, count++)   
	P[i].Left= P[i].Right= 0;
      
      do
	{
	  for(i=1+N_gas, NumForceUpdate= 0, NumSphUpdate=0; i<=NumPart; i++) 
	    { 
	      if( P[i].NgbVelDisp < (All.DesNumNgb-All.MaxNumNgbDeviation) || 
		  P[i].NgbVelDisp > (All.DesNumNgb+All.MaxNumNgbDeviation))
		{
		  if(P[i].Left>0 && P[i].Right>0)
		    if((P[i].Right-P[i].Left) < 1.0e-3 * P[i].Left)
		      continue;

		  if(NumForceUpdate==0)
		    IndFirstUpdate= i;
		  else
		    P[last].ForceFlag= i;

		  NumForceUpdate++;
		  last=i;
		  
		  if(P[i].NgbVelDisp < (All.DesNumNgb-All.MaxNumNgbDeviation))
		    P[i].Left= dmax(P[i].HsmlVelDisp, P[i].Left);
		  else
		    if(P[i].Right!=0)
		      {
			if(P[i].HsmlVelDisp<P[i].Right)
			  P[i].Right= P[i].HsmlVelDisp;
		      }
		    else
		      P[i].Right= P[i].HsmlVelDisp;
		}
	    }
	  
	  MPI_Allreduce(&NumForceUpdate, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	  
	  if(ntot>0)
	    {
	      if(ThisTask==0)
		printf("ngb iteration %d.  still %d particles\n", iter, ntot);

	      for(i=IndFirstUpdate, count=0; count<NumForceUpdate; i=P[i].ForceFlag, count++)   
		{
		  if(iter >= 20)
		    {
		      printf("i=%d ID=%d Hsml=%g Left=%g Right=%g Ngbs=%d Right-Left=%g\n   pos=(%g|%g|%g)\n",
			     i, P[i].ID, P[i].HsmlVelDisp, P[i].Left, P[i].Right, P[i].NgbVelDisp, P[i].Right-P[i].Left,
			     P[i].PosPred[0], P[i].PosPred[1], P[i].PosPred[2]);
		    }
		  
		  if(iter == MAXITER)
		    {
		      printf("ThisTask=%d Mi=(%g|%g|%g) Ma=(%g|%g|%g)\n", ThisTask,
			     DomainMin[P[i].Type][0], DomainMin[P[i].Type][1], DomainMin[P[i].Type][2],
			     DomainMax[P[i].Type][0], DomainMax[P[i].Type][1], DomainMax[P[i].Type][2]);
		      printf("i=%d ID=%d coord=(%g|%g|%g)\n", i, P[i].ID, P[i].PosPred[0], P[i].PosPred[1], P[i].PosPred[2]);
		      printf("ngb_treefind= %g\n", sqrt(ngb_treefind(P[i].PosPred , All.DesNumNgb, 0, P[i].Type, &ngblist, &r2list)));
		    }
		  
		  if(P[i].Left==0 || P[i].Right==0) 
		    {
		      if(P[i].Right==0 && P[i].NgbVelDisp<15 && NtypeLocal[P[i].Type]>All.DesNumNgb)
			{
			  P[i].HsmlVelDisp= sqrt(ngb_treefind(P[i].PosPred , All.DesNumNgb, 0, P[i].Type, &ngblist, &r2list));  
			}
		      else
			{
			  P[i].HsmlVelDisp=  P[i].HsmlVelDisp*( 0.5 + 0.5*pow(P[i].NgbVelDisp/((double)All.DesNumNgb), -1.0/3));
			}
		    }
		  else
		    {
		      P[i].HsmlVelDisp=  0.5*(P[i].Left + P[i].Right);
		    }
		}
	      
	      veldisp();

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


#endif









