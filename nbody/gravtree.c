#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"



/*  This function computes the gravitational forces for all active particles.
 *  A new tree is constructed, if the number of force computations since
 *  it's last construction exceeds some fraction of the total
 *  particle number, otherwise tree nodes are dynamically updated if needed.
 */
void gravity_tree(void)
{
  int     i,j,k;
  int     level,sendTask,recvTask;
  double  s_a,s_a_inverse,dt,dt_h0;
  int     ntot,ntotleft,npleft,nthis;     
  int     nstart,nbuffer,ncount,nchunk;
  int     timelinecounter,startcounter;   
  int     *nrecv,*noffset;
  double  tstart,tend,timetree=0;
  double  fac1,fac2,fac3,a,a2;  
  static int numnodes;
  MPI_Status status;
#ifdef DIAG
  int     *numnodeslist,maxnumnodes;
  int     costtotal,*costtreelist,costtotal_quadru,*costtreelist_quadru;
  double  maxt,sumt,*timetreelist;
  double  fac,plb,plb_max;
  int     tot_nodeupdates,tot_nodeupdate_particles;
#endif



  if(All.ComovingIntegrationOn)
    {
      set_softenings();   /* set new softening lengths */

      s_a_inverse=1/(All.Hubble*sqrt(All.Omega0 + All.Time*(1-All.Omega0-All.OmegaLambda)+ All.Time*All.Time*All.Time*All.OmegaLambda));
    }
  else
    s_a_inverse=1;
    



  /* NumForceUpdate = number of particles on this processor that want a force update */
  
  MPI_Allreduce(&NumForceUpdate, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  All.NumForcesSinceLastDomainDecomp     += ntot;
  All.NumForcesSinceLastTreeConstruction += ntot;



  tstart=second();
  if(All.NumForcesSinceLastTreeConstruction >= All.TreeUpdateFrequency*All.TotNumPart)
    {
      if(NoCostFlag==0)   
	force_costevaluate(); /* assign accumulated cost to particles */ 

      if(ThisTask==0)
	printf("Tree construction.\n");

      predict_collisionless_only(All.Time); 
      /* need to do the prediction here since its usually not done for all particles  */
      /* but: SPH particles have already been predicted, and if the density was
         computed before, some of the SPH rates have already been updated */

      numnodes= force_treebuild();

      All.NumForcesSinceLastTreeConstruction=0;

      if(ThisTask==0)
	printf("Tree construction done.\n");
	
      NoCostFlag=0;
    }
  else
    {
#ifdef VELDISP
      predict_collisionless_only(All.Time); 
#endif
#ifdef SIDM
      // Is this necessary? What is it's computational cost?
      predict_collisionless_only(All.Time); 
#endif

      ngb_update_nodes();  /* need to make sure that we have outer boundaries for all nodes */
    }

  tend=second();
  All.CPU_TreeConstruction+= timediff(tstart,tend);


  
  ntotleft= ntot;              /* particles left for all tasks together */
  npleft= NumForceUpdate;      /* particles left for this task */

  nstart= IndFirstUpdate;      /* first particle for this task */
  startcounter=0;


  nthis = ((double)NumForceUpdate*(All.BunchSizeForce-NTask))/ntot + 1; /* maximum size of chunck coming from this task */

  if(NumForceUpdate>0) 
    {
      nchunk= NumForceUpdate/nthis;                             /* number of chunks needed */  
      if((NumForceUpdate%nthis)>0) nchunk+=1;
    }
  else /* note: it could happen that this processor does have NumForceUpdate==0, while ntot>0 */
    nchunk=0;


  nrecv=   malloc(sizeof(int)*NTask);  /* list of particle numbers that constituate current bunch */
  noffset= malloc(sizeof(int)*NTask);  /* offsets of bunches in common list */


  force_resetcost();  /* resets counters for statistics of number of particle-node interactions */
  
  while(ntotleft>0)
    {
      if(nthis>npleft)
	nthis=npleft;


      for(i=nstart, ncount=0, timelinecounter=startcounter; 
	  timelinecounter<NumForceUpdate && ncount<nthis; 
	  i=P[i].ForceFlag, timelinecounter++, ncount++);

	
      MPI_Allgather(&ncount, 1, MPI_INT, nrecv, 1, MPI_INT, MPI_COMM_WORLD);


      for(i=nbuffer=0; i<NTask; i++)  /* compute length of common list */
	nbuffer += nrecv[i];

      for(i=1,noffset[0]=0; i<NTask; i++) /* set-up offset-table */
	noffset[i]=noffset[i-1]+nrecv[i-1];


      /* fill in the own particles at the right place */
      for(i=nstart, ncount=0, timelinecounter=startcounter; 
	  timelinecounter<NumForceUpdate && ncount<nthis; 
	  i=P[i].ForceFlag, timelinecounter++, ncount++)
	{
	  dt = (All.Time - P[i].CurrentTime);  
	  dt_h0 = dt*s_a_inverse;  
	  
	  /* need to predict the particle, since it may not have been done yet */
	  /* note: the predicted velocity will be needed in the comoving integration (below) */

	  for(k=0;k<3;k++)
	    {
	      GravDataIn[noffset[ThisTask] + ncount].Pos[k]= P[i].PosPred[k] = P[i].Pos[k] + P[i].Vel[k]*dt_h0;  
	      P[i].VelPred[k]= P[i].Vel[k] + P[i].Accel[k]*dt;
	    }
	  GravDataIn[noffset[ThisTask]+ncount].Type= P[i].Type;
	  GravDataIn[noffset[ThisTask]+ncount].OldAcc= P[i].OldAcc;
	}

      /* now start big communication */

      tstart=second();      
      for(level=1;level<NTask;level++)
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ level; 

	  MPI_Sendrecv(&GravDataIn[noffset[sendTask]], nrecv[sendTask]*sizeof(struct gravdata_in), MPI_BYTE, recvTask, TAG_ANY,
		       &GravDataIn[noffset[recvTask]], nrecv[recvTask]*sizeof(struct gravdata_in), MPI_BYTE, recvTask, TAG_ANY ,
		       MPI_COMM_WORLD,&status);
	}
      tend=second();
      All.CPU_CommSum+= timediff(tstart,tend);



      /* ok, all preparations are done. Now force evaluation. */

      tstart=second();

      for(i=0; i<nbuffer; i++)
	force_treeevaluate(i, s_a_inverse);
      
      tend=second();

      All.CPU_TreeWalk+= timediff(tstart, tend);
      timetree += timediff(tstart, tend);


      tstart=second();
      MPI_Barrier(MPI_COMM_WORLD);
      tend=second();
      All.CPU_Imbalance+= timediff(tstart, tend);



      /* now communicate force contributions, and sum up */

      tstart=second();
      for(level=1; level<NTask; level++)
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ level; 

	  MPI_Sendrecv(&GravDataResult[noffset[recvTask]], nrecv[recvTask]*sizeof(struct gravdata_out), MPI_BYTE, recvTask, TAG_ANY,
		       &GravDataPartialResult[0],  nrecv[sendTask]*sizeof(struct gravdata_out), MPI_BYTE, recvTask, TAG_ANY,
		       MPI_COMM_WORLD, &status);
	
	  for(i=0; i<nrecv[ThisTask]; i++)
	    {
	      for(k=0;k<3;k++)
		GravDataResult[noffset[ThisTask] + i].Acc[k] += GravDataPartialResult[i].Acc[k];
	    }
	}
      tend=second();
      All.CPU_CommSum+= timediff(tstart,tend);



      /* transfer the result to the particles */

      for(i=nstart, ncount=0, timelinecounter=startcounter; 
	  timelinecounter<NumForceUpdate && ncount<nrecv[ThisTask]; 
	  i=P[i].ForceFlag, timelinecounter++, ncount++)
	{
	  for(k=0;k<3;k++)
	    {
	      P[i].Accel[k] = GravDataResult[noffset[ThisTask]+ncount].Acc[k];
	    }
	}

      nstart=i;                  /* continue at this point in the timeline in next iteration */
      startcounter=timelinecounter;


      npleft-= nrecv[ThisTask];
      ntotleft-= nbuffer;
    }



  /* now add things for comoving integration */

  if(All.ComovingIntegrationOn)
    {
      if(All.TypeOfOpeningCriterion==1) 
	{
	  fac3= 0.5*All.Hubble*All.Hubble*All.Omega0 / All.G;
	  for(i=IndFirstUpdate,timelinecounter=0; 
	      timelinecounter<NumForceUpdate;
	      i=P[i].ForceFlag, timelinecounter++)
	    {
	      for(j=0,a2=0;j<3;j++)
		{
#ifdef PERIODIC
		  a = P[i].Accel[j];
#else
		  a = P[i].Accel[j]+ fac3*P[i].PosPred[j];
#endif
		  a2 += a*a;
		}
	      P[i].OldAcc= sqrt(a2);
	    }
	}



      s_a=sqrt(All.Omega0 + All.Time*(1-All.Omega0-All.OmegaLambda)+ All.Time*All.Time*All.Time*All.OmegaLambda);

      fac1= All.G/ ( All.Hubble * All.Time * All.Time* s_a );

      fac2= -1.5/All.Time;

      fac3= 0.5*All.Hubble * All.Omega0 / ( All.Time * All.Time * s_a );
	
      for(i=IndFirstUpdate,timelinecounter=0; 
	  timelinecounter<NumForceUpdate;
	  i=P[i].ForceFlag, timelinecounter++)
	{
	  for(j=0;j<3;j++)
#ifdef PERIODIC
	    P[i].Accel[j] = fac1*P[i].Accel[j] 
	                   + fac2*P[i].VelPred[j];
#else
	    P[i].Accel[j] = fac1*P[i].Accel[j] 
	                   + fac2*P[i].VelPred[j] 
      	                   + fac3*P[i].PosPred[j];  
#endif
	}
    }
  else  /* else muliply by G */
    {
      if(All.TypeOfOpeningCriterion==1) 
	for(i=IndFirstUpdate,timelinecounter=0; 
	    timelinecounter<NumForceUpdate;
	    i=P[i].ForceFlag, timelinecounter++)
	  {
	    P[i].OldAcc= sqrt(P[i].Accel[0]*P[i].Accel[0] + 
			      P[i].Accel[1]*P[i].Accel[1] +
			      P[i].Accel[2]*P[i].Accel[2] );
	  }


      fac1= All.OmegaLambda*All.Hubble*All.Hubble;  
            /* this factor allows a computation of cosmological simulation 
               with vacuum energy in physical coordinates */

      for(i=IndFirstUpdate,timelinecounter=0; 
	  timelinecounter<NumForceUpdate;
	  i=P[i].ForceFlag,timelinecounter++)
	{
	    for(j=0;j<3;j++)
	      P[i].Accel[j] = All.G*P[i].Accel[j] + fac1*P[i].PosPred[j];
	    
	}
    }
  


  /* Now the force computation is finished */


#ifdef DIAG
  /*  gather some diagnostic information */

  tstart=second();

  timetreelist=malloc(sizeof(double)*NTask);
  costtreelist=malloc(sizeof(int)*NTask);
  costtreelist_quadru=malloc(sizeof(int)*NTask);
  numnodeslist=malloc(sizeof(int)*NTask);

  costtotal= force_getcost_single();
  costtotal_quadru= force_getcost_quadru();

  MPI_Gather(&costtotal, 1, MPI_INT, costtreelist, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&costtotal_quadru, 1, MPI_INT, costtreelist_quadru, 1, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Gather(&timetree, 1, MPI_DOUBLE, timetreelist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&NumPart, 1, MPI_INT, nrecv, 1, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Gather(&numnodes, 1, MPI_INT, numnodeslist, 1, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Reduce(&Num_nodeupdates, &tot_nodeupdates, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&Num_nodeupdate_particles, &tot_nodeupdate_particles, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);


  if(ThisTask==0)
    { 
      All.TotNumOfForces+= ntot;

      fprintf(FdTimings,"Step= %d  t= %g  dt= %g \n",All.NumCurrentTiStep,All.Time,All.TimeStep);
      fprintf(FdTimings,"Nf= %d  total-Nf= %d\n", ntot, All.TotNumOfForces);
      fprintf(FdTimings,"nodeupdates: %d (%d part)\n",tot_nodeupdates,tot_nodeupdate_particles);


      fac=NTask/((double) All.TotNumPart);

      for(i=0, maxt=timetreelist[0], sumt=0, plb_max=0,maxnumnodes=0, costtotal=costtotal_quadru=0; i<NTask; i++)
	{
	  costtotal+= costtreelist[i];
          costtotal_quadru+= costtreelist_quadru[i];

	  if(maxt<timetreelist[i])
	    maxt=timetreelist[i];
	  sumt+=timetreelist[i];

	  plb = nrecv[i] * fac;

	  if(plb>plb_max)
	    plb_max=plb;

	  if(numnodeslist[i]>maxnumnodes)
	    maxnumnodes=numnodeslist[i];
	}
      fprintf(FdTimings,"work-load balance: %g  max=%g avg=%g PE0=%g\n", maxt/(sumt/NTask), maxt, sumt/NTask,  timetreelist[0]);
      fprintf(FdTimings,"particle-load balance: %g\n", plb_max);
      fprintf(FdTimings,"max. nodes: %d, filled: %g\n", maxnumnodes, maxnumnodes/(All.TreeAllocFactor*All.MaxPart));
      fprintf(FdTimings,"part/sec=%g | %g  ia/part=%g (QP-frac %g)\n", 
                        ntot/(sumt+1.0e-20) , ntot/(maxt*NTask), 
                        ((double)(costtotal+costtotal_quadru))/ntot, costtotal_quadru/((double)(costtotal+costtotal_quadru)));
      fprintf(FdTimings,"\n");

      fflush(FdTimings);
    }
 

  free(numnodeslist); 
  free(costtreelist);
  free(costtreelist_quadru);
  free(timetreelist);

  tend=second();
  All.CPU_Diagnostic+= timediff(tstart, tend);

#endif


  free(noffset);
  free(nrecv);


#ifdef DEB
  if(ThisTask==0)
    { 	
      FdDEB=fopen("gravtree.deb","w");
      fprintf(FdDEB," %d \n",All.NumCurrentTiStep);
      fclose(FdDEB);
    }
#endif
}






/*  This function sets the (comoving) softening length of 
 *  all particle species in the table All.SofteningTable[...]
 *  We check here that the proper softening length is bounded
 *  by the ..MaxPhys values.
 */
void set_softenings(void)
{
  if(All.SofteningGas*All.Time > All.SofteningGasMaxPhys)
    All.SofteningTable[0] = All.SofteningGasMaxPhys/All.Time;
  else
    All.SofteningTable[0] = All.SofteningGas;
  
  if(All.SofteningHalo*All.Time > All.SofteningHaloMaxPhys)
    All.SofteningTable[1] = All.SofteningHaloMaxPhys/All.Time;
  else
    All.SofteningTable[1] = All.SofteningHalo;
  
  if(All.SofteningDisk*All.Time > All.SofteningDiskMaxPhys)
    All.SofteningTable[2] = All.SofteningDiskMaxPhys/All.Time;
  else
    All.SofteningTable[2] = All.SofteningDisk;

  if(All.SofteningBulge*All.Time > All.SofteningBulgeMaxPhys)
    All.SofteningTable[3] = All.SofteningBulgeMaxPhys/All.Time;
  else
    All.SofteningTable[3] = All.SofteningBulge;

  if(All.SofteningStars*All.Time > All.SofteningStarsMaxPhys)
    All.SofteningTable[4] = All.SofteningStarsMaxPhys/All.Time;
  else
    All.SofteningTable[4] = All.SofteningStars;

  All.MinGasHsml= All.MinGasHsmlFractional*All.SofteningTable[0]; 
}






