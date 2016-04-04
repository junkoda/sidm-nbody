#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"
#include "domain.h"

#define TOLERANCE      0.001    /* try to load balance up to this level */
#define MAX_ITERATIONS 128     /* maxmimum number of iterations */
#define REDUC_FAC      0.995


static int    ntype[5], ntottype[5];
static int    groupShift;




/* This is the main routine for the domain decomposition. The code will 
 * decompose each particle type separately. It will try to balance the 
 * work-load as defined by the sum of the P[i]-GravCost in each domain.
 * The decomposition will respect the maximum memory-imbalance given
 * by the value of PartAllocFactor.
 *
 * The decomposition is done by performing an orthogonal recursive bisection,
 * using a hypercube communication model.
 */
void DomainDecomposition(void)
{
  int    i, type;
  double t0, t1;
  
  t0=second();
  
  if(ThisTask==0)
    printf("domain decomposition... \n");
  
  for(i=0; i<5; i++)
    ntype[i]=0;
  
  for(i=1; i<=NumPart; i++)
    ntype[P[i].Type]++;
  


  MPI_Allreduce(&ntype, &ntottype, 5, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  for(type=0; type<5; type++)
    {
      if(ntottype[type]>0)
	decomposeType(type);
    }

  if(ThisTask==0)
    printf("domain decomposition done. \n");

  t1=second();
  All.CPU_Domain+= timediff(t0,t1);

#ifdef DEB
  if(ThisTask==0)
    { 	
      FdDEB=fopen("domain.deb","w");
      fprintf(FdDEB," %d \n", All.NumCurrentTiStep);
      fclose(FdDEB);
    }
#endif
}



/*  This function does the domain decomposition for a single
 *  particle type 
 */
void decomposeType(int type) 
{
  int    level;
  int    axis, repflag;
  double posMin[3], posMax[3];
  double xsplit;
  int    sendTask, recvTask;
  int    imbalance;
  int    redo, redosum;

  if(ThisTask==0)
    printf("domain decompos. for type=%d\n", type);

  for(level=NTask,axis=0; level>1; level>>=1,axis++) 
    {
      if(axis>2)
	axis=0;

      findExtent(type, level, axis, &posMin[axis], &posMax[axis]);

      xsplit= findSplitPoint(type, level, posMin[axis], posMax[axis], axis, &imbalance);

      groupShift=0;

      do
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ ((level>>1) + groupShift); 

	  do  /* send and receive particles from other side  of domain */
	    {
	      if(sendTask < recvTask) 
		{
		  repflag = exchangeParticles_A(type, recvTask, axis, xsplit);     
		
		  posMax[axis]=xsplit;
		}
	      else 
		{
		  repflag = exchangeParticles_B(type, recvTask, axis, xsplit);
		  
		  posMin[axis]=xsplit;
		}
	      
	      if(repflag==0)
		break;
	    }
	  while(repflag==1);


	  if(repflag==2)
	    redo=1;
	  else
	    redo=0;

	  MPI_Allreduce(&redo, &redosum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      
	  groupShift++;
	  
	  if(ThisTask==0 && redosum>0)
	    printf("repeating with groupShift=%d\n",groupShift);
	}
      while(redosum>0 && groupShift<(level>>1));
    
      MPI_Barrier(MPI_COMM_WORLD);
      fflush(stdout);
    }
}



/*  This routine finds the split point for a processor group of
 *  size `level'. The split point is found by bisection, under 
 *  the constraint that the memory imbalance does not exceed
 *  a certain level.
 */
double findSplitPoint(int type, int level, double xleft, double xright, int axis, int *imbalanceflags)
{
  int    iter,flagmemoryimbalance,rightflag,leftflag;
  double xguess, balance, eps, xmaxleftofguess, xminrightofguess;
  int    pabove, pbelow, plimit;

  xguess = 0.5 * (xright + xleft);
  iter = 0;

  do
    {
      balance = combinedWork(type, level, axis, xguess, 
			     &xmaxleftofguess, &xminrightofguess, &flagmemoryimbalance,
			     &pabove, &pbelow, &plimit);
      
      if((pabove+pbelow)==0) /* nothing to split */
	break;

      eps = fabs(((double)(pbelow-plimit))/(pabove+pbelow));
      if(eps < TOLERANCE)
	{
	  if(pbelow<=plimit)
	    {
	      if(balance<(0.5-TOLERANCE/2))  /* still not enough work below . memory imbalance limits things */
		{
		  if(ThisTask == (ThisTask/level)*level)
		    {
		      /*
			printf("\nMemory imbalance limits work-load balancing. Task=%d\n", ThisTask);
			printf("\nbelow=%d above=%d limit=%d eps=%g balance=%g \n",pbelow,pabove,plimit,eps,balance);
			*/
		    }	     
		  break;
		}
	    }
	}

      eps = fabs(((double)(pabove-plimit))/(pabove+pbelow));
      if(eps < TOLERANCE)
	{
	  if(pabove<=plimit)
	    {
	      if(balance>(0.5+TOLERANCE/2))  /* still not enough work above . memory imbalance limits things */
		{
		  if(ThisTask == (ThisTask/level)*level)
		    {
		      /*
			printf("\nMemory imbalance limits work-load balancing. Task=%d\n", ThisTask);
			printf("\n below=%d above=%d limit=%d eps=%g balance=%g \n",
      			pbelow,pabove,plimit,eps,balance);
			*/
		    }	     
		  break;
		}
	    }
	}
      
      eps = fabs(balance - 0.5);

      if(pbelow>plimit)
	eps=2*TOLERANCE;   

      if(pabove>plimit)
	eps=2*TOLERANCE;   

      if(eps<TOLERANCE) /* work-load balance achieved */
	break;    
      
      if(xmaxleftofguess<=xleft &&  xminrightofguess>=xright) /* we are already inbetween two particles */
	{
	  /*
	    printf("\nno further improvement possible because we are inbetween two particles. iterations=%d balance=%g  task=%d\n",  
	             iter, balance, ThisTask);  
	  */
 	  break; 
	}
      
      if(flagmemoryimbalance&12)
	{
	  if(flagmemoryimbalance&4) /* lower side has too many particles */
	    {
	      xright = xguess;
	      rightflag=1;
	    }
	  
	  if(flagmemoryimbalance&8) /* upper side has too many particles */
	    {
	      xleft = xguess;
	      leftflag=1;
	    }
	}
      else
	{
	  if(balance > 0.5)
	    {
	      xright = xguess;
	      rightflag=0;
	    }
	  else 
	    {
	      xleft = xguess;
	      leftflag=0;
	    }
	}

      xguess = 0.5*(xleft + xright);

      iter++;
    }
  while(iter< MAX_ITERATIONS);
  
  if(iter>=MAX_ITERATIONS) 
    {
      if(ThisTask == (ThisTask/level)*level)
	{
	  printf("\nWarning. No convergence in `findSplitPoint' achieved.\nBalance: %g  Iterations: %d  Task: %d\n", 
		 balance, iter, ThisTask);
	  printf("flagmemoryimbalance=%d \n",flagmemoryimbalance);
	  printf("xmaxleftofguess=%g  xminrightofguess=%g   xleft=%g xright=%g\n",
		 xmaxleftofguess, xminrightofguess, xleft, xright);
	}
    }

  *imbalanceflags = flagmemoryimbalance;

  return xguess;
}



/*  This routine sums the work on both sides of the proposed split for 
 *  the particle group at the given level of the recursive bisection.
 *  The routine also reports the numbers of particles on both sides of the
 *  split and the closest particle coordinates to it.
 */
double combinedWork(int type, int level, int axis, double xsplit, double *xmaxleftofguess, 
		    double *xminrightofguess, int *flagmemoryimbalance,
		    int    *pabove, int *pbelow, int *plimit)
{
  int    i, j, masterTask;
  double work[2];
  double totworkabove=0, totworkbelow=0, workabove, workbelow;
  int    particlesabove, particlesbelow, par[3];
  int    parti[2], partj[2];
  int    flaglimit, sumlimit;
  int    startdummy=2;
  double xmaxleft, xminright;
  double maxmin[2];
  MPI_Status status;


  sumlimit=(level/2)*(1 + (int)(((double)ntottype[type])/NTask * (1.0 + (All.PartAllocFactor-1.0)*pow(REDUC_FAC,level/2))));

  /* `level' is equal to the number of Tasks in the group that is about to be split */
 
  masterTask = (ThisTask/level)*level;  /* master task in this group */
	
 
  getWork(type, axis, xsplit, &workabove, &workbelow, 
	  &particlesabove, &particlesbelow, &xmaxleft, &xminright);  /* Determine the work above 
                                                                            and below for this task */
  
  if(ThisTask == masterTask)  /* sum up the work of the group */ 
    {
      for(i=0, flaglimit=0; i<level/2; i++) 
	{
	  j = (i+masterTask) ^ (level>>1); /* partner task in particle exchange */

	  if(i==0)
	    {
	      totworkabove = workabove;  
	      totworkbelow = workbelow;
	      maxmin[0]= xmaxleft;
	      maxmin[1]= xminright;
	      parti[0]= particlesabove;
	      parti[1]= particlesbelow;
	    }
	  else
	    {
	      MPI_Send(&startdummy,   1, MPI_INT, masterTask+i, TAG_N, MPI_COMM_WORLD);

	      MPI_Recv(&work[0],   2, MPI_DOUBLE, masterTask+i, TAG_WORK, MPI_COMM_WORLD, &status);
	      MPI_Recv(&maxmin[0], 2, MPI_DOUBLE, masterTask+i, TAG_MAXMIN, MPI_COMM_WORLD, &status);
	      MPI_Recv(&parti[0],  2, MPI_INT   , masterTask+i, TAG_ANY, MPI_COMM_WORLD, &status);

	      totworkabove +=  work[0];
	      totworkbelow +=  work[1];
	      particlesabove+= parti[0];
	      particlesbelow+= parti[1];
	      if(xmaxleft < maxmin[0])   xmaxleft =maxmin[0];
	      if(xminright> maxmin[1])   xminright=maxmin[1];
	    }

	  MPI_Send(&startdummy,   1, MPI_INT, j, TAG_N, MPI_COMM_WORLD);

	  MPI_Recv(&work[0],   2, MPI_DOUBLE, j, TAG_WORK, MPI_COMM_WORLD, &status);
	  MPI_Recv(&maxmin[0], 2, MPI_DOUBLE, j, TAG_MAXMIN, MPI_COMM_WORLD, &status);
	  MPI_Recv(&partj[0], 2,  MPI_INT   , j, TAG_ANY, MPI_COMM_WORLD, &status);

	  totworkabove +=  work[0];
	  totworkbelow +=  work[1];
	  particlesabove+= partj[0];
	  particlesbelow+= partj[1];
	  if(xmaxleft < maxmin[0])   xmaxleft =maxmin[0];
	  if(xminright> maxmin[1])   xminright=maxmin[1];
	}


      if((particlesabove+particlesbelow) >= 2*sumlimit )
	{
	  printf("Both sides start out with more particles than sumlimit.\n");
	  printf("task: %d  level: %d  below: %d   above: %d  sumlimit: %d\n",
		 ThisTask,level, particlesbelow, particlesabove, sumlimit);

	  sumlimit=(particlesabove+particlesbelow+2)/2;

	  printf("I'll change sumlimit accordingly: new value=%d\n",sumlimit);
	  flaglimit|=16;
	  fflush(stdout);
	}
      

      if(particlesbelow > sumlimit) /* lower side has too many particles */
	flaglimit|=4;
	
      if(particlesabove > sumlimit) /* higher side has too many particles */
	flaglimit|=8;


      work[0]=   totworkabove;
      work[1]=   totworkbelow;

      maxmin[0]= xmaxleft;
      maxmin[1]= xminright;

      par[0]= particlesabove;
      par[1]= particlesbelow;
      par[2]= sumlimit;

      
      for(i=1; i<level; i++)
	{
	  MPI_Send(&work[0],   2, MPI_DOUBLE, masterTask+i, TAG_WORK, MPI_COMM_WORLD);
	  MPI_Send(&maxmin[0], 2, MPI_DOUBLE, masterTask+i, TAG_MAXMIN, MPI_COMM_WORLD);
	  MPI_Send(&flaglimit, 1, MPI_INT   , masterTask+i, TAG_ANY, MPI_COMM_WORLD);
	  MPI_Send(&par[0],    3, MPI_INT   , masterTask+i, TAG_N, MPI_COMM_WORLD);
  	}
    }
  else 
    {
      work[0]= workabove; 
      work[1]= workbelow;

      maxmin[0]= xmaxleft;
      maxmin[1]= xminright;

      parti[0]= particlesabove;
      parti[1]= particlesbelow; 

      MPI_Recv(&startdummy, 1, MPI_INT, masterTask, TAG_N, MPI_COMM_WORLD, &status);

      MPI_Send(&work[0],   2, MPI_DOUBLE, masterTask, TAG_WORK, MPI_COMM_WORLD);
      MPI_Send(&maxmin[0], 2, MPI_DOUBLE, masterTask, TAG_MAXMIN, MPI_COMM_WORLD);
      MPI_Send(&parti[0], 2,  MPI_INT,    masterTask, TAG_ANY, MPI_COMM_WORLD);

      MPI_Recv(&work[0], 2, MPI_DOUBLE, masterTask, TAG_WORK, MPI_COMM_WORLD, &status);
      MPI_Recv(&maxmin[0], 2, MPI_DOUBLE, masterTask, TAG_MAXMIN, MPI_COMM_WORLD, &status);
      MPI_Recv(&flaglimit, 1, MPI_INT, masterTask, TAG_ANY, MPI_COMM_WORLD, &status);
      MPI_Recv(&par[0],    3, MPI_INT   , masterTask, TAG_N, MPI_COMM_WORLD, &status);
       
      totworkabove=work[0];
      totworkbelow=work[1];
    }

 
  *xmaxleftofguess=maxmin[0];
  *xminrightofguess=maxmin[1];
  *flagmemoryimbalance=flaglimit;

  *pabove=par[0];
  *pbelow=par[1];
  *plimit=par[2];

  if((totworkabove + totworkbelow)>0)
    return(totworkbelow/(totworkabove + totworkbelow));
  else
    return 0.5;
}




/*  This function computes the work for the local domain on both sides of the proposed
 *  split, as well as the particle numbers, and closest coordinates.
 */
void getWork(int type, int axis, double xsplit, double *workabove, double *workbelow, 
	     int *particlesabove,int *particlesbelow, double *xmaxleft, double *xminright)
{
  int i;
  double above,below;
  int    partabove,partbelow;

  *xminright= MAX_REAL_NUMBER;
  *xmaxleft =-MAX_REAL_NUMBER;

  for(i=1, below=above=0, partbelow=partabove=0; i<=NumPart; i++)
    {
      if(P[i].Type == type)
	{
	  if(P[i].Pos[axis] < xsplit)
	    {
	      below+= P[i].GravCost; 
	      partbelow++;
	      
	      if(*xmaxleft <  P[i].Pos[axis])
		*xmaxleft = P[i].Pos[axis];
	    }
	  else
	    {
	      above+= P[i].GravCost; 
	      partabove++;

	      if(*xminright >  P[i].Pos[axis])
		*xminright = P[i].Pos[axis];
	    }
	}
    }
     
  *workabove=above;
  *workbelow=below;

  *particlesabove=partabove;
  *particlesbelow=partbelow;
}


/*  Ok, we have settled on a split. This function tries to get rid of all its 
 *  particles that are to the RIGHT of the split. It communicates with 
 *  another processor which runs exchangeParticles_B(). A is sending
 *  first to B, then it is receiving from it.
 *  The number of particles exchanged is negotiated such that the transfer
 *  is possible within the given memory constraints. 
 */
int exchangeParticles_A(int type, int recvTask, int axis, double xsplit)
{
  MPI_Status status;
  int N,Nrecv,i,ntobesent,maxsend;
  int ntobesentB,maxsendB,numpartB,numpartsphB;
  int ntype;
  int maxsend_old,maxsendB_old;
 

  /* all particles from this task which have a coordinate > xsplit
     are sent to the 'recvTask' */
  

  /* first we decide on the size of the exchanged chunks */

  for(i=1,ntobesent=ntype=0; i<=NumPart; i++)
    {
      if(P[i].Type == type)
	{
	  ntype++;

	  if(P[i].Pos[axis] > xsplit)
	    ntobesent++;
	}
    }

  MPI_Recv(&numpartB,   1, MPI_INT, recvTask, TAG_N, MPI_COMM_WORLD, &status);
  MPI_Recv(&ntobesentB, 1, MPI_INT, recvTask, TAG_ANY, MPI_COMM_WORLD, &status);
  MPI_Recv(&numpartsphB,1, MPI_INT, recvTask, TAG_N, MPI_COMM_WORLD, &status);

  MPI_Send(&NumPart,    1, MPI_INT, recvTask, TAG_N, MPI_COMM_WORLD);
  MPI_Send(&ntobesent,  1, MPI_INT, recvTask, TAG_ANY, MPI_COMM_WORLD);
  MPI_Send(&N_gas,      1, MPI_INT, recvTask, TAG_N, MPI_COMM_WORLD);

  
  maxsend=  imin(ntobesent,  All.BunchSizeDomain);
  maxsendB= imin(ntobesentB, All.BunchSizeDomain);

  do 
    {
      maxsend_old=  maxsend;
      maxsendB_old= maxsendB;

      maxsend=  imin(All.MaxPart - numpartB + ntobesentB,  maxsend);
      maxsendB= imin(All.MaxPart - NumPart  + ntobesent ,  maxsendB);
    }
  while((maxsend!=maxsend_old) || (maxsendB!=maxsendB_old));


  /* now make also sure that there is enough space for SPH particeles */
  if(type==0)
    {
      do 
	{
	  maxsend_old=  maxsend;
	  maxsendB_old= maxsendB;
	  
	  maxsend=  imin(All.MaxPartSph - numpartsphB + ntobesentB,  maxsend);
	  maxsendB= imin(All.MaxPartSph - N_gas       + ntobesent ,  maxsendB);
	}
      while((maxsend!=maxsend_old) || (maxsendB!=maxsendB_old));
    }


  for(i=1, N=0; N<maxsend && i<=NumPart; )
    {
      if(P[i].Type == type)
	{
	  if(P[i].Pos[axis] > xsplit)
	    {
	      if(type==0) /* special reorder routine for SPH particles (need to stay at beginning) */
		{
		  DomainPartBuf[N]= P_data[i-1];  /* copy particle and collect in contiguous memory */
		  DomainSphBuf[N]=  SphP_data[i-1];
		  
		  P_data[i-1]= P_data[N_gas-1];
		  P_data[N_gas-1]= P_data[NumPart-1];
		      
		  SphP_data[i-1]= SphP_data[N_gas-1];
		  
		  N_gas--;
		}
	      else
		{
		  DomainPartBuf[N]=P_data[i-1];  /* copy particle and collect in contiguous memory */
		  P_data[i-1]=P_data[NumPart-1];
		}
		  
	      N++;
	      NumPart--;
	    }
	  else
	    {
	      i++;
	    }
	}
      else
	{
	  i++;
	}
    }


  
  /* transmit */

  MPI_Send(&N, 1, MPI_INT, recvTask, TAG_N, MPI_COMM_WORLD);
  
  if(N>0)
    {
      MPI_Ssend(&DomainPartBuf[0], N*sizeof(struct particle_data), MPI_BYTE, recvTask, TAG_PDATA, MPI_COMM_WORLD);
      
      if(type==0)
	MPI_Ssend(&DomainSphBuf[0], N*sizeof(struct sph_particle_data), MPI_BYTE, recvTask, TAG_ANY, MPI_COMM_WORLD);
    }
  

  MPI_Recv(&Nrecv, 1, MPI_INT, recvTask, TAG_N, MPI_COMM_WORLD, &status);
  
  if(Nrecv>0)
    {
      if(type==0)
	{
	  for(i=0; i<Nrecv; i++)
	    {
	      P_data[NumPart+i]= P_data[N_gas+i];
	    }
	  MPI_Recv(&P_data[N_gas], Nrecv*sizeof(struct particle_data), MPI_BYTE, recvTask, TAG_PDATA, MPI_COMM_WORLD, &status);
	  MPI_Recv(&SphP_data[N_gas], Nrecv*sizeof(struct sph_particle_data), MPI_BYTE, recvTask, TAG_ANY, MPI_COMM_WORLD, &status);

	  N_gas+= Nrecv;
	}
      else
	{
	  MPI_Recv(&P_data[NumPart], Nrecv*sizeof(struct particle_data), MPI_BYTE, recvTask, TAG_PDATA, MPI_COMM_WORLD, &status);

	}

      NumPart+=Nrecv;
    }


  if(N<ntobesent || Nrecv<ntobesentB)
    {
      if(N>0 || Nrecv>0)
	return 1;
      else
	return 2;
    }
  else
    {
      return 0;
    }
}



/*  Ok, we have settled on a split. This function tries to get rid of all its 
 *  particles that are to the LEFT of the split. It communicates with 
 *  another processor which runs exchangeParticles_A(). A is sending
 *  first to B, so B first receives and then sends itself.
 *  The number of particles exchanged is negotiated such that the transfer
 *  is possible within the given memory constraints. 
 */
int exchangeParticles_B(int type, int recvTask, int axis, double xsplit)
{
  MPI_Status status;
  int N,Nrecv,i,ntobesent,maxsend;
  int ntobesentA,maxsendA,numpartA,numpartsphA;
  int ntype;
  int maxsend_old,maxsendA_old;


  /* all particles from this task which have a coordinate > xsplit
     are sent to the 'recvTask' */
 

  /* first we decide on the size of the exchanged chunks */

  for(i=1,ntobesent=ntype=0; i<=NumPart; i++)
    {
      if(P[i].Type == type)
	{
	  ntype++;

	  if(P[i].Pos[axis] < xsplit)
	    ntobesent++;
	}
    }

  MPI_Send(&NumPart,    1, MPI_INT, recvTask, TAG_N, MPI_COMM_WORLD);
  MPI_Send(&ntobesent,  1, MPI_INT, recvTask, TAG_ANY, MPI_COMM_WORLD);
  MPI_Send(&N_gas,      1, MPI_INT, recvTask, TAG_N, MPI_COMM_WORLD);

  MPI_Recv(&numpartA,   1, MPI_INT, recvTask, TAG_N, MPI_COMM_WORLD, &status);
  MPI_Recv(&ntobesentA, 1, MPI_INT, recvTask, TAG_ANY, MPI_COMM_WORLD, &status);
  MPI_Recv(&numpartsphA,1, MPI_INT, recvTask, TAG_N, MPI_COMM_WORLD, &status);

  
  maxsend=  imin(ntobesent,  All.BunchSizeDomain);
  maxsendA= imin(ntobesentA, All.BunchSizeDomain);

  do 
    {
      maxsend_old=maxsend;
      maxsendA_old=maxsendA;

      maxsend=  imin(All.MaxPart - numpartA + ntobesentA,  maxsend);
      maxsendA= imin(All.MaxPart - NumPart  + ntobesent ,  maxsendA);
    }
  while((maxsend!=maxsend_old) || (maxsendA!=maxsendA_old));

  /* now make also sure that there is enough space for SPH particeles */
  if(type==0)
    {
      do 
	{
	  maxsend_old=  maxsend;
	  maxsendA_old= maxsendA;
	  
	  maxsend=  imin(All.MaxPartSph - numpartsphA + ntobesentA,  maxsend);
	  maxsendA= imin(All.MaxPartSph - N_gas       + ntobesent ,  maxsendA);
	}
      while((maxsend!=maxsend_old) || (maxsendA!=maxsendA_old));
    }

 
  for(i=1, N=0; N<maxsend && i<=NumPart; )
    {
      if(P[i].Type == type)
	{
	  if(P[i].Pos[axis] < xsplit)
	    {
	      if(type==0) /* special reorder routine for SPH particles (need to stay at beginning) */
		{
		  DomainPartBuf[N]= P_data[i-1];  /* copy particle and collect in contiguous memory */
		  DomainSphBuf[N]=  SphP_data[i-1];
		  
		  P_data[i-1]= P_data[N_gas-1];
		  P_data[N_gas-1]= P_data[NumPart-1];
		  
		  SphP_data[i-1]= SphP_data[N_gas-1];
		  
		  N_gas--;
		}
	      else
		{
		  DomainPartBuf[N]=P_data[i-1];  /* copy particle and collect in contiguous memory */
		  P_data[i-1]=P_data[NumPart-1];
		}

	      N++;
	      NumPart--;
	    }
	  else
	    {
	      i++;
	    }
	}
      else
	{
	  i++;
	}
    }


  
  /* transmit */

  MPI_Recv(&Nrecv, 1, MPI_INT, recvTask, TAG_N, MPI_COMM_WORLD, &status);
  
  if(Nrecv>0)
    {
      if(type==0)
	{
	  for(i=0; i<Nrecv; i++)
	    {
	      P_data[NumPart+i]= P_data[N_gas+i];
	    }
	  MPI_Recv(&P_data[N_gas], Nrecv*sizeof(struct particle_data), MPI_BYTE, recvTask, TAG_PDATA, MPI_COMM_WORLD, &status);
	  MPI_Recv(&SphP_data[N_gas], Nrecv*sizeof(struct sph_particle_data), MPI_BYTE, recvTask, TAG_ANY, MPI_COMM_WORLD, &status);
	  
	  N_gas+= Nrecv;
	}
      else
	{
	  MPI_Recv(&P_data[NumPart], Nrecv*sizeof(struct particle_data), MPI_BYTE, recvTask, TAG_PDATA, MPI_COMM_WORLD, &status);
	}

      NumPart+=Nrecv;
    }


  MPI_Send(&N, 1, MPI_INT, recvTask, TAG_N, MPI_COMM_WORLD);
  
  if(N>0)
    {
      MPI_Ssend(&DomainPartBuf[0], N*sizeof(struct particle_data), MPI_BYTE, recvTask, TAG_PDATA, MPI_COMM_WORLD);

      if(type==0)
	MPI_Ssend(&DomainSphBuf[0], N*sizeof(struct sph_particle_data), MPI_BYTE, recvTask, TAG_ANY, MPI_COMM_WORLD);
    }


  if(N<ntobesent || Nrecv<ntobesentA)
    {
      if(N>0 || Nrecv>0)
	return 1;
      else
	return 2;
    }
  else
    {
      return 0;
    }
}



/* This routine finds the extent (i.e. an enclosing rectangle) for the 
 * particles stored by a group of processors at a given level in the hierarchical
 * domain decomposition
 */
void findExtent(int type, int level, int axis,double *xleft, double *xright)
{
  int      i, masterTask;
  int      startdummy=3;
  double   xmin,xmax;
  double   xx[2];
  MPI_Status status;

  xmin= MAX_REAL_NUMBER;
  xmax=-MAX_REAL_NUMBER;


  for(i=1; i<=NumPart; i++)
    if(P[i].Type == type)
      {
	if(xmin > P[i].Pos[axis])
	  xmin=P[i].Pos[axis];
	
	if(xmax < P[i].Pos[axis])
	  xmax=P[i].Pos[axis];
      }

  if(level>1)
    {
      masterTask = (ThisTask/level)*level;  /* master task in this group */
  
      if(ThisTask==masterTask)
	{
	  for(i=1;i<level; i++) 
	    {
	      MPI_Send(&startdummy,   1, MPI_INT, i+masterTask,TAG_ANY, MPI_COMM_WORLD);	  
	      
	      MPI_Recv(&xx[0],   2, MPI_DOUBLE, i+masterTask, TAG_WORK, MPI_COMM_WORLD, &status);
	      if(xx[0]<xmin)
		xmin=xx[0];
	      if(xx[1]>xmax)
		xmax=xx[1];
	    }
	  
	  xx[0]=xmin;
	  xx[1]=xmax;

	  for(i=1;i<level; i++) 
	    {
	      MPI_Send(&xx[0],   2, MPI_DOUBLE, i+masterTask, TAG_N, MPI_COMM_WORLD);
	    }
	}
      else
	{
	  xx[0]=xmin;
	  xx[1]=xmax;

	  MPI_Recv(&startdummy,   1, MPI_INT, masterTask,  TAG_ANY, MPI_COMM_WORLD,&status);
	  
	  MPI_Send(&xx[0],   2, MPI_DOUBLE, masterTask,TAG_WORK, MPI_COMM_WORLD);
	  MPI_Recv(&xx[0],   2, MPI_DOUBLE, masterTask,  TAG_N, MPI_COMM_WORLD,&status);
	  
	  xmin=xx[0];
	  xmax=xx[1];
	}
    }

  *xleft=xmin;
  *xright=xmax;
}





































