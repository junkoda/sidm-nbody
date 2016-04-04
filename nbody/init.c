#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"

/*
 *  init() reads in the initial conditions,
 *  and allocates storage for the tree(s).
 *  An intial domain decomposition and force
 *  computation is done, followed by a second
 *  domain decomposition based on the initial
 *  work estimates. Then the first particle timesteps
 *  are determined. The simulation is set up for
 *  the timestep iteration in run().
 */
void init(void)
{
  int i,j;

  All.Time = All.TimeBegin;   

  switch(All.ICFormat)
    {
    case 1:
      read_ic(All.InitCondFile); 
      break;
    case 3:
      if(RestartFlag==2)
	read_ic(All.InitCondFile); 
      else
	read_ic_cluster(All.InitCondFile);
      break;
    default:
      printf("ICFormat=%d not supported.\n", All.ICFormat);
      endrun(0);
    }

  All.SofteningTable[0] = All.SofteningGas;
  All.SofteningTable[1] = All.SofteningHalo;
  All.SofteningTable[2] = All.SofteningDisk;
  All.SofteningTable[3] = All.SofteningBulge;
  All.SofteningTable[4] = All.SofteningStars;
  All.MinGasHsml= All.MinGasHsmlFractional*All.SofteningTable[0]; 



  All.NumCurrentTiStep=0;    /* setup some counters */
  All.SnapshotFileCount=0;
  if(RestartFlag==2)
    All.SnapshotFileCount= atoi(All.InitCondFile+strlen(All.InitCondFile)-3)+1;
 
  All.TotNumOfForces=0;
  All.NumForcesSinceLastDomainDecomp=0;

  if(All.ComovingIntegrationOn)
    if(All.PeriodicBoundariesOn==1)
      check_omega();

  All.TimeLastStatistics= All.TimeBegin-All.TimeBetStatistics;

  if(RestartFlag==2)
    {
      if(All.ComovingIntegrationOn)
	All.TimeOfFirstSnapshot = All.Time*All.TimeBetSnapshot;
      else
	All.TimeOfFirstSnapshot = All.Time+All.TimeBetSnapshot;
    }

#ifdef SIDM
  vmax= getvmax();
#endif

  for(i=1; i<=NumPart; i++)    /*  start-up initialization */
    {
      for(j=0;j<3;j++)
	{
	  P[i].PosPred[j]=P[i].Pos[j];
	  P[i].VelPred[j]=P[i].Vel[j];
	  P[i].Accel[j]=0;
#ifdef SIDM
	  P[i].dVel[j]= 0; // Oct 9, 2005
#endif
	}
      P[i].OldAcc=0;
      P[i].GravCost=1; 
      P[i].Potential=0; 
#ifdef STELLARAGE
      if(RestartFlag==0)
	P[i].MeanStellarAge= 0;
#endif
    }
    
  for(i=1; i<=N_gas; i++)     /* initialize sph_properties */
    {
      SphP[i].EgySpecPred  = SphP[i].EgySpec;
      SphP[i].DtEgySpec=0;
      SphP[i].DtDensity = 0;
      SphP[i].DtHsml = 0;
      if(RestartFlag==0)
	{
	  SphP[i].Hsml=0;
	  SphP[i].Density= 0;
#ifdef COOLING
	  SphP[i].Ne= 1.0;
#endif
#ifdef SFR
	  SphP[i].FormedStellarMass= 0;
#ifdef CLOUDS
	  SphP[i].CloudMass= 0;
#endif
#endif
	}
    }
    
  if(All.ComovingIntegrationOn==0)
    compute_global_quantities_of_system();    
    

  force_treeallocate(All.TreeAllocFactor*All.MaxPart, All.MaxPart);

  DomainDecomposition();  /* do initial domain decomposition (gives equal numbers of particles) */

  for(i=1; i<=NumPart; i++)  
    {
      P[i].CurrentTime=All.TimeBegin;
      P[i].ForceFlag=i+1;
    }
  P[NumPart].ForceFlag=1; IndFirstUpdate=1; NumForceUpdate=NumPart; NumSphUpdate=N_gas;

  ngb_treeallocate(MAX_NGB);
  ngb_treebuild(); 
  determine_interior();

  setup_smoothinglengths(All.DesNumNgb);

#ifdef VELDISP
  for(i=1; i<=NumPart; i++)  
   {
      P[i].CurrentTime=All.TimeBegin;
      P[i].ForceFlag=i+1;
    }
  P[NumPart].ForceFlag=1; IndFirstUpdate=1; NumForceUpdate=NumPart; NumSphUpdate=N_gas;

  setup_smoothinglengths_veldisp(All.DesNumNgb);
#endif

#ifdef SIDM
  for(i=1; i<=NumPart; i++)  
   {
      P[i].CurrentTime=All.TimeBegin;
      P[i].ForceFlag=i+1;
    }
  P[NumPart].ForceFlag=1; IndFirstUpdate=1; NumForceUpdate=NumPart; NumSphUpdate=N_gas;

  setup_smoothinglengths_sidm(All.DesNumNgb);
#endif

  for(i=1; i<=NumPart; i++) /* redo time-line */
    {
      P[i].CurrentTime= P[i].MaxPredTime= All.TimeBegin;
      P[i].ForceFlag=i+1;
    }
  P[NumPart].ForceFlag=1; IndFirstUpdate=1; NumForceUpdate=NumPart; NumSphUpdate=N_gas;



  All.NumForcesSinceLastTreeConstruction= All.TreeUpdateFrequency*All.TotNumPart ; /* ensures that new tree will be constructed */
  NoCostFlag=1;

  compute_accelerations(1);       /* ... compute accelerations */


  force_costevaluate();

  DomainDecomposition();         /* try to balance the initial work */

  All.NumForcesSinceLastTreeConstruction= All.TreeUpdateFrequency*All.TotNumPart ; /* ensures that new tree will be constructed */
  NoCostFlag=1;


  for(i=1; i<=NumPart; i++) 
      P[i].ForceFlag=i+1;
  P[NumPart].ForceFlag=1; IndFirstUpdate=1; NumForceUpdate=NumPart; NumSphUpdate=N_gas;

  find_timesteps(2);              /* ... set particels to initial timesteps   */
                                  /* ordered timeline will be constructed there */

  compute_global_quantities_of_system();
  SysStateAtStart = SysState;     /* ... remember global initial state */
}


/* This routine computes the mass content of the box and
 * compares it to the specified value of Omega.
 * If discrepant, the run is terminated.
 */
void check_omega(void)
{
  double mass=0, masstot, omega;
  int    i;

  for(i=1; i<=NumPart; i++)
    mass+= P[i].Mass;


  MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  omega= masstot/(All.BoxSize*All.BoxSize*All.BoxSize)/ (3*All.Hubble*All.Hubble/(8*PI*All.G));

  if(fabs(omega-All.Omega0) > 1.0e-3)
    {
      if(ThisTask==0)
	{
	  printf("\n\nI've found something odd!\n");
	  printf("The mass content accounts only for Omega=%g,\nbut you specified Omega=%g in the parameterfile.\n", 
		 omega, All.Omega0);
	  printf("\nI better stop.\n");
	}
      endrun(0);
    }
}



/*
 *  This function is used to find an initial smoothing length for each SPH
 *  particle. It guarantees that the number of neighbours will be
 *  between desired_ngb-MAXDEV and desired_ngb+MAXDEV
 */
void setup_smoothinglengths(int desired_ngb)
{
  int    i, ntot, last=0;
  float *r2list;
  int    *ngblist, count;
  int    iter=0;

#ifdef VELDISP
#define PPP P
#else
#ifdef SIDM
#define PPP P
#else
#define PPP SphP
#endif
#endif


  for(i=1; i<=N_gas; i++) 
    { 
      SphP[i].Hsml= sqrt(ngb_treefind( P[i].PosPred, desired_ngb, 0, 0, &ngblist, &r2list));  
      PPP[i].Right= PPP[i].Left=0;
    } 

  density();

  /* we will here store lower
   * and upper boundary of intervals for smoothing range bisection
   * in SphP[i].Left and SphP[i].Right 
   */

  do
    {
      for(i=1, NumForceUpdate=NumSphUpdate=0; i<=N_gas; i++) 
	{ 
	  if(SphP[i].NumNgb < (All.DesNumNgb-All.MaxNumNgbDeviation ) || SphP[i].NumNgb > (All.DesNumNgb+All.MaxNumNgbDeviation ))
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
		PPP[i].Left= SphP[i].Hsml;

	      if(SphP[i].NumNgb > (All.DesNumNgb+All.MaxNumNgbDeviation))
		PPP[i].Right= SphP[i].Hsml;
	    }
	}
      
      
      MPI_Allreduce(&NumSphUpdate, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      
      if(ntot>0)
	{
	  if(ThisTask==0)
	    {
	      printf("\nNGB iteration %d.  still %d particles\n", iter, ntot);
	      fflush(stdout);
	    }     

	  for(i=IndFirstUpdate,count=0; count<NumForceUpdate; i=P[i].ForceFlag, count++)   
	    {
	      if(PPP[i].Left==0 || PPP[i].Right==0) 
		SphP[i].Hsml=  SphP[i].Hsml*( 0.5 + 0.5*pow(SphP[i].NumNgb/((double)All.DesNumNgb), -1.0/3));
	      else
		SphP[i].Hsml=  0.5*(PPP[i].Left + PPP[i].Right);

	      if(iter>25)
		printf("i=%d ID=%d task=%d left=%g right=%g ngb=%d\n", i, P[i].ID, ThisTask, PPP[i].Left, PPP[i].Right, SphP[i].NumNgb);

	    }
	  
	  density();
	  iter++;
	  if(iter>60)
	    {
	      endrun(1155);
	    }
	}
    }
  while(ntot>0);


  for(i=1;i<=N_gas;i++)
    {
      PPP[i].Left= PPP[i].Right= 0;

      if(SphP[i].Hsml < All.MinGasHsml)
	SphP[i].Hsml = All.MinGasHsml;
    }
}




#ifdef VELDISP 
/*
 *  This function is used to find an initial smoothing length for each 
 *  dark matter particle. It guarantees that the number of neighbours will be
 *  between desired_ngb-MAXDEV and desired_ngb+MAXDEV
 */
void setup_smoothinglengths_veldisp(int desired_ngb)
{
  int    i, ntot, last=0;
  float *r2list;
  int    *ngblist, count;
  int    iter=0;

  for(i=1+N_gas; i<=NumPart; i++) 
    { 
      P[i].HsmlVelDisp= sqrt(ngb_treefind( P[i].PosPred, desired_ngb, 0, P[i].Type, &ngblist, &r2list));  
      P[i].Left= P[i].Right= 0;
   } 

  veldisp();

  /* we will here store lower
   * and upper boundary of intervals for smoothing range bisection
   * in P[i].Left and P[i].Right
   */

  do
    {
      for(i=1+N_gas, NumForceUpdate=0, NumSphUpdate=0; i<=NumPart; i++) 
	{ 
	  if(P[i].NgbVelDisp < (All.DesNumNgb-All.MaxNumNgbDeviation ) 
	     || P[i].NgbVelDisp > (All.DesNumNgb+All.MaxNumNgbDeviation ))
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
		P[i].Left= P[i].HsmlVelDisp;

	      if(P[i].NgbVelDisp > (All.DesNumNgb+All.MaxNumNgbDeviation))
		P[i].Right= P[i].HsmlVelDisp;
	    }
	}

      MPI_Allreduce(&NumForceUpdate, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      
      if(ntot>0)
	{
	  if(ThisTask==0)
	    {
	      printf("\nNGB iteration %d.  still %d particles\n", iter, ntot);
	      fflush(stdout);
	    }     

	  for(i=IndFirstUpdate,count=0; count<NumForceUpdate; i=P[i].ForceFlag, count++)   
	    {
	      if(P[i].Left==0 || P[i].Right==0) 
		P[i].HsmlVelDisp=  P[i].HsmlVelDisp*( 0.5 + 0.5*pow(P[i].NgbVelDisp/((double)All.DesNumNgb), -1.0/3));
	      else
		P[i].HsmlVelDisp=  0.5*(P[i].Left + P[i].Right);

	      if(iter>25)
		printf("i=%d ID=%d task=%d left=%g right=%g ngb=%d\n", i, P[i].ID, ThisTask, P[i].Left, P[i].Right, P[i].NgbVelDisp);

	    }
	  
	  veldisp();
	  iter++;
	  if(iter>60)
	    {
	      endrun(1155);
	    }
	}
    }
  while(ntot>0);

}
#endif


#ifdef SIDM 
/*
 *  This function is used to find an initial smoothing length for each 
 *  dark matter particle. It guarantees that the number of neighbours will be
 *  between desired_ngb-MAXDEV and desired_ngb+MAXDEV
 */
void setup_smoothinglengths_sidm(int desired_ngb)
{
  int    i, ntot, last=0;
  float *r2list;
  int    *ngblist, count;
  int    iter=0;

  for(i=1+N_gas; i<=NumPart; i++) 
    { 
      P[i].HsmlVelDisp= sqrt(ngb_treefind( P[i].PosPred, desired_ngb, 0, P[i].Type, &ngblist, &r2list));  
      P[i].Left= P[i].Right= 0;
   } 

  setup_nbr_sidm(); /* SIDM */

  /* we will here store lower
   * and upper boundary of intervals for smoothing range bisection
   * in P[i].Left and P[i].Right
   */

  do
    {
      for(i=1+N_gas, NumForceUpdate=0, NumSphUpdate=0; i<=NumPart; i++) 
	{ 
	  if(P[i].NgbVelDisp < (All.DesNumNgb-All.MaxNumNgbDeviation ) 
	     || P[i].NgbVelDisp > (All.DesNumNgb+All.MaxNumNgbDeviation ))
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
		P[i].Left= P[i].HsmlVelDisp;

	      if(P[i].NgbVelDisp > (All.DesNumNgb+All.MaxNumNgbDeviation))
		P[i].Right= P[i].HsmlVelDisp;
	    }
	}

      MPI_Allreduce(&NumForceUpdate, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      
      if(ntot>0)
	{
	  if(ThisTask==0)
	    {
	      printf("\nNGB iteration %d.  still %d particles\n", iter, ntot);
	      fflush(stdout);
	    }     

	  for(i=IndFirstUpdate,count=0; count<NumForceUpdate; i=P[i].ForceFlag, count++)   
	    {
	      if(P[i].Left==0 || P[i].Right==0) 
		P[i].HsmlVelDisp=  P[i].HsmlVelDisp*( 0.5 + 0.5*pow(P[i].NgbVelDisp/((double)All.DesNumNgb), -1.0/3));
	      else
		P[i].HsmlVelDisp=  0.5*(P[i].Left + P[i].Right);

	      if(iter>25)
		printf("i=%d ID=%d task=%d left=%g right=%g ngb=%d\n", i, P[i].ID, ThisTask, P[i].Left, P[i].Right, P[i].NgbVelDisp);

	    }
	  
	  setup_nbr_sidm();
	  iter++;
	  if(iter>60)
	    {
	      endrun(1155);
	    }
	}
    }
  while(ntot>0);

}
#endif



