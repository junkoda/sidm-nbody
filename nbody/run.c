#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"


/* This routine contains the main simulation loop that iterates over
 * the single timesteps. The loop terminates when the cpu-time
 * limit is reached, when a `stop' file is found in the output 
 * directory, or when the simulation ends because we arrived
 * at TimeMax
 */
void run(void)    
{
  FILE   *fd;
  int    i, stopflag=0;
  double savetime;
  char   buf[200],stopfname[200];
  double t0,t1;
#ifdef SFR
  double t0_sfr,t1_sfr;
  int    starscreated;
#endif
#ifdef SIDM
  double tstart, tend;
#endif

  sprintf(stopfname,"%sstop",All.OutputDir);

  do   /* main loop */
    {
      t0=second();

      find_next_time();   /* increases Time to smallest max-prediction time and
                             determines which particles are grouped together
                             for force evaluation */

      if(All.Time > All.TimeMax)
	All.Time = All.TimeMax;
      
#ifdef COOLING
      IonizeParams();
#endif
      
      every_timestep_stuff();   /* write some info to log-files */

      if((All.Time-All.TimeLastStatistics)>=All.TimeBetStatistics) /* check whether we want a full statistics */
	{
	  savetime=All.Time;
	  predict(All.Time=All.TimeLastStatistics+All.TimeBetStatistics);
	  compute_potential();
	  energy_statistics(); /* compute and output energy statistics */
	  All.TimeLastStatistics += All.TimeBetStatistics;
	  All.Time=savetime;
	}

      if((All.Time-All.TimeOfFirstSnapshot)>=0) /* check whether it's time for a snapshot file */
	{ 
	  savetime=All.Time;
	  predict(All.Time=All.TimeOfFirstSnapshot);
	  savepositions(All.SnapshotFileCount++);   /* write snapshot file */
	  if(All.OutputListOn)
	    All.TimeOfFirstSnapshot= find_next_outputtime(savetime);
	  else
	    if(All.ComovingIntegrationOn)
	      All.TimeOfFirstSnapshot *= All.TimeBetSnapshot;
	    else
	      All.TimeOfFirstSnapshot += All.TimeBetSnapshot;
	  All.Time=savetime;
	}
	
      predict_sph_particles(All.Time);  /* SPH particles are alwyas predicted, while
                                           this is done for collisionless particles
                                           either before the tree construction,
                                           or on the fly while the tree is walked */

      compute_accelerations(0);
      /* compute accelerations for 
	 the particles that are to be advanced 
	 Note: the particle positions and tree nodes are 
	 predicted for the current time during the force 
	 computation */ 

      advance();                  /* and advance the active particles  */
      /*
#ifdef SIDM
      tstart=second();
      determine_interior();
  
      //if(mode == 0){
      sidm();
      sidm_ensure_neighbours(0);
      //}
      
      tend=second();
      All.CPU_EnsureNgb+= timediff(tstart, tend);
#endif
      */


#ifdef REFLECTIONBOUNDARY
      reflect();
#endif

#ifdef SIDM
      // update force tree that is affected by the scattering
      update_node_sidm();
#endif

      find_timesteps(0);           /* compute new timesteps for these particles, 
                                     and set their maximum prediction times, 
                                     and reinsert them into the timeline as needed */

      /* Check whether it is time for a new domain decomposition*/
      if(All.NumForcesSinceLastDomainDecomp > All.TotNumPart*All.DomainUpdateFrequency) 
	{
#ifdef SIDM
	  vmax= getvmax();
#endif
	  force_costevaluate();   /* assign accumulated cost to particles */
#ifdef SFR 
	  t0_sfr=second();

	  starscreated= dissolvegas();             /* check whether gas particles are to be converted
 					      fully into collisionless stars */
	  t1_sfr=second();
	  All.CPU_Sfr+= t1_sfr-t0_sfr;
#endif

#ifdef PERIODIC
	  do_box_wrapping();      /* map the particles back onto the box */
#endif	  
	  DomainDecomposition();  /* do domain decomposition */

	  for(i=1;i<=NumPart;i++)   
	    P[i].GravCost*= 0.5;  /* reset cost of particles to half their original value */

	  construct_timetree();   /* construct new timeline */

	  All.NumForcesSinceLastDomainDecomp=0;
	  All.NumForcesSinceLastTreeConstruction= All.TreeUpdateFrequency*All.TotNumPart ; /* ensures that new tree will be constructed */
	  NoCostFlag=1;  /* flags that force_costevaluate() needs not to be called before next tree construction */
	}


      All.NumCurrentTiStep++;

      if(ThisTask==0)  /* Check whether we need to interrupt the run */
	{
	  if((fd=fopen(stopfname,"r")))   /* Is the stop-file present? If yes, interrupt run. */
	    {
	      fclose(fd);
	      stopflag=1;
	      sprintf(buf,"rm -f %s",stopfname);
	      //system(buf);
	    }

	  if(CPUThisRun > 0.85*All.TimeLimitCPU)   /* are we running out of CPU-time ? If yes, interrupt run. */ 
	    {
	      printf("reaching time-limit. stopping.\n");
	      stopflag=2;
	    }
	}

      MPI_Bcast(&stopflag, 1, MPI_INT, 0, MPI_COMM_WORLD);

      if(stopflag)
	{
	  restart(0);  /* write restart file */
	  MPI_Barrier(MPI_COMM_WORLD);

	  if(stopflag==2 && All.ResubmitOn && ThisTask==0)
	    {
	      close_outputfiles();
	      sprintf(buf,"%s", All.ResubmitCommand);
	      //system(buf); 
	    }
	  return;
	}

      if(ThisTask==0) /* is it time to write a restart-file? (for security) */
	{
	  if((CPUThisRun - All.TimeLastRestartFile)>= All.CpuTimeBetRestartFile) 
	    {
	      All.TimeLastRestartFile= CPUThisRun;
	      stopflag= 3;
	    }
	  else 
	    stopflag= 0;
	}

      MPI_Bcast(&stopflag, 1, MPI_INT, 0, MPI_COMM_WORLD);
      
      if(stopflag==3)
	{
	  restart(0);                        /* write an occasional restart file */
	  stopflag=0;
	}

      t1=second();

      All.CPU_Total+= timediff(t0,t1);
      CPUThisRun+=    timediff(t0,t1);
    }
  while(All.Time < All.TimeMax);

  restart(0);  /* write restart file to allow a continuation of the run beyond TimeMax */


  /* write a last snapshot file at final time (will be overwritten if 
     All.TimeMax is increased and the run is continued) */

  savetime=All.Time;
  predict(All.Time=All.TimeMax);
  savepositions(All.SnapshotFileCount++);   /* write snapshot file */
  All.Time=savetime;
}



/* This routine writes one line for every timestep to two log-files.
 * In FdInfo, we just list the timesteps that have been done,
 * while in FdCPU the cumulative cpu-time consumption in various parts
 * of the code is stored.
 */
void every_timestep_stuff(void)
{
  double z;

  if(ThisTask==0)
    {
      if(All.ComovingIntegrationOn)
	{
	  z=1.0/(All.Time)-1;
	  fprintf(FdInfo,"\nBegin Timestep %d, Time: %g, Redshift: %g, Timestep: %g\n",All.NumCurrentTiStep,All.Time,z,All.TimeStep);
	  printf(        "\nBegin Timestep %d, Time: %g, Redshift: %g, Timestep: %g\n",All.NumCurrentTiStep,All.Time,z,All.TimeStep);
	  fflush(FdInfo);
	}
      else
	{
	  fprintf(FdInfo,"\nBegin Timestep %d, Time: %g, Timestep: %g\n",All.NumCurrentTiStep,All.Time,All.TimeStep);
	  printf(        "\nBegin Timestep %d, Time: %g, Timestep: %g\n",All.NumCurrentTiStep,All.Time,All.TimeStep);
	  fflush(FdInfo);
	}
    
      fprintf(FdCPU,"Timestep %d, Time: %g\n",All.NumCurrentTiStep,All.Time);
#ifdef SFR
      fprintf(FdCPU,"%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n",
	      All.CPU_Total,
	      All.CPU_Gravity,
	      All.CPU_Hydro,
	      All.CPU_Domain,
	      All.CPU_Potential,
              All.CPU_Predict,
              All.CPU_TimeLine,
	      All.CPU_Snapshot,
	      All.CPU_TreeWalk,
	      All.CPU_TreeConstruction,
	      All.CPU_CommSum,
	      All.CPU_Imbalance,
	      All.CPU_EnsureNgb,
              All.CPU_Diagnostic,
              All.CPU_Sfr);
#else
      fprintf(FdCPU,"%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n",
	      All.CPU_Total,
	      All.CPU_Gravity,
	      All.CPU_Hydro,
	      All.CPU_Domain,
	      All.CPU_Potential,
              All.CPU_Predict,
              All.CPU_TimeLine,
	      All.CPU_Snapshot,
	      All.CPU_TreeWalk,
	      All.CPU_TreeConstruction,
	      All.CPU_CommSum,
	      All.CPU_Imbalance,
	      All.CPU_EnsureNgb,
              All.CPU_Diagnostic);
#endif
      fflush(FdCPU);
    }
}


/* This routine first calls a computation of various global quantities
 * of the particle distribution, and then writes some statistics
 * about the energies in the various particle components to the 
 * file FdEnergy.
 */
void energy_statistics(void)
{
  compute_global_quantities_of_system();

  if(ThisTask==0)
    {
      fprintf(FdEnergy,"%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	      All.Time,
	      SysState.EnergyInt,
	      SysState.EnergyPot,
	      SysState.EnergyKin,
	      SysState.EnergyIntComp[0],
	      SysState.EnergyPotComp[0],
	      SysState.EnergyKinComp[0],
	      SysState.EnergyIntComp[1],
	      SysState.EnergyPotComp[1],
	      SysState.EnergyKinComp[1],
	      SysState.EnergyIntComp[2],
	      SysState.EnergyPotComp[2],
	      SysState.EnergyKinComp[2],
	      SysState.EnergyIntComp[3],
	      SysState.EnergyPotComp[3],
	      SysState.EnergyKinComp[3],
	      SysState.EnergyIntComp[4],
	      SysState.EnergyPotComp[4],
	      SysState.EnergyKinComp[4],
	      SysState.MassComp[0],
	      SysState.MassComp[1],
	      SysState.MassComp[2],
	      SysState.MassComp[3],
	      SysState.MassComp[4]
	      ); 

      fflush(FdEnergy);
    }
}
































