#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"



/*
 *  This function initializes the MPI communication packages,
 *  and sets cpu-time counters to 0.
 *  Then begrun() is called, which sets up the simulation
 *  either from IC's or from restart files.
 *  Finally, run() is started, the main simulation loop,
 *  which iterates over the timesteps.
 */
int main(int argc, char **argv)
{
  double t0, t1;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);

  // SIDM PARAMETER OUTPUT
  if(ThisTask == 0){
    fprintf(stdout, "GadgetSidm3 %s version\n", VERSION);
#if NOFORCE
    fprintf(stdout, "NOFORCE is ON.\n");
#endif
#if CLOUDCLOUD
    fprintf(stdout, "Use CLOUDCLOUD formula for scattering.\n");
#endif
  }

  if(NTask<=1)
    {
      if(ThisTask==0)
    fprintf(stdout, "Number of processors MUST be a larger than 1.\n");
      endrun(0);
    }

  for(PTask=0; NTask>(1<<PTask); PTask++);

  if(NTask!=(1<<PTask))
    {
      if(ThisTask==0)
    fprintf(stdout, "Number of processors MUST be a power of 2.\n");
      endrun(0);
    }

  if(argc<2)
    {
      if(ThisTask==0)
    {
      fprintf(stdout, "Parameters are missing.\n");
      fprintf(stdout, "Call with <ParameterFile> [<RestartFlag>]\n");
    }
      endrun(0);
    }

  strcpy(ParameterFile, argv[1]);

  if(argc>=3)
    RestartFlag= atoi(argv[2]);
  else
    RestartFlag= 0;

  
  All.CPU_TreeConstruction=All.CPU_TreeWalk=All.CPU_Gravity=All.CPU_Potential=All.CPU_Domain=
    All.CPU_Snapshot=All.CPU_Total=All.CPU_CommSum=All.CPU_Imbalance=All.CPU_Hydro=All.CPU_EnsureNgb=
       All.CPU_Predict= All.CPU_TimeLine= All.CPU_Diagnostic= 0;
#ifdef SFR
  All.CPU_Sfr=0;
#endif

  CPUThisRun=0;


  t0=second();

  begrun();       /* set-up run  */

  if(ThisTask == 0) {
#if (CROSS_SECTION_TYPE == 0)
    fprintf(stdout, "CROSS_SECTION_TYPE=0. Velocity Independent.\n");  
#elif (CROSS_SECTION_TYPE == 1)
    fprintf(stdout, "CROSS_SECTION_TYPE=1. Maxwellian Molecule.\n");
#elif (CROSS_SECTION_TYPE == 2)
    fprintf(stdout, "CROSS_SECTION_TYPE=2. Yukawa-like Interaction.\n");
#elif (CROSS_SECTION_TYPE == 3)
    fprintf(stdout, "CROSS_SECTION_TYPE=3. Power-Law n=%f, vel=%f.\n",
	    All.CrossSectionPowLaw, All.CrossSectionVelScale);
#elif (CROSS_SECTION_TYPE == 4)
    fprintf(stdout, "CROSS_SECTION_TYPE=4. Yukawa Interaction.\n");
#endif
  }

  t1=second();
  CPUThisRun+= timediff(t0,t1);
  All.CPU_Total+= timediff(t0,t1);



  run();          /* main simulation loop */

  MPI_Finalize(); /* clean up & finalize MPI */

  return 0;
}







