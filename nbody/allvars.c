/* This file declares all global variables. Further variables should be added here, 
   and declared as 'extern'. The actual existence of these variables is provided by
   the file 'allvars.c'. To produce 'allvars.c' from 'allvars.h', do the following:

   1.) Erase all #define's
   2.) add #include "allvars.h"
   3.) delete all keywords 'extern'
   4.) delete all struct definitions enclosed in {...}, e.g.
       "extern struct global_data_all_processes {....} All;"
       becomes "struct global_data_all_processes All;"
*/

#include "allvars.h"
#include <stdio.h>


#ifdef SIDM
#if (CROSS_SECTION_TYPE == 2 || CROSS_SECTION_TYPE == 4)
 double vc; /* Yukawa velocity scale */
#endif

 double vmax;
 long   iseed;
 int    n_scat_particles;
 int    scat_particles[MAX_SCAT];
#endif
 

 int    ThisTask;       /* the local processors  */
 int    NTask, PTask;   /* note: NTask = 2^PTask */

 double CPUThisRun;

 int    NumForceUpdate, IndFirstUpdate;
 int    NumSphUpdate;
 int    TimeTreeRoot;
 int    RestartFlag;
 int    NoCostFlag;
 int    Num_nodeupdates, Num_nodeupdate_particles;

 int    NumPart;        /* Note: this is the LOCAL processor value */
 int    N_gas;          /* Note: this is the LOCAL processor value */
 int    Ntype[6], NtypeLocal[6];

 float  DomainMin[6][3], DomainMax[6][3];  /* spatial extension of SPH particles locally */
 float  *InteriorMin, *InteriorMax;  /* ... and for all the other processors */



/* variables for input/output ,  usually only used on process 0 
 */
  char   ParameterFile[100];
  FILE  *FdInfo,
        *FdEnergy,
        *FdTimings,
        *FdCPU,
        *FdDEB;
#ifdef SFR
  FILE  *FdSfr;
#endif


/* tabulated SPH-smoothing kernel */

 double  Kernel[KERNEL_TABLE+2],
               KernelDer[KERNEL_TABLE+2],
               KernelRad[KERNEL_TABLE+2];

#ifdef SIDM
// double  CosTable[KERNEL_TABLE+2],
// 	     SinTable[KERNEL_TABLE+2];
#endif

 char *CommBuffer;  /* communication buffer, used at a number of places */



/* this structure contains data which is the SAME for all 
 * tasks (mostly code parameters read from the parameter file). 
 * Holding this data in a structure is convenient for writing/reading
 * the restart file, and it allows the introduction of new global
 * variables in a simple way. The only thing to do is to introduce them
 * into this structure.
 */
 struct global_data_all_processes  
 All;



/* The following structure holds all the information that is
 * stored for each particle of the simulation.
 */
 struct particle_data 
 *P,*P_data;



/* the following struture holds data that is stored for each SPH particle
 * in addition to the collisionless variables.
 */
 struct sph_particle_data
 *SphP,*SphP_data;


/* this structure holds nodes for the ordered binary tree of the timeline.
 */
 struct timetree_data 
 *PTimeTree;   



/* global state of system 
*/
 struct state_of_system 
 SysState, SysStateAtStart, SysStateAtEnd;






/* Various structure for communication during the gravity 
 * computation.
 */
 struct gravdata_in
 *GravDataIn;	 

 struct gravdata_out
 *GravDataResult, *GravDataPartialResult;	 

 double *GravDataPotential, *GravDataPartialPotential;




/* Various structure for communication during the density
 * computation.
 */
 struct densdata_in
 *DensDataIn;	 

 struct densdata_out
 *DensDataResult, *DensDataPartialResult;	 



/* Various structures for communication during the 
 * computation of hydrodynamical forces.
 */
 struct hydrodata_in
 *HydroDataIn;	 

 struct hydrodata_out
 *HydroDataResult, *HydroDataPartialResult;	 

#ifdef VELDISP
 struct veldispdata_in
 *VelDispDataIn;	 

 struct veldispdata_out
 *VelDispDataResult, *VelDispDataPartialResult;	 
#endif

#ifdef SIDM
 struct sidmdata_in
 *SidmDataIn;
 
 struct sidmdata_out
 *SidmDataResult, *SidmDataPartialResult;
 
 struct sidmdata_confirm
 *SidmDataConfirm;
 
 int *SidmTarget;
#endif


#ifdef SFR
 struct weightdata_in
 *WeightDataIn;	 

 struct weightdata_out
 *WeightDataResult, *WeightDataPartialResult;	 

 struct dissolvedata_in
 *DissolveDataIn;	 

 struct dissolvedata_out
 *DissolveDataResult, *DissolveDataPartialResult;
#endif


/* Structures for communication during the 
 * domain decomposition.
 */
  struct particle_data *DomainPartBuf;
  struct sph_particle_data *DomainSphBuf;




/* Header for the standard file format.
 */
 struct io_header_1
 header1;

// struct scatlog ScatLog;








