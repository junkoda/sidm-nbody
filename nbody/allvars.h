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


#include <stdio.h>
#include "tags.h"


#ifdef T3E
  typedef short int int4byte;   /* Note: int has 8 Bytes on the T3E ! */
#else
  typedef int int4byte;
#endif


/* SIDM DEBUG OPTIONS */
#define  VERSION "Feb 3, 2008"

#ifndef NOFORCE
#define  NOFORCE 0
#endif
  /* i.e. NO VELOCITY UPDATE */

#define  CLOUDCLOUD 0
  /* i.e USE overlapping length */

// Moved Makefile OPT
//#define  REFLECTIONBOUNDARY 0
//#define  SCATTERLOG         1
//#define  FINDNBRLOG         1
  //#define  RANDOM_GSL           // don't define when you don't use GSL

#define  SCATKERNELFACTOR   1.0
#define  SAFEFACTOR         1.0

//#ifndef CROSS_SECTION_TYPE
//#define  CROSS_SECTION_TYPE 0 // hard sphere 0, Maxwellian 1
// #endif

/* ... often used constants (cgs units) */

#define  MAX_REAL_NUMBER  1e37
#define  MIN_REAL_NUMBER  1e-37

#define  THIRD            (1.0/3.0)
#define  PI               3.14159265358979323846 
#define  PI_INV           (1/PI)
#define  LN2              0.69314718

#define  GRAVITY     6.672e-8
#define  SOLAR_MASS  1.989e33
#define  SOLAR_LUM   3.826e33
#define  RAD_CONST   7.565e-15
#define  AVOGADRO    6.0222e23
#define  BOLTZMANN   1.3806e-16
#define  GAS_CONST   8.31425e7
#define  C           2.9979e10
#define  PLANCK      6.6262e-27
#define  CM_PER_MPC  3.085678e24
#define  PROTONMASS  1.6726e-24

#define  HUBBLE      3.2407789e-18   /* in h/sec */

#define  SEC_PER_MEGAYEAR   3.155e13
#define  SEC_PER_YEAR       3.155e7

#define  GAMMA         (5.0/3)
#define  GAMMA_MINUS1  (GAMMA-1)
 
#define KERNEL_TABLE 1000

#define MAX_NGB  5000000  /* defines maximum length of neighbour list */

#define  MAXLEN_OUTPUTLIST 350 /* maxmimum number of entries in output list*/

#define  TIMESTEP_INCREASE_FACTOR 1.3


#define  HYDROGEN_MASSFRAC 0.76

#define  DISSOLVE_FRAC     0.10 
#define  STARFORM_FRAC     0.333 

#ifdef SIDM

#define  BALLINVERSE       (3./4./PI)
#define  OVERLAP           0
#define  MAX_SCAT           10000
 /* number of scattered particles per timestep.
    Used for scat_particle list only */

#if (CROSS_SECTION_TYPE == 2 || CROSS_SECTION_TYPE == 4)
extern double vc;
#endif
extern double vmax;
extern long   iseed;
extern int    n_scat_particles;
extern int4byte scat_particles[MAX_SCAT];

#endif

extern int    ThisTask;       /* the local processors  */
extern int    NTask, PTask;   /* note: NTask = 2^PTask */

extern double CPUThisRun;

extern int    NumForceUpdate, IndFirstUpdate;
extern int    NumSphUpdate;
extern int    TimeTreeRoot;
extern int    RestartFlag;
extern int    NoCostFlag;
extern int    Num_nodeupdates, Num_nodeupdate_particles;

extern int    NumPart;        /* Note: this is the LOCAL processor value */
extern int    N_gas;          /* Note: this is the LOCAL processor value */
extern int    Ntype[6], NtypeLocal[6];


extern float  DomainMin[6][3], DomainMax[6][3];  /* spatial extension of SPH particles locally */
extern float  *InteriorMin, *InteriorMax;  /* ... and for all the other processors */



/* variables for input/output ,  usually only used on process 0 
 */
extern  char   ParameterFile[100];
extern  FILE  *FdInfo,
              *FdEnergy,
              *FdTimings,
              *FdCPU,
              *FdDEB;

#ifdef SFR
extern  FILE  *FdSfr;
#endif


/* tabulated SPH-smoothing kernel */

extern double  Kernel[KERNEL_TABLE+2],
               KernelDer[KERNEL_TABLE+2],
               KernelRad[KERNEL_TABLE+2];

#ifdef SIDM
//extern double  CosTable[KERNEL_TABLE+2],
//               SinTable[KERNEL_TABLE+2];

#endif

extern char *CommBuffer;  /* communication buffer, used at a number of places */



/* this structure contains data which is the SAME for all 
 * tasks (mostly code parameters read from the parameter file). 
 * Holding this data in a structure is convenient for writing/reading
 * the restart file, and it allows the introduction of new global
 * variables in a simple way. The only thing to do is to introduce them
 * into this structure.
 */
extern struct global_data_all_processes  
{
  int    TotNumPart,    /*  particle numbers (total)  Note: these are global values */
         TotN_gas,
         TotN_halo,
         TotN_disk,
         TotN_bulge,
         TotN_stars;

  int    MaxPart;        /* This gives the maxmimum number of particles that can be 
                stored on one processor. */
  int    MaxPartSph;     /* same for SPH */


  double MinGasMass;     /* This block contains variables needed in the  */
  double EgySpecCold;    /* star formation and feedback sector */
  double EgySpecSN;
  double OverDensThresh;
  double PhysDensThresh;
  double FeedbackEnergy;
  double TempSupernova;
  double TempClouds;
  double CritOverDensity;
  double CritPhysDensity;
  double FactorSFR;
  double FactorSN;
  double FactorEVP;


  int    ICFormat;       /* selects different versions of IC file-format */

  int    NumFilesPerSnapshot; 
  int    NumFilesWrittenInParallel; 

  int    BufferSize;         /* size of communication buffer in MB */
  int    BunchSizeForce;     /* number of particles fitting into the buffer */
  int    BunchSizeDensity;
  int    BunchSizeHydro;
  int    BunchSizeDomain;
#ifdef VELDISP
  int    BunchSizeVelDisp;
#endif
#ifdef SIDM
  int    BunchSizeSidm;
#endif
#ifdef SFR
  int    BunchSizeWeight;
  int    BunchSizeDissolve;
#endif

  double PartAllocFactor;  /* in order to maintain work-load balance,  
                  the particle load will usually NOT be balanced.
                  Each processor allocates memory for PartAllocFactor times  
                  the average number of particles to allow for that */
                           
  double TreeAllocFactor;  /* similarly for the tree:
                  each processor allocates a number of nodes which is 
                  TreeAllocFactor times the maximum(!) number of particles. 
                  Note: A typical local tree for N particles needs 
                  usually ~0.65*N nodes */

  /* some SPH parameters */

  int    DesNumNgb;
  int    MaxNumNgbDeviation;
  double ArtBulkViscConst;
  double InitGasTemp;     /* may be used to set the temperature in the IC's */   
  double MinGasTemp;      /* may be used to set a floor for the gas temperature */
  double MinEgySpec;


  /* some force counters  */

  int    TotNumOfForces;   /* counts total number of force computations  */

  int    NumForcesSinceLastDomainDecomp;
  int    NumForcesSinceLastTreeConstruction;

  /* system of units  */

  double UnitTime_in_s,
         UnitMass_in_g,
         UnitVelocity_in_cm_per_s,
         UnitLength_in_cm,
         UnitPressure_in_cgs,
         UnitDensity_in_cgs,
         UnitCoolingRate_in_cgs,
         UnitEnergy_in_cgs,
         UnitTime_in_Megayears,
         GravityConstantInternal,
         G;

  /* Cosmology */

  double Hubble;
  double BoxSize, BoxHalf;
  double Omega0,        
         OmegaLambda,
         OmegaBaryon,
         HubbleParam;   /* little `h', i.e. Hubble constant in units of 100 km/s/Mpc. 
                         * Only needed to get absolute physical values 
             * for cooling physics 
             */    

  /* Code options */

  int    ComovingIntegrationOn;   /* enables comoving integration */
  int    PeriodicBoundariesOn;
  int    ResubmitOn;
  int    TypeOfOpeningCriterion;
  int    TypeOfTimestepCriterion;
  int    OutputListOn;
  int    CoolingOn;
  int    StarformationOn;
  int    MultiPhaseModelOn;


  /* parameters determining output frequency */

  int    SnapshotFileCount;
  double TimeBetSnapshot,
         TimeOfFirstSnapshot,
         CpuTimeBetRestartFile,
         TimeLastRestartFile,
         TimeBetStatistics,
         TimeLastStatistics;


  /* Current time of the simulation, global step, and end of simulation */

  int     NumCurrentTiStep;

  double  Time, 
          TimeBegin,
          TimeStep,
          TimeMax;   /* marks end of the simulation */


 /* variables that keep track of cumulative CPU consumption */

  double  TimeLimitCPU;
  double  CPU_TreeConstruction;
  double  CPU_TreeWalk;
  double  CPU_Gravity;
  double  CPU_Potential;
  double  CPU_Domain;
  double  CPU_Snapshot;
  double  CPU_Total;
  double  CPU_CommSum;
  double  CPU_Imbalance;
  double  CPU_Hydro;
  double  CPU_EnsureNgb;
  double  CPU_Predict;
  double  CPU_TimeLine;
  double  CPU_Diagnostic;
#ifdef SFR
  double  CPU_Sfr;
#endif

  /* tree code opening criterion */

  double  ErrTolTheta;        /* BH-opening angle */
  double  ErrTolForceAcc;     /* for new opening criterion */

  /* adjusts accuracy of time-integration */

  double  ErrTolIntAccuracy;  /* for 1/a^{1/2} collisionless timestep criterion */
  double  ErrTolDynamicalAccuracy; /* March 24 */
  double  ErrTolVelScale;     /* for 1/a  collisionless timestep criterion */
  double  MinSizeTimestep,
          MaxSizeTimestep;

  double  CourantFac;      /* SPH-Courant factor */
 

  /* frequency of tree reconstruction/domain decomposition */

  double  MaxNodeMove;
  double  TreeUpdateFrequency;
  double  DomainUpdateFrequency;


  /* gravitational and hydrodynamical softening lengths 
   * (given in terms of an `equivalent' Plummer softening length) 
   *
   * five groups of particles are supported 
   * 0=gas,1=halo,2=disk,3=bulge,4=stars 
   */
  double  MinGasHsmlFractional, MinGasHsml;

  double  SofteningGas,
          SofteningHalo,
          SofteningDisk,
          SofteningBulge,
          SofteningStars;

  double  SofteningGasMaxPhys,
          SofteningHaloMaxPhys,
          SofteningDiskMaxPhys,
          SofteningBulgeMaxPhys,
          SofteningStarsMaxPhys;

  double  SofteningTable[6];    
  double  SofteningTableMaxPhys[6];


  /* If particle masses are all equal for one type, 
   *  the corresponding entry in MassTable is set to this value,
   * allowing the size of the snapshot files to be reduced
   */
  double  MassTable[6];   
  

  /* some filenames */
  char    InitCondFile[100],
          OutputDir[100],
          SnapshotFileBase[100],
          EnergyFile[100],
          CpuFile[100],
          InfoFile[100],
          TimingsFile[100],
          RestartFile[100],
          ResubmitCommand[100],
          OutputListFilename[100];

  double  OutputListTimes[MAXLEN_OUTPUTLIST]; /* was 200 in earlier version */
  int     OutputListLength;
  
#ifdef SIDM
  double  CrossSection;            /* hinverse cm2/g */
      /* CROSS_SECTION_TYPE=1 => CrossSection=sigma*(v km/s) hinv cm2/g */
  double  CrossSectionInternal;
#if (CROSS_SECTION_TYPE == 2 || CROSS_SECTION_TYPE == 4)
  double  YukawaVelocity;
#elif (CROSS_SECTION_TYPE == 3)
  double CrossSectionPowLaw, CrossSectionVelScale;
#endif
  int     Seed1, Seed2;
  double  ProbabilityTol;
#endif

#ifdef REFLECTIONBOUNDARY
  double  ReflectionRadius;
#endif

} All;



/* The following structure holds all the information that is
 * stored for each particle of the simulation.
 */
extern struct particle_data 
{
  float     Pos[3];       /* particle position at its current time */  
  float     Vel[3];       /* particle velocity at its current time */  
  float     Mass;         /* particle mass */   
  int4byte  ID;           /* unique particle identifier */  
  int4byte  Type;         /* flags particle type. 0=gas, 1=halo, 2=disk, 3=bulge, 4=stars */

  float     CurrentTime;  /* current time of the particle */
  float     MaxPredTime;  /* current time plus half the particles allowed timestep */
  float     PosPred[3];   /* particle position at the global prediction time */   
  float     VelPred[3];   /* particle velocity at the global prediction time */

  float     Accel[3];     /* particle acceleration */
  float     Potential;    /* particle potential */

  float     GravCost;     /* weight factor used to balance the work-load */
  float     OldAcc;       /* magnitude of old force. Used in new relative opening criterion */

  int4byte  ForceFlag;    /* points to next active particle */

#ifdef VELDISP  
  float     Left, Right;
  float     VelDisp;
  int4byte  NgbVelDisp;
  float     HsmlVelDisp;
  float     DensVelDisp;
#endif
#ifdef SIDM
  float     Left, Right;
  int4byte  NgbVelDisp;
  float     HsmlVelDisp;
  float     dVel[3];
  // int4byte  ForceFlagBackup;
#endif
#ifdef STELLARAGE
  float     MeanStellarAge;
#endif
} *P,*P_data;



/* the following struture holds data that is stored for each SPH particle
 * in addition to the collisionless variables.
 */
extern struct sph_particle_data
{
  float  Density;         /* particle density at its current time */ 
  float  DtDensity;       /* rate of change of density */
  float  DensityPred;     /* predicted particle density */

  float  EgySpec;         /* internal energy per unit mass */
  float  DtEgySpec;       /* rate of change of the internal energy */
  float  EgySpecPred;     /* predicted internal energy per unit mass */

  float  Pressure;        /* pressure */

  float  Hsml;            /* smoothing length */
  float  DtHsml;          /* rate of change of smoothing length */

  int    NumNgb;          /* number of SPH neighbours */

  float  DivVel;          /* local velocity divergence */
  float  CurlVel;         /* local velocity curl */

#ifndef VELDISP  
  float    Left, Right;
#endif

#ifdef COOLING
  float  Ne;              /* electron fraction. Gives indirectly ionization state 
                             and mean molecular weight. */
#endif

#ifdef SFR
  float FormedStellarMass;
#ifdef CLOUDS
  float CloudMass;
#endif
#endif

} *SphP,*SphP_data;


/* this structure holds nodes for the ordered binary tree of the timeline.
 */
extern struct timetree_data 
{
  int4byte left, right;
} *PTimeTree;   



/* global state of system 
*/
extern struct state_of_system 
{
  double  Mass,
          EnergyKin,
          EnergyPot,
          EnergyInt,
          EnergyTot,
          Momentum[4],
          AngMomentum[4],
          CenterOfMass[4],

          MassComp[5],
          EnergyKinComp[5],
          EnergyPotComp[5],
          EnergyIntComp[5],
          EnergyTotComp[5],
          MomentumComp[5][4],
          AngMomentumComp[5][4],
          CenterOfMassComp[5][4];

} SysState, SysStateAtStart, SysStateAtEnd;






/* Various structure for communication during the gravity 
 * computation.
 */
extern struct gravdata_in
{
  float Pos[3];
  float OldAcc;
  int   Type;

} *GravDataIn;   

extern struct gravdata_out
{
  double Acc[3];

} *GravDataResult, *GravDataPartialResult;   

extern double *GravDataPotential, *GravDataPartialPotential;




/* Various structure for communication during the density
 * computation.
 */
extern struct densdata_in
{
  float Pos[3];
  float Vel[3];
  float Hsml;

} *DensDataIn;   

extern struct densdata_out
{
  float Rho;
  int   Ngb;
  float Div, Rot[3];

} *DensDataResult, *DensDataPartialResult;   


#ifdef VELDISP
extern struct veldispdata_in
{
  float Pos[3];
  float Hsml;
  int   Type;
  
} *VelDispDataIn;    

extern struct veldispdata_out
{
  int   Ngb;
  float Rho;
  float Vsum[3];
  float V2sum[3];

} *VelDispDataResult, *VelDispDataPartialResult;     
#endif

#ifdef SIDM
extern struct sidmdata_in
{
  float Pos[3];
  float Vel[3];
  float Mass;
  float dt;
  int4byte ID;
  float Hsml;
  int   Type;
  
} *SidmDataIn;   

extern struct sidmdata_out
{
  int   Ngb;
  float dv[3];

} *SidmDataResult, *SidmDataPartialResult;

extern struct sidmdata_confirm
{
  int Task;
} *SidmDataConfirm;

extern int4byte *SidmTarget;

void   initsidm(int seed);
double getvmax(void);
void   sidm(void);
void   sidm_ensure_neighbours(int mode);
void   setup_smoothinglengths_sidm(int desired_ngb);
void   setup_nbr_sidm(void);
//double ran2(long *idum);
void   update_node_of_scat_particle(int i);
void   update_node_sidm(void);

#endif


#ifdef SFR
/* Various structure for communication during the actual spawning
 * of star particles
 */
extern struct weightdata_in
{
  float Pos[3];
  float Hsml;
} *WeightDataIn;     

extern struct weightdata_out
{
  double WeightSum;
} *WeightDataResult, *WeightDataPartialResult;   

/* Various structures for dissolving gas particles.
 */
extern struct dissolvedata_in
{
  float Pos[3];
  float Vel[3];
  float Hsml;
  float Mass; 
  double WeightSum;
  float FormedStellarMass;
#ifdef CLOUDS
  float CloudMass;
#endif
  float U;
} *DissolveDataIn;   

extern struct dissolvedata_out
{
  float Pos[3];
  float Vel[3];
  float Mass;
} *DissolveDataResult, *DissolveDataPartialResult;

#endif



/* Various structures for communication during the 
 * computation of hydrodynamical forces.
 */
extern struct hydrodata_in
{
  float Pos[3];
  float Vel[3];
  float Hsml;
  float Mass;
  float DensityPred;
  float Pressure;
  float Current;
  float F1;
#ifdef SFR
  float FormedStellarMass;
#ifdef CLOUDS
  float CloudMass;
#endif
#endif
} *HydroDataIn;  

extern struct hydrodata_out
{
  float Acc[3];
  float DtEgySpec;
} *HydroDataResult, *HydroDataPartialResult;     



/* Structures for communication during the 
 * domain decomposition.
 */
extern  struct particle_data *DomainPartBuf;
extern  struct sph_particle_data *DomainSphBuf;




/* Header for the standard file format.
 */
extern struct io_header_1
{
  int4byte npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int4byte flag_sfr;
  int4byte flag_feedback;
  int4byte npartTotal[6];
  int4byte flag_cooling;
  int4byte num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam;
  int4byte flag_multiphase;
  int4byte flag_stellarage;
  int4byte flag_sfrhistogram;
  char     fill[84];  /* fills to 256 Bytes */
} header1;

/*
extern struct scatlog
{
  float    time;
  int4byte id1;
  int4byte id2;
  float    Hsml1, Hsml2;
  float    x1[3], x2[3];
  float    v1[3], v2[3];
  float    dv[3];
} ScatLog;
*/
