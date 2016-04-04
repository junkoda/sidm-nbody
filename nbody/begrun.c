#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"
#include "sidm_rand.h"

/*
 *  Does the initial set-up of the simulation.
 *  Reading the parameterfile, setting units, 
 *  getting IC's/restart files, etc.
 */
void begrun(void)
{
  struct global_data_all_processes all;

  read_parameter_file(ParameterFile);   /* ... read in parameters for this run */

  set_sph_kernel();                    

  allocate_commbuffers();               /* ... allocate buffer-memory for particle 
                                           exchange during force computation */

  set_units();

#ifdef COOLING
  All.Time = All.TimeBegin;   
  InitCool();
#endif

#ifdef PERIODIC
  All.BoxHalf= All.BoxSize/2;
  ewald_init();
#endif
  
  open_outputfiles(); 

  srand48(42);                          /* put a seed for the random number generator */

#ifdef SIDM
  init_rand(All.Seed1 + All.Seed2*ThisTask, ThisTask);
#endif

  All.TimeLastRestartFile= CPUThisRun;

  if(RestartFlag==0 || RestartFlag==2)
    {
      init();    /* ... read in initial model */
    }
  else
    {
      all=All;  /* save global variables. (will be read from restart file) */

      restart(RestartFlag);  /* ... read restart file. Note: This also resets 
			        all variables in the struct `All'. 
				However, during the run, some variables in the parameter
				file are allowed to be changed, if desired. These need to 
                                copied in the way below.
                                Note:  All.PartAllocFactor is treated in restart() separately.  
				*/

      All.TimeMax=                  all.TimeMax;
      All.MinSizeTimestep=          all.MinSizeTimestep;
      All.MaxSizeTimestep=          all.MaxSizeTimestep;
      All.TreeAllocFactor=          all.TreeAllocFactor;
      All.BufferSize=               all.BufferSize;
      All.BunchSizeForce=           all.BunchSizeForce;
      All.BunchSizeDensity=         all.BunchSizeDensity;
      All.BunchSizeHydro=           all.BunchSizeHydro;
      All.BunchSizeDomain=          all.BunchSizeDomain;
      All.TimeLimitCPU=             all.TimeLimitCPU;
      All.ResubmitOn=               all.ResubmitOn;
      All.TimeBetSnapshot=          all.TimeBetSnapshot;
      All.TimeBetStatistics=        all.TimeBetStatistics;
      All.CpuTimeBetRestartFile=    all.CpuTimeBetRestartFile;
      All.ErrTolIntAccuracy=        all.ErrTolIntAccuracy;
      All.ErrTolDynamicalAccuracy=  all.ErrTolDynamicalAccuracy;      
      All.ErrTolVelScale=           all.ErrTolVelScale;
      All.ErrTolTheta=              all.ErrTolTheta;
      All.ErrTolForceAcc=           all.ErrTolForceAcc;
      All.TypeOfTimestepCriterion=  all.TypeOfTimestepCriterion;
      All.TypeOfOpeningCriterion=   all.TypeOfOpeningCriterion;
      All.NumFilesWrittenInParallel=all.NumFilesWrittenInParallel;
      All.DomainUpdateFrequency=    all.DomainUpdateFrequency;
      All.TreeUpdateFrequency=      all.TreeUpdateFrequency;
      All.MaxNodeMove=              all.MaxNodeMove;
      All.OutputListOn=             all.OutputListOn;

      All.OutputListLength=       all.OutputListLength;
      memcpy(All.OutputListTimes, all.OutputListTimes, sizeof(double)*All.OutputListLength);

      strcpy(All.ResubmitCommand, all.ResubmitCommand);
      strcpy(All.OutputListFilename, all.OutputListFilename);
      strcpy(All.OutputDir,   all.OutputDir);
      strcpy(All.RestartFile, all.RestartFile);
      strcpy(All.EnergyFile,  all.EnergyFile);
      strcpy(All.InfoFile,    all.InfoFile);
      strcpy(All.CpuFile,     all.CpuFile);
      strcpy(All.TimingsFile, all.TimingsFile);
      strcpy(All.SnapshotFileBase, all.SnapshotFileBase);

#ifdef SIDM
      vmax= getvmax();
#endif
      force_treeallocate(All.TreeAllocFactor*All.MaxPart, All.MaxPart);

      /* ensures that new tree will be constructed */
      All.NumForcesSinceLastTreeConstruction= All.TreeUpdateFrequency*All.TotNumPart ;
      NoCostFlag=1;

      ngb_treeallocate(MAX_NGB);

      construct_timetree();  /* reconstruct timeline */
    }
  
  if(All.OutputListOn)
    All.TimeOfFirstSnapshot= find_next_outputtime(All.Time);

  All.TimeLastRestartFile= CPUThisRun;
}




/*
 *  Compute conversion factors between internal code units
 *  and the cgs-system.  
 */
void set_units(void)   
{
#ifdef SFR
  double feedbackenergyinergs;
#endif

  All.UnitTime_in_s= All.UnitLength_in_cm / All.UnitVelocity_in_cm_per_s;
  All.UnitTime_in_Megayears= All.UnitTime_in_s/SEC_PER_MEGAYEAR;

  if(All.GravityConstantInternal==0)
    All.G= GRAVITY/pow(All.UnitLength_in_cm,3)*All.UnitMass_in_g*pow(All.UnitTime_in_s,2);
  else
    All.G= All.GravityConstantInternal;

  All.UnitDensity_in_cgs= All.UnitMass_in_g/pow(All.UnitLength_in_cm,3);
  All.UnitPressure_in_cgs= All.UnitMass_in_g/All.UnitLength_in_cm/pow(All.UnitTime_in_s,2);
  All.UnitCoolingRate_in_cgs= All.UnitPressure_in_cgs/All.UnitTime_in_s;
  All.UnitEnergy_in_cgs= All.UnitMass_in_g * pow(All.UnitLength_in_cm,2) / pow(All.UnitTime_in_s,2);

  /* convert some physical input parameters to internal units */

  All.Hubble = HUBBLE * All.UnitTime_in_s;

#ifdef SIDM
  All.CrossSectionInternal = All.CrossSection*All.UnitMass_in_g/pow(All.UnitLength_in_cm,2);
#endif

  if(ThisTask==0)
    {
      fprintf(stdout, "\nHubble (internal units) = %g\n",All.Hubble);
      fprintf(stdout, "G (internal units) = %g\n",All.G);
      fprintf(stdout, "UnitMass_in_g = %g \n", All.UnitMass_in_g);
      fprintf(stdout, "UnitTime_in_s = %g \n", All.UnitTime_in_s);
      fprintf(stdout, "UnitVelocity_in_cm_per_s = %g \n", All.UnitVelocity_in_cm_per_s);
      fprintf(stdout, "UnitDensity_in_cgs = %g \n",All.UnitDensity_in_cgs);
      fprintf(stdout, "UnitEnergy_in_cgs = %g \n", All.UnitEnergy_in_cgs);
#ifdef SIDM
      fprintf(stdout, "CrossSection = %g \n", All.CrossSection);
      fprintf(stdout, "CrossSectionInternal = %g \n", All.CrossSectionInternal);
#if (CROSS_SECTION_TYPE == 2)
     fprintf(stdout, "YukawaVelocity vc = %g \n", All.YukawaVelocity); 
#endif
#endif
    }
  
  All.MinEgySpec  =(1.0/GAMMA_MINUS1)*(BOLTZMANN/PROTONMASS)*All.MinGasTemp;
  All.MinEgySpec *= All.UnitMass_in_g/All.UnitEnergy_in_cgs;


#ifdef SFR
  All.EgySpecCold= (1.0/GAMMA_MINUS1)*(BOLTZMANN/PROTONMASS) *All.TempClouds; 
  All.EgySpecCold *= All.UnitMass_in_g/All.UnitEnergy_in_cgs;
  
  All.EgySpecSN= (1.0/GAMMA_MINUS1)*(BOLTZMANN/PROTONMASS)*All.TempSupernova;
  All.EgySpecSN *= All.UnitMass_in_g/All.UnitEnergy_in_cgs;

  All.OverDensThresh= All.CritOverDensity * All.OmegaBaryon*3*All.Hubble*All.Hubble/(8*PI*All.G);
  All.PhysDensThresh= All.CritPhysDensity * PROTONMASS / HYDROGEN_MASSFRAC / All.UnitDensity_in_cgs;

  All.FeedbackEnergy= All.FactorSN/(1-All.FactorSN)*All.EgySpecSN;

  feedbackenergyinergs= All.FeedbackEnergy/All.UnitMass_in_g*(All.UnitEnergy_in_cgs*SOLAR_MASS);

  if(ThisTask==0)
    {
      printf("Feedback energy per formed solar mass in stars= %g  ergs\n",feedbackenergyinergs);
      printf("OverDensThresh= %g\nPhysDensThresh= %g\n", All.OverDensThresh, All.PhysDensThresh);
    }
#endif
}




/*
 *  Opens various log-files. On restart, the code 
 *  will append to the files.
 */
void open_outputfiles(void)
{
  char mode[2], buf[200];

  if(ThisTask!=0)  /* only the root processors writes to the log files */
    return;

 
  if(RestartFlag==0)
    strcpy(mode, "w");
  else
    strcpy(mode, "a");


  sprintf(buf,"%s%s", All.OutputDir, All.CpuFile);
  if(!(FdCPU= fopen(buf, mode)))
    {
      fprintf(stdout, "error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf,"%s%s", All.OutputDir, All.InfoFile);
  if(!(FdInfo= fopen(buf, mode)))
    {
      fprintf(stdout,"error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, All.EnergyFile);
  if(!(FdEnergy= fopen(buf, mode)))
    {
      fprintf(stdout, "error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, All.TimingsFile);
  if(!(FdTimings= fopen(buf, mode)))
    {
      fprintf(stdout, "error in opening file '%s'\n", buf);
      endrun(1);
    }

#ifdef SFR
  sprintf(buf, "%s%s", All.OutputDir, "sfr.txt");
  if(!(FdSfr= fopen(buf, mode)))
    {
      fprintf(stdout, "error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif
}


void close_outputfiles(void)
{
  if(ThisTask!=0)  /* only the root processors writes to the log files */
    return;

  fclose(FdCPU);
  fclose(FdInfo);
  fclose(FdEnergy);
  fclose(FdTimings);
}




/*
 *  This function parses the parameterfile in a simple way.
 *  Each paramater is defined by a keyword (`tag'), and can be
 *  either of type douple, int, or character string.
 *  The routine makes sure that each parameter appears 
 *  exactly once in the parameterfile.
 */
void read_parameter_file(char *fname)
{
#define DOUBLE 1
#define STRING 2
#define INT 3
#define MAXTAGS 300

  FILE *fd, *fdout;
  char buf[200], buf1[200], buf2[200], buf3[200];
  int  i, j, nt;
  int  id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  int  pnum, errorFlag=0;

  All.StarformationOn = All.MultiPhaseModelOn= 0; /* defaults */

  if(ThisTask==0)  /* read parameter file on process 0 */
    {
      nt=0;
    
      strcpy(tag[nt],"InitCondFile"); 
      addr[nt]=All.InitCondFile;
      id[nt++]=STRING;
      
      strcpy(tag[nt],"OutputDir"); 
      addr[nt]=All.OutputDir;
      id[nt++]=STRING;

      strcpy(tag[nt],"SnapshotFileBase"); 
      addr[nt]=All.SnapshotFileBase;
      id[nt++]=STRING;

      strcpy(tag[nt],"EnergyFile"); 
      addr[nt]=All.EnergyFile;
      id[nt++]=STRING;

      strcpy(tag[nt],"CpuFile"); 
      addr[nt]=All.CpuFile;
      id[nt++]=STRING;

      strcpy(tag[nt],"InfoFile"); 
      addr[nt]=All.InfoFile;
      id[nt++]=STRING;

      strcpy(tag[nt],"TimingsFile"); 
      addr[nt]=All.TimingsFile;
      id[nt++]=STRING;

      strcpy(tag[nt],"RestartFile"); 
      addr[nt]=All.RestartFile;
      id[nt++]=STRING;

      strcpy(tag[nt],"ResubmitCommand");
      addr[nt]=All.ResubmitCommand;
      id[nt++]=STRING;

      strcpy(tag[nt],"OutputListFilename");
      addr[nt]=All.OutputListFilename;
      id[nt++]=STRING;

      strcpy(tag[nt],"OutputListOn"); 
      addr[nt]=&All.OutputListOn;
      id[nt++]=INT;

      strcpy(tag[nt],"Omega0"); 
      addr[nt]=&All.Omega0;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"OmegaBaryon"); 
      addr[nt]=&All.OmegaBaryon;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"OmegaLambda"); 
      addr[nt]=&All.OmegaLambda;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"HubbleParam"); 
      addr[nt]=&All.HubbleParam;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"BoxSize"); 
      addr[nt]=&All.BoxSize;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"PeriodicBoundariesOn"); 
      addr[nt]=&All.PeriodicBoundariesOn;
      id[nt++]=INT;

      strcpy(tag[nt],"TimeOfFirstSnapshot"); 
      addr[nt]=&All.TimeOfFirstSnapshot;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"CpuTimeBetRestartFile"); 
      addr[nt]=&All.CpuTimeBetRestartFile;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"TimeBetStatistics"); 
      addr[nt]=&All.TimeBetStatistics;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"TimeBegin"); 
      addr[nt]=&All.TimeBegin;
      id[nt++]=DOUBLE;      

      strcpy(tag[nt],"TimeMax"); 
      addr[nt]=&All.TimeMax;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"TimeBetSnapshot"); 
      addr[nt]=&All.TimeBetSnapshot;
      id[nt++]=DOUBLE;
 
      strcpy(tag[nt],"UnitVelocity_in_cm_per_s"); 
      addr[nt]=&All.UnitVelocity_in_cm_per_s;
      id[nt++]=DOUBLE;
  
      strcpy(tag[nt],"UnitLength_in_cm"); 
      addr[nt]=&All.UnitLength_in_cm;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"UnitMass_in_g"); 
      addr[nt]=&All.UnitMass_in_g;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"MaxNodeMove"); 
      addr[nt]=&All.MaxNodeMove;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"TreeUpdateFrequency"); 
      addr[nt]=&All.TreeUpdateFrequency;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"ErrTolIntAccuracy"); 
      addr[nt]=&All.ErrTolIntAccuracy;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"ErrTolDynamicalAccuracy"); 
      addr[nt]=&All.ErrTolDynamicalAccuracy;
      id[nt++]=DOUBLE;



      strcpy(tag[nt],"ErrTolVelScale"); 
      addr[nt]=&All.ErrTolVelScale;
      id[nt++]=DOUBLE;
       
      strcpy(tag[nt],"ErrTolTheta"); 
      addr[nt]=&All.ErrTolTheta;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"ErrTolForceAcc"); 
      addr[nt]=&All.ErrTolForceAcc;
      id[nt++]=DOUBLE;
      
      strcpy(tag[nt],"MinGasHsmlFractional"); 
      addr[nt]=&All.MinGasHsmlFractional;
      id[nt++]=DOUBLE;
      
      strcpy(tag[nt],"MaxSizeTimestep"); 
      addr[nt]=&All.MaxSizeTimestep;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"MinSizeTimestep"); 
      addr[nt]=&All.MinSizeTimestep;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"ArtBulkViscConst"); 
      addr[nt]=&All.ArtBulkViscConst;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"CourantFac"); 
      addr[nt]=&All.CourantFac;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"DesNumNgb"); 
      addr[nt]=&All.DesNumNgb;
      id[nt++]=INT;

      strcpy(tag[nt],"MaxNumNgbDeviation"); 
      addr[nt]=&All.MaxNumNgbDeviation;
      id[nt++]=INT;

      strcpy(tag[nt],"ComovingIntegrationOn"); 
      addr[nt]=&All.ComovingIntegrationOn;
      id[nt++]=INT;

      strcpy(tag[nt],"ICFormat"); 
      addr[nt]=&All.ICFormat;
      id[nt++]=INT;

      strcpy(tag[nt],"NumFilesPerSnapshot"); 
      addr[nt]=&All.NumFilesPerSnapshot;
      id[nt++]=INT;

      strcpy(tag[nt],"NumFilesWrittenInParallel"); 
      addr[nt]=&All.NumFilesWrittenInParallel;
      id[nt++]=INT;

      strcpy(tag[nt],"ResubmitOn"); 
      addr[nt]=&All.ResubmitOn;
      id[nt++]=INT;

      strcpy(tag[nt],"CoolingOn"); 
      addr[nt]=&All.CoolingOn;
      id[nt++]=INT;

      strcpy(tag[nt],"TypeOfTimestepCriterion"); 
      addr[nt]=&All.TypeOfTimestepCriterion;
      id[nt++]=INT;

      strcpy(tag[nt],"TypeOfOpeningCriterion"); 
      addr[nt]=&All.TypeOfOpeningCriterion;
      id[nt++]=INT;

      strcpy(tag[nt],"TimeLimitCPU"); 
      addr[nt]=&All.TimeLimitCPU;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"DomainUpdateFrequency"); 
      addr[nt]=&All.DomainUpdateFrequency;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"SofteningHalo"); 
      addr[nt]=&All.SofteningHalo;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"SofteningDisk"); 
      addr[nt]=&All.SofteningDisk;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"SofteningBulge"); 
      addr[nt]=&All.SofteningBulge;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"SofteningGas"); 
      addr[nt]=&All.SofteningGas;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"SofteningStars"); 
      addr[nt]=&All.SofteningStars;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"SofteningHaloMaxPhys"); 
      addr[nt]=&All.SofteningHaloMaxPhys;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"SofteningDiskMaxPhys"); 
      addr[nt]=&All.SofteningDiskMaxPhys;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"SofteningBulgeMaxPhys"); 
      addr[nt]=&All.SofteningBulgeMaxPhys;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"SofteningGasMaxPhys"); 
      addr[nt]=&All.SofteningGasMaxPhys;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"SofteningStarsMaxPhys"); 
      addr[nt]=&All.SofteningStarsMaxPhys;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"BufferSize"); 
      addr[nt]=&All.BufferSize;
      id[nt++]=INT;

      strcpy(tag[nt],"PartAllocFactor"); 
      addr[nt]=&All.PartAllocFactor;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"TreeAllocFactor"); 
      addr[nt]=&All.TreeAllocFactor;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"GravityConstantInternal"); 
      addr[nt]=&All.GravityConstantInternal;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitGasTemp"); 
      addr[nt]=&All.InitGasTemp;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"MinGasTemp"); 
      addr[nt]=&All.MinGasTemp;
      id[nt++]=DOUBLE;
      
#ifdef SIDM
      strcpy(tag[nt],"CrossSection"); 
      addr[nt]=&All.CrossSection;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"RandomSeed1"); 
      addr[nt]=&All.Seed1;
      id[nt++]=INT;

      strcpy(tag[nt],"RandomSeed2"); 
      addr[nt]=&All.Seed2;
      id[nt++]=INT;

      strcpy(tag[nt],"ProbabilityTol"); 
      addr[nt]=&All.ProbabilityTol;
      id[nt++]=DOUBLE;

#if (CROSS_SECTION_TYPE == 2 || CROSS_SECTION_TYPE == 4)
      strcpy(tag[nt],"YukawaVelocity"); 
      addr[nt]=&All.YukawaVelocity;
      id[nt++]=DOUBLE;
#elif (CROSS_SECTION_TYPE == 3)
      strcpy(tag[nt], "CrossSectionVelScale");
      addr[nt]=&All.CrossSectionVelScale;
      id[nt++]=DOUBLE;
      
      strcpy(tag[nt], "CrossSectionPowLaw");
      addr[nt]=&All.CrossSectionPowLaw;
      id[nt++]=DOUBLE;
#endif

#endif

#ifdef REFLECTIONBOUNDARY
      strcpy(tag[nt],"ReflectionBoundary"); 
      addr[nt]=&All.ReflectionRadius;
      id[nt++]=DOUBLE;
#endif

#ifdef MOREPARAMS
      strcpy(tag[nt],"StarformationOn"); 
      addr[nt]=&All.StarformationOn;
      id[nt++]=INT;

      strcpy(tag[nt],"MultiPhaseModelOn"); 
      addr[nt]=&All.MultiPhaseModelOn;
      id[nt++]=INT;

      strcpy(tag[nt],"FactorSFR"); 
      addr[nt]=&All.FactorSFR;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"FactorSN"); 
      addr[nt]=&All.FactorSN;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"FactorEVP"); 
      addr[nt]=&All.FactorEVP;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"TempSupernova"); 
      addr[nt]=&All.TempSupernova;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"TempClouds"); 
      addr[nt]=&All.TempClouds;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"CritOverDensity"); 
      addr[nt]=&All.CritOverDensity;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"CritPhysDensity"); 
      addr[nt]=&All.CritPhysDensity;
      id[nt++]=DOUBLE;
#endif


      if((fd=fopen(fname,"r")))
	{
	  sprintf(buf,"%s%s", fname, "-usedvalues");
	  if(!(fdout=fopen(buf,"w")))
	    {
	      fprintf(stdout,"error opening file '%s' \n", buf);
 	      errorFlag=1;
	    }
	  else
	    {
	      while(!feof(fd))
		{
		  fgets(buf,200,fd);
		  if(sscanf(buf,"%s%s%s", buf1, buf2, buf3)<2)
		    continue;
		  
		  if(buf1[0]=='%')
		    continue;
		  
		  for(i=0,j=-1;i<nt;i++)
		    if(strcmp(buf1,tag[i])==0)
		      {
			j=i;
			tag[i][0]=0;
			break;
		      }
		  
		  if(j>=0)
		    {
		      switch(id[j])
			{
			case DOUBLE:
			  *((double*)addr[j])=atof(buf2); 
			  fprintf(fdout,"%-35s%g\n",buf1,*((double*)addr[j]));
			  break;
			case STRING:
			  strcpy(addr[j],buf2);
			  fprintf(fdout,"%-35s%s\n",buf1,buf2);
			  break;
			case INT:
			  *((int*)addr[j])=atoi(buf2);
			  fprintf(fdout,"%-35s%d\n",buf1,*((int*)addr[j]));
			  break;
			}
		    }
		  else
		    {
		      fprintf(stdout,"Error in file %s:   Tag '%s' not allowed or multiple defined.\n",fname,buf1);
		      errorFlag=1;
		    }
		}
	      fclose(fd);
	      fclose(fdout);
	      
	      i= strlen(All.OutputDir);
	      if(i>0)
		if(All.OutputDir[i-1]!='/')
		  strcat(All.OutputDir, "/");

	      sprintf(buf1, "%s%s", fname, "-usedvalues");
	      sprintf(buf2, "%s%s", All.OutputDir, "parameters_out");
	      rename(buf1, buf2);
	    }
	}
      else
	{
	  fprintf(stdout,"Parameter file %s not found.\n",fname);
	  errorFlag=1;
	}
      
      
      for(i=0;i<nt;i++)
	{
	  if(*tag[i])
	    {
	      fprintf(stdout,"Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
	      errorFlag=1;
	    }
	}

      if(All.OutputListOn && errorFlag==0)
	errorFlag+= read_outputlist(All.OutputListFilename);
      else
	All.OutputListLength=0;
    }

  MPI_Bcast(&errorFlag, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if(errorFlag)
    {
      MPI_Finalize();
      exit(0);
    }

  /* now communicate the relevant parameters to the other processes */
  MPI_Bcast(&All, sizeof(struct global_data_all_processes), MPI_BYTE, 0, MPI_COMM_WORLD);


  for(pnum=0; All.NumFilesWrittenInParallel>(1<<pnum); pnum++);

  if(All.NumFilesWrittenInParallel!=(1<<pnum))
    {
      if(ThisTask==0)
	fprintf(stdout,"NumFilesWrittenInParallel MUST be a power of 2\n");
      endrun(0);
    }

  if(All.NumFilesWrittenInParallel>NTask)
    {
      if(ThisTask==0)
	fprintf(stdout,"NumFilesWrittenInParallel MUST be smaller than number of processors\n");
      endrun(0);
    }

#ifdef PERIODIC
  if(All.PeriodicBoundariesOn==0)
    {
      if(ThisTask==0)
	{
	  fprintf(stdout,"Code was compiled with periodic boundary conditions switched on.\n");
	  fprintf(stdout,"You must set `PeriodicBoundariesOn=1', or recompile the code.\n");
	}	  
      endrun(0);
    }
#else
  if(All.PeriodicBoundariesOn==1)
    {
      if(ThisTask==0)
	{
	  fprintf(stdout,"Code was compiled with periodic boundary conditions switched off.\n");
	  fprintf(stdout,"You must set `PeriodicBoundariesOn=0', or recompile the code.\n");
	}	  
      endrun(0);
    }
#endif

#ifdef COOLING
  if(All.CoolingOn==0)
    {
      if(ThisTask==0)
	{
	  fprintf(stdout,"Code was compiled with cooling switched on.\n");
	  fprintf(stdout,"You must set `CoolingOn=1', or recompile the code.\n");
	}	  
      endrun(0);
    }
#else
  if(All.CoolingOn==1)
    {
      if(ThisTask==0)
	{
	  fprintf(stdout,"Code was compiled with cooling switched off.\n");
	  fprintf(stdout,"You must set `CoolingOn=0', or recompile the code.\n");
	}	  
      endrun(0);
    }
#endif

#ifdef VELDISP
  if(All.TypeOfTimestepCriterion < 2 || All.TypeOfTimestepCriterion > 4)
    {
      if(ThisTask==0)
	{
	  fprintf(stdout,"Code was compiled with computation of dark matter\n");
	  fprintf(stdout,"velocity dispersion. However, the selected timestep\n");
	  fprintf(stdout,"criterion does not use it..! So, either use another\n");
	  fprintf(stdout,"timestep criterion, or switch of VELDISP\n");
	}
      endrun(0);
    }
#else
  if(All.TypeOfTimestepCriterion >= 2 && All.TypeOfTimestepCriterion <= 4)
    {
      if(ThisTask==0)
	{
	  fprintf(stdout,"The specified timestep criterion\n");
	  fprintf(stdout,"requires that the code is compiled with the\n");
	  fprintf(stdout,"VELDISP option.\n");
	}
      endrun(0);
    }
#endif

#ifdef SFR
 
#ifndef MOREPARAMS
  if(ThisTask==0)
    {
      fprintf(stdout,"Code was compiled with SFR, but not with MOREPARAMS.\n");
      fprintf(stdout,"This is not allowed.\n");
    }	  
  endrun(0);
#endif

  if(All.StarformationOn==0)
    {
      if(ThisTask==0)
	{
	  fprintf(stdout,"Code was compiled with star formation switched on.\n");
	  fprintf(stdout,"You must set `StarformationOn=1', or recompile the code.\n");
	}	  
      endrun(0);
    }
  if(All.CoolingOn==0)
    {
      if(ThisTask==0)
	{
	  fprintf(stdout,"You try to use the code with star formation enabled,\n");
	  fprintf(stdout,"but you did not switch on cooling.\nThis mode is not supported.\n");
	}	  
      endrun(0);
    }
#else
  if(All.StarformationOn==1)
    {
      if(ThisTask==0)
	{
	  fprintf(stdout,"Code was compiled with star formation switched off.\n");
	  fprintf(stdout,"You must set `StarformationOn=0', or recompile the code.\n");
	}	  
      endrun(0);
    }
#endif

#ifdef CLOUDS
  if(All.MultiPhaseModelOn==0)
   { 
     if(ThisTask==0)
       {
	 fprintf(stdout,"Code was compiled with multi-phase model enabled,\n");
	 fprintf(stdout,"but you did not switch it on in the parameterfile.\n");
	 fprintf(stdout,"Tou must set `MultiPhaseModel=1' or recompile the code\n");
       }	  
     endrun(0);
   }
  if(All.StarformationOn==0)
    {
      if(ThisTask==0)
	{
	  fprintf(stdout,"Code was compiled with multi-phase model enabled,\n");
	  fprintf(stdout,"but you did not switch on star formation\n");
	  fprintf(stdout,"This combination is not supported.\n");
	}	  
      endrun(0);
    }
#else
  if(All.MultiPhaseModelOn==1)
   { 
     if(ThisTask==0)
       {
	 fprintf(stdout,"Code was compiled with multi-phase model switched off.\n");
	 fprintf(stdout,"To use this option, you need to recompile.\n");
       }	  
     endrun(0);
   }
#endif


#undef DOUBLE 
#undef STRING 
#undef INT 
#undef MAXTAGS
}


/* this function reads a table with a list of desired output
 * times. The table does not have to be ordered in any way,
 * but may not contain more than MAXLEN_OUTPUTLIST entries.
 */
int read_outputlist(char *fname)
{
  FILE *fd;
  
  if(!(fd=fopen(fname,"r")))
    {
      fprintf(stdout, "can't read output list in file '%s'\n", fname);
      return 1;
    }
  
  All.OutputListLength= 0;
  do
    {
      if(fscanf(fd, " %lg ", &All.OutputListTimes[All.OutputListLength])==1)
	All.OutputListLength++;
      else
	break;
    }
  while(All.OutputListLength < MAXLEN_OUTPUTLIST);

  fclose(fd);

  printf("\nfound %d times in output-list.\n", All.OutputListLength);

  return 0;
}

/* this function returns the next output time
 * in case a table of output times is used.
 */
double find_next_outputtime(double time)
{
  int i;
  double next;
  
  for(i=0, next=MAX_REAL_NUMBER; i<All.OutputListLength; i++)
    {
      if(All.OutputListTimes[i] > time)
	if(next > All.OutputListTimes[i])
	  next= All.OutputListTimes[i];
    }
  
  return next;
}



/* Here the lookup table for the kernel of the SPH calculation
 * is initialized. 
 */ 
void set_sph_kernel(void)
{
  int i;
#ifdef SIDM
  double theta;
#endif

  for(i=0;i<=KERNEL_TABLE+1;i++)
    KernelRad[i] = ((double)i)/KERNEL_TABLE;

  Kernel[KERNEL_TABLE+1] = KernelDer[KERNEL_TABLE+1]= 0;

  for(i=0;i<=KERNEL_TABLE;i++)
    {
      if(KernelRad[i]<=0.5)
	  {
	    Kernel[i] = 8/PI *(1-6*KernelRad[i]*KernelRad[i]*(1-KernelRad[i]));
	    KernelDer[i] = 8/PI *( -12*KernelRad[i] + 18*KernelRad[i]*KernelRad[i]);
	  }
      else
	  {
	    Kernel[i] = 8/PI * 2*(1-KernelRad[i])*(1-KernelRad[i])*(1-KernelRad[i]);
	    KernelDer[i] = 8/PI *( -6*(1-KernelRad[i])*(1-KernelRad[i]));
	  }
    }

#ifdef SIDM
  // Lookup Table for SIDM scattering
  /*
  for(i=0 ; i<=KERNEL_TABLE ; i++){
    theta = 2*PI*i/KERNEL_TABLE;
    CosTable[i] = cos(theta);
    SinTable[i] = sin(theta);
  }
  */
#endif
}








