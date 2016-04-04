#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"



 
/* This function reads initial conditions that are in the default file
 * format of Gadget, i.e. snapshot files can be used as input files.
 * However, when a snapshot file is used as input, not all the
 * information in the header is used: THE STARTING TIME NEEDS TO BE
 * SET IN THE PARAMETERFILE.
 *
 * Alternatively, the code can be started with restartflag==2, then
 * snapshots from the code can be used as initial conditions-files 
 * without having to change the parameterfile.
 *
 * For gas particles, only the internal energy is read, the density
 * and mean molecular weight will be recomputed by the code.
 *
 * Only Task=0 reads the file(s), and distributes the particles evenly
 * onto the processors.
 *
 * When InitGasTemp>0 is given, the gas temperature will be initialzed
 * to this value assuming a mean colecular weight of 1.  */
void read_ic(char *fname)  
{
#define BLOCKSIZE 10000
#define SKIP  {my_fread(&blksize,sizeof(int4byte),1,fd);printf("blksize= %d\n", blksize);}
  FILE *fd=0;
  int   i,k,n,files,ntot_withmasses;
  float *fp,dummy[BLOCKSIZE][3];
  double u_init;
#ifdef SFR
  double original_gas_mass, mass, masstot;
#endif
  int   pc,pc_here,nall;
  int   type,pr,left,groupsize,posvel;
  MPI_Status status;
  int4byte   blksize;
  char  buf[100];
  int   n_for_this_task, n_in_file;
  int4byte *ip;
  
  pc= 1;
  NumPart= 0;
  N_gas= 0;
      
  for(files=0; files<All.NumFilesPerSnapshot; files++)   
    {
      if(ThisTask==0)
	{
	  if(All.NumFilesPerSnapshot>1)
	    sprintf(buf,"%s.%d", fname, files); 
	  else
	    sprintf(buf,"%s", fname); 
 
	  if(!(fd=fopen(buf,"r")))
	    {
	      printf("can't open file `%s' for reading initial conditions.\n",buf);
	      endrun(123);
	    }
	  
	  printf("Reading file '%s' ...\n", buf);

	  SKIP;
	  if(blksize!=256)
	    {
	      printf("incorrect header format (1)\n");
	      fflush(stdout);
	      endrun(888);
	    }
	  my_fread(&header1, sizeof(header1), 1, fd);
	  SKIP;
	  if(blksize!=256)
	    {
	      printf("incorrect header format (2)\n");
	      fflush(stdout);
	      endrun(889);
	    }

	  if(files==0)
	    {
	      if(All.NumFilesPerSnapshot==1)
		for(i=0;i<6;i++)
		  header1.npartTotal[i] = header1.npart[i];

	      All.TotN_gas=  header1.npartTotal[0];
	      All.TotN_halo= header1.npartTotal[1];
	      All.TotN_disk= header1.npartTotal[2];
	      All.TotN_bulge=header1.npartTotal[3];
	      All.TotN_stars=header1.npartTotal[4];

	      for(i=0;i<6;i++)
		All.MassTable[i] = header1.mass[i];

	      All.TotNumPart= All.TotN_gas+ All.TotN_halo+ 
		All.TotN_disk+ All.TotN_bulge+ All.TotN_stars;
	    }

	  printf("file `%s' contains %d particles.\n", 
		 buf, header1.npart[0]+header1.npart[1]+header1.npart[2]+header1.npart[3]+header1.npart[4]);
	  printf("Type 0 (gas):   %d\nType 1 (halo):  %d\nType 2 (disk):  %d\nType 3 (bulge): %d\nType 4 (stars): %d\n",
		 All.TotN_gas, All.TotN_halo, All.TotN_disk, All.TotN_bulge, All.TotN_stars);
	  printf("\nMassTable:\nGas: %g\nHalo: %g\nDisk: %g\nBulge: %g\nStars: %g\n",
		 All.MassTable[0], All.MassTable[1], All.MassTable[2], All.MassTable[3], All.MassTable[4]);
	}

      MPI_Bcast(&header1, sizeof(struct io_header_1), MPI_BYTE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&All, sizeof(struct global_data_all_processes), MPI_BYTE, 0, MPI_COMM_WORLD);

      if(RestartFlag==2)
	{
	  All.Time = All.TimeBegin = header1.time;
	}


      for(i=0, ntot_withmasses=0; i<5; i++)
	{
	  if(All.MassTable[i]==0)
	    ntot_withmasses+=  header1.npart[i];
	}

      if(files==0)
	{
	  All.MaxPart=    All.PartAllocFactor * (All.TotNumPart/NTask);    /* sets the maximum number of particles that may */
 	  All.MaxPartSph= All.PartAllocFactor * (All.TotN_gas/NTask);      /* sets the maximum number of particles that may 
									      reside on a processor */
	  allocate_memory();
	}

      if(files>0)
	{
	  /* to collect the gas particles all at the beginning, we move the collisionless
             particles such that a gap of the right size is created */

	  for(type=0,nall=0; type<5; type++)
	    {
	      n_in_file=header1.npart[type];
	      
	      n_for_this_task = n_in_file / NTask;
	      if(ThisTask <  (n_in_file % NTask))
		n_for_this_task++;
	      
	      nall+= n_for_this_task;
	    }
	  
	  memmove(&P[1+N_gas+nall], &P[1+N_gas], (NumPart-N_gas)*sizeof(struct particle_data)); 
	}


  /*  posvel=  0 for positions, 1 for velocities, 
       *  2 for ID, 3 for masses, 4 for internal energy, 
       *  5 for density  
       *  6 for electron abundance
       *  7 for neutral hydrogen abundance
       *  8 for formed stellar mass
       *  9 for for cold cloud mass (only for CLOUDS)
       *  10 for smoothing length
       *  11 for current star formation rate
       *  12 for mean stellar age  (only for STELLARAGE)
       */	

      for(posvel=0; posvel<=12; posvel++)  
	{
	  if(RestartFlag==0 && posvel>4)
	    continue;

	  if(posvel==3 && ntot_withmasses==0) /* no mass-block in file */
	    continue;
	  
	  if(posvel==4 && header1.npart[0]==0) /* no energy-block in file */
	    continue;

	  if(posvel==5 && header1.npart[0]==0) /* no density-block in file */
	    continue;

	  if(All.CoolingOn==0 && (posvel==6 || posvel==7))
	    continue;

          if(posvel==10 && header1.npart[0]==0) /* no hsml-block in file */
            continue;

	  if(All.StarformationOn==0 && (posvel==8 || posvel==9 || posvel==11 || 
					posvel==12 || posvel==13 || posvel==14))
	    continue;

	  if(All.MultiPhaseModelOn==0 && posvel==9)
	    continue;

#ifndef STELLARAGE
	  if(posvel==12)
	    continue;
#endif
	  if(ThisTask==0)
	    SKIP;
	    
	  pc_here=pc;

	  for(type=0; type<5; type++)
	    {
	      n_in_file=header1.npart[type];

	      n_for_this_task = n_in_file / NTask;
	      if(ThisTask <  (n_in_file % NTask))
		n_for_this_task++;

	      if(posvel>=4 && type>0 && type<4) 
		{
		  pc_here+= n_for_this_task;
		  continue;
		}

	      if(posvel>=4 && posvel<12 && type==4) 
		{
		  pc_here+= n_for_this_task;
		  continue;
		}



	      for(pr=0; pr<NTask; pr++) /* go through all processes, note: pr is the receiving process */
		{
		  if(ThisTask==0 || ThisTask==pr)
		    {
		      n = n_for_this_task;    /* number of particles for this process */
		      
		      if(ThisTask==0 && pr>0)
			MPI_Recv(&n, 1, MPI_INT, pr , TAG_N, MPI_COMM_WORLD, &status);
		  
		      if(ThisTask==pr && pr>0)
			MPI_Send(&n, 1, MPI_INT, 0 , TAG_N, MPI_COMM_WORLD);

		  
		      left=n;
		      
		      while(left>0)
			{
			  if(left>BLOCKSIZE)   /* restrict transmission size to buffer length */
			    groupsize=BLOCKSIZE;
			  else
			    groupsize=left;
			  
			  if(ThisTask==0)
			    {
			      if(posvel<2)
				my_fread(&dummy[0][0],sizeof(float), 3*groupsize ,fd);
			      else 
				{
				  if(posvel!=3) /* ID, or internal energy, or everything else */
				    my_fread(&dummy[0][0],sizeof(int4byte),  groupsize ,fd);
				  else /* masses */
				    {
				      if(All.MassTable[type]==0) 
					my_fread(&dummy[0][0],sizeof(float),  groupsize ,fd);
				      else
					{
					  for(i=0, fp=&dummy[0][0]; i<groupsize; i++, fp++)
					    *fp = All.MassTable[type];
					}
				    }
				}
			    }
		      
			  if(ThisTask==0 && pr!=0)
			    {
			      if(posvel<2)
				MPI_Send(&dummy[0][0], 3*groupsize, MPI_FLOAT, pr, TAG_ANY, MPI_COMM_WORLD);
			      else
				MPI_Send(&dummy[0][0], groupsize*sizeof(int4byte), MPI_BYTE, pr, TAG_ANY, MPI_COMM_WORLD);
			    }
		      
			  if(ThisTask!=0 && pr!=0)
			    {
			      if(posvel<2)
				MPI_Recv(&dummy[0][0], 3*groupsize, MPI_FLOAT, 0 , TAG_ANY, MPI_COMM_WORLD, &status);
			      else
				MPI_Recv(&dummy[0][0], groupsize*sizeof(int4byte), MPI_BYTE, 0 , TAG_ANY, MPI_COMM_WORLD, &status);
			    }
		      
			  if(ThisTask==pr)
			    {
			      for(i=0,ip=(int4byte *)&dummy[0][0], fp=&dummy[0][0]; i<groupsize; i++)
				{
				  switch(posvel)
				    {
				    case 0:
				      P[pc_here].Type=type;
				      for(k=0;k<3;k++)
					P[pc_here].Pos[k]=dummy[i][k]; 
				      pc_here++;
				      break;
				    case 1:
				      for(k=0;k<3;k++)
					P[pc_here].Vel[k]=dummy[i][k]; 
				      pc_here++;
				      break;
				    case 2:
				      P[pc_here].ID  =ip[i];    /* now set ID */
				      pc_here++;
				      break;
				    case 3:
				      P[pc_here].Mass= fp[i];    /* now set mass */
				      pc_here++;
				      break;
				    case 4:
				      SphP[pc_here].EgySpec= fp[i];    /* now set internal energy for gas particle */
				      pc_here++;
				      break;
				    case 5:
				      SphP[pc_here].Density= SphP[pc_here].DensityPred= fp[i];
				      pc_here++;
				      break;
#ifdef COOLING
				    case 6:
				      SphP[pc_here].Ne = fp[i];
				      pc_here++;
				      break;
				    case 7:
				      /* skip neutral hydrogen fraction */
				      pc_here++;
				      break;
#endif
				      
#ifdef SFR
				    case 8:
				      SphP[pc_here].FormedStellarMass=  fp[i];
				      pc_here++;
				      break;
#ifdef CLOUDS
				    case 9:
				      SphP[pc_here].CloudMass=   fp[i];
				      pc_here++;
				      break;
#endif
#endif
				    case 10:
				      SphP[pc_here].Hsml= fp[i];
				      pc_here++;
				      break;
#ifdef SFR
				    case 11:
				      /* skip star formation rate */
				      pc_here++;
				      break;
#endif
				      
#ifdef STELLARAGE
				    case 12:
				      if(type==0)
					{
					  P[pc_here].MeanStellarAge= fp[i]*SphP[pc_here].FormedStellarMass;
					  pc_here++;
					}
				      if(type==4)
					{
					  P[pc_here].MeanStellarAge= fp[i]*P[pc_here].Mass; 
					  pc_here++;
					}
				      break;
#endif
				    }
				}
			    }		      
			  
			  left-=groupsize;
			}
		    }	      
		  			      
		  MPI_Barrier(MPI_COMM_WORLD);
		}
	    }

	  if(ThisTask==0)
	    SKIP;
	}
      
      for(type=0; type<5; type++)
	{
	  n_in_file=header1.npart[type];
	  
	  n_for_this_task = n_in_file / NTask;
	  if(ThisTask <  (n_in_file % NTask))
	    n_for_this_task++;
	  
	  NumPart  += n_for_this_task;

	  if(type==0)
	    {
	      pc +=    n_for_this_task;
	      N_gas += n_for_this_task;
	    }
	}

   
      if(ThisTask==0)
	fclose(fd);
    }
  

  /* this makes sure that masses are initialized in the case that the mass-block
     is completely empty */
  for(i=1;i<=NumPart;i++)  
    {
      if(All.MassTable[P[i].Type]!=0)
	P[i].Mass= All.MassTable[P[i].Type];
    } 

#ifdef SFR
  if(RestartFlag==0)
    {
      if(All.MassTable[4]==0 && All.MassTable[0]>0)
	{
	  All.MinGasMass= DISSOLVE_FRAC * All.MassTable[0];
	  All.MassTable[0]= 0;
	  All.MassTable[4]= 0;
	}
      else
	{
	  for(i=1, mass=0; i<=NumPart; i++)  
	    {
	      if(All.MassTable[P[i].Type]==0 || All.MassTable[P[i].Type]==4)
		mass+= P[i].Mass;
	    } 
      
	  MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	  original_gas_mass= masstot/(All.TotN_gas +  All.TotN_stars);
	  All.MinGasMass= DISSOLVE_FRAC * original_gas_mass;
	}
    }
  if(RestartFlag==2)
    {
      for(i=1, mass=0; i<=NumPart; i++)  
	{
	  if(All.MassTable[P[i].Type]==0 || All.MassTable[P[i].Type]==4)
	    mass+= P[i].Mass;
	} 
      
      MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      original_gas_mass= masstot/(All.TotN_gas +  All.TotN_stars);
      if(ThisTask==0)
	printf("original_gas_mass= %g\n", original_gas_mass);
      All.MinGasMass= DISSOLVE_FRAC * original_gas_mass;
    }
#endif

  if(RestartFlag==0)
    {
      if(All.InitGasTemp>0)
	{
	  u_init = (1.0/GAMMA_MINUS1)*(BOLTZMANN/PROTONMASS)*All.InitGasTemp;
	  u_init*= All.UnitMass_in_g/All.UnitEnergy_in_cgs;      /* unit conversion */  
	  
	  for(i=1;i<=N_gas;i++) 
	    {
	      if(SphP[i].EgySpec==0)
		SphP[i].EgySpec= u_init;
	    }
	}
    }

  MPI_Barrier(MPI_COMM_WORLD);
  if(ThisTask==0)
    {
      fprintf(stdout,"reading done.\n");
      fflush(stdout);
    }

  if(ThisTask==0)
    {
      printf("Total number of particles :  %d\n\n", All.TotNumPart);
      fflush(stdout);
    }
}










