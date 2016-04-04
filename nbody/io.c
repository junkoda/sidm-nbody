#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"
#include "tags.h"



/* This wrapper function select the desired output 
 * routine for snapshot files.
 */
void savepositions(int num)
{
  double t0,t1;

  t0=second();

  if(ThisTask==0)
    printf("\nwriting snapshot file... \n");
 
  savepositions_ioformat1(num);

  if(ThisTask==0)
    printf("done with snapshot.\n");
  
  t1=second();

  All.CPU_Snapshot+= timediff(t0,t1);
}





/* This function writes a snapshot of the particle ditribution to
 * one or several files using Gadget's default file format.
 * If NumFilesPerSnapshot>1, the snapshot is distributed into several files,
 * which are written simultaneously. Each file contains data from a group
 * of processors of size NTasl/NumFilesPerSnapshot.
 *
 * Each snapshot file contains a header first, then particle positions, 
 * velocities and ID's.
 * Then particle masses are written for those particle types with zero entry in
 * MassTable.
 * After that, first the internal energies u, and then the density is written
 * for the SPH particles.
 * Finally, if cooling is enabled, the mean molecular weight is written for the gas
 * particles. 
 */
void savepositions_ioformat1(int num)
{
#define BLOCKSIZE 10000
  FILE  *fd=0;
  char  buf[300];
  int   i,p,k,ntot=0,type,ntot_withmasses=0;
  float      block[BLOCKSIZE][3];
  float      blockmass[BLOCKSIZE];
  float      blocku[BLOCKSIZE];
  int4byte   blockid[BLOCKSIZE];
  int   pc;
  int   flag;
  int   posvel;
  MPI_Status status;
  int   n_type[5], ntot_type[5], ntot_type_all[5], nn[5];
  int4byte dummy;
  int   nprocgroup;
  int   masterTask;
  double a3inv=1.0;
#ifdef COOLING
  double ne, nh0;
#endif

  if(num<0)
    num=1000+num;

  nprocgroup=NTask/All.NumFilesPerSnapshot;
  
  if((NTask % nprocgroup))
    {
      printf("Fatal error.\nNumber of processors must be a multiple of All.NumFilesPerSnapshot.\n");
      endrun(213);
    }

  if(All.ComovingIntegrationOn)
    a3inv=  1/(All.Time*All.Time*All.Time);

  masterTask=(ThisTask/nprocgroup)*nprocgroup;


  if(ThisTask==masterTask)
    {
      if(All.NumFilesPerSnapshot>1)
	sprintf(buf,"%s%s_%03d.%d",All.OutputDir,All.SnapshotFileBase,num,ThisTask/nprocgroup); 
      else
	sprintf(buf,"%s%s_%03d",All.OutputDir,All.SnapshotFileBase,num); 
    
      if(!(fd=fopen(buf,"w")))
	{ 
	  fprintf(stdout,"Error. Can't write in file '%s'\n",buf); 
	  endrun(10); 
	}
    }  

  for(i=0;i<5;i++)
    n_type[i]=0;

  for(i=1;i<=NumPart;i++)
    n_type[P[i].Type]++;

  MPI_Reduce(n_type, ntot_type_all, 5, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Bcast(ntot_type_all, 5, MPI_INT, 0, MPI_COMM_WORLD);

  if(All.NumFilesPerSnapshot==1)
    {
      MPI_Reduce(n_type, ntot_type, 5, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Bcast(ntot_type, 5, MPI_INT, 0, MPI_COMM_WORLD);

      for(i=0,ntot=0,ntot_withmasses=0;i<5;i++)
	{
	  ntot+=ntot_type[i];
	  if(All.MassTable[i]==0)
	    ntot_withmasses+=ntot_type[i];
	}
    }
  else
    {
      if(ThisTask==masterTask)
	{
	  for(i=0;i<5;i++)
	    ntot_type[i]=n_type[i];

	  for(p=1;p<nprocgroup;p++)
	    {
	      MPI_Recv(&nn[0], 5, MPI_INT, masterTask+p, TAG_ANY, MPI_COMM_WORLD, &status);
	      for(i=0;i<5;i++)
		ntot_type[i]+=nn[i];
	    }	  
  
	  for(i=0,ntot=0,ntot_withmasses=0;i<5;i++)
	    {
	      ntot+=ntot_type[i];
	      if(All.MassTable[i]==0)
		ntot_withmasses+=ntot_type[i];
	    }
	}
      else
	{
	  MPI_Send(&n_type[0], 5, MPI_INT, masterTask, TAG_ANY, MPI_COMM_WORLD);
	}
    }



  if(ThisTask==masterTask)
    {
      for(i=0;i<5;i++)
	{
	  header1.npart[i]= ntot_type[i];
	  header1.npartTotal[i]= ntot_type_all[i];
	}
      header1.npart[5]=0;
      header1.npartTotal[5]=0;

      for(i=0;i<6;i++)
	header1.mass[i]=All.MassTable[i];

      header1.time=All.Time;

      if(All.ComovingIntegrationOn)
	header1.redshift=1.0/All.Time -1;
      else
	header1.redshift=0;  

      header1.flag_sfr= 0;
      header1.flag_feedback= 0;
      header1.flag_cooling= 0;
      header1.flag_multiphase= 0;
      header1.flag_stellarage= 0;

#ifdef COOLING
      header1.flag_cooling= 1;
#endif
#ifdef SFR
      header1.flag_sfr= 1;
      header1.flag_feedback= 1;
#ifdef CLOUDS
      header1.flag_multiphase= 1;
#endif
#ifdef STELLARAGE
      header1.flag_stellarage= 1;
#endif
#ifdef SFRHISTOGRAM
      header1.flag_sfrhistogram= 1;
#endif
#endif

      header1.num_files= All.NumFilesPerSnapshot;
      header1.BoxSize= All.BoxSize;
      header1.Omega0= All.Omega0;
      header1.OmegaLambda= All.OmegaLambda;
      header1.HubbleParam= All.HubbleParam;
      
      dummy=sizeof(header1);
      my_fwrite(&dummy, sizeof(dummy), 1, fd);
      my_fwrite(&header1, sizeof(header1), 1, fd);
      my_fwrite(&dummy, sizeof(dummy), 1, fd);
    }

  for(posvel=0; posvel<=12; posvel++)  
    {
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

      if(All.CoolingOn==0 && (posvel==6 || posvel==7))
	continue;

      if(All.StarformationOn==0 && (posvel==8 || posvel==9 || posvel==11 || posvel==12))
	continue;

      if(All.MultiPhaseModelOn==0 && posvel==9)
	continue;

#ifndef STELLARAGE
      if(posvel==12)
	continue;
#endif

      if(ThisTask==masterTask)
	{
	  switch(posvel)
	    {
	    case 0: case 1:
	      dummy= sizeof(float)*3*ntot; 
	      break;
	    case 2:
	      dummy= sizeof(int4byte)*ntot; 
	      break;
	    case 3:
	      dummy= sizeof(float)*ntot_withmasses; 
	      break;
	    case 4: case 5: case 6: case 7: case 8: case 9: case 10: case 11:
	      dummy= sizeof(float)*ntot_type[0];   
	      break;
	    case 12: 
	      dummy= sizeof(float)*(ntot_type[0]+ntot_type[4]);  
	      break;
	    }
	  if(dummy>0)
	    my_fwrite(&dummy, sizeof(dummy), 1, fd);
	}


      for(type=0;type<5;type++)
	{
	  for(i=1,pc=0;i<=NumPart;i++)
	    {
	      if(P[i].Type == type)
		{
		  if(posvel==0)
		    {
		      for(k=0; k<3; k++) 
			block[pc][k] = P[i].PosPred[k]; 
#ifdef PERIODIC		      
		      for(k=0; k<3; k++) 
			{
			  while(block[pc][k] < 0)
			    block[pc][k]+= All.BoxSize;
			  while(block[pc][k] > All.BoxSize)
			    block[pc][k]-= All.BoxSize;
			}
#endif
		    }
		  if(posvel==1)
		    {
		      for(k=0;k<3;k++) 
			block[pc][k] = P[i].VelPred[k]; 
		    } 
		  if(posvel==2)
		    blockid[pc]  = P[i].ID;

		  if(posvel==3)
		    blockmass[pc]= P[i].Mass; 

		  if(posvel==4 && type==0)
		    blocku[pc]= SphP[i].EgySpecPred; 
		    
		  if(posvel==5 && type==0)
		    blocku[pc]= SphP[i].DensityPred; 

#ifdef COOLING
		  if(posvel==6 && type==0)
		    blocku[pc]= SphP[i].Ne;

		  if(posvel==7 && type==0)
		    {
		      ne= SphP[i].Ne;

		      AbundanceRatios(SphP[i].EgySpecPred, SphP[i].DensityPred*a3inv,
				      &ne, &nh0);

		      blocku[pc]= nh0; /* neutral hydrogen fraction */
		    }
#endif

#ifdef SFR
		  if(posvel==8 && type==0)
		    {
		      blocku[pc]=  SphP[i].FormedStellarMass;
		    }
#ifdef CLOUDS
		  if(posvel==9 && type==0)
		    {
		      blocku[pc]=  SphP[i].CloudMass; 
		    }
#endif

#endif
		  if(posvel==10 && type==0)
		    blocku[pc]= SphP[i].Hsml; 
#ifdef SFR
		  if(posvel==11 && type==0)
		    blocku[pc]= get_starformation_rate(i);
#endif

#ifdef STELLARAGE
		  if(posvel==12 && (type==0 || type==4))
		    {
		      if(type==0)
			{
			  if(SphP[i].FormedStellarMass>0)
			    blocku[pc]=  P[i].MeanStellarAge/SphP[i].FormedStellarMass; 
			  else
			    blocku[pc]= 0;
			}
		      if(type==4)
			{
			  if(P[i].Mass>0)
			    blocku[pc]=  P[i].MeanStellarAge/P[i].Mass; 
			  else
			    blocku[pc]= 0;
			}
		    }
#endif
		  pc++;
		}
	      if(pc==BLOCKSIZE)
		{
		  if(ThisTask==masterTask)
		    {
		      if(posvel<=1)
			my_fwrite(&block[0][0],sizeof(float), 3*pc ,fd); 
		      else
			{
			  if(posvel==2)
			    my_fwrite(&blockid[0],sizeof(int4byte), pc ,fd); 
			  else
			    {
			      if(posvel==3)
				{
				  if(All.MassTable[type]==0)
				    my_fwrite(&blockmass[0],sizeof(float), pc ,fd); 
				}
			      else
				{ 
				  if(posvel>=4 && posvel<=11)
				    {
				      if(type==0)  
					my_fwrite(&blocku[0], sizeof(float), pc ,fd); 
				    }
				  else  /* 12 */
				    {
				      if(type==0 || type==4)  
					{
					  if(posvel==12)
					    my_fwrite(&blocku[0], sizeof(float), pc ,fd); 
					}
				    }
				}
			    }
			}
		    }
		  else
		    {
		      MPI_Send(&pc, 1, MPI_INT, masterTask , TAG_N, MPI_COMM_WORLD);
		      MPI_Recv(&flag, 1, MPI_INT, masterTask , TAG_FLAG, MPI_COMM_WORLD,&status); 
		      /* the receive statements protects against overflows of the send buffer */
		      if(posvel<=1)
			MPI_Send(&block[0][0], 3*pc, MPI_FLOAT, masterTask , TAG_ANY, MPI_COMM_WORLD);
		      else
			{
			  if(posvel==2)
			    MPI_Send(&blockid[0], 4*pc, MPI_BYTE, masterTask , TAG_ANY, MPI_COMM_WORLD);
			  else
			    {
			      if(posvel==3)
				MPI_Send(&blockmass[0], pc, MPI_FLOAT, masterTask , TAG_ANY, MPI_COMM_WORLD);
			      else
				{
				  if(posvel>=4 && posvel<=11)
				    {
				      if(type==0)
					MPI_Send(&blocku[0],    pc, MPI_FLOAT, masterTask , TAG_ANY, MPI_COMM_WORLD);
				    }
				  else /* 12 */
				    {
				      if(type==0 || type==4)
					{
					  if(posvel==12)
					    MPI_Send(&blocku[0],    pc, MPI_FLOAT, masterTask , TAG_ANY, MPI_COMM_WORLD);
					}
				    }
				}

			    }
			}
		    }
		  pc=0;
		}
	    }
	  
	  if(pc>0)
	    {
	      if(ThisTask==masterTask)
		{
		  if(posvel<=1)
		    my_fwrite(&block[0][0],sizeof(float), 3*pc ,fd); 
		  else
		    {
		      if(posvel==2)
			my_fwrite(&blockid[0],sizeof(int4byte), pc ,fd); 
		      else
			{
			  if(posvel==3)
			    {
			      if(All.MassTable[type]==0)
				my_fwrite(&blockmass[0],sizeof(float), pc ,fd); 
			    }
			  else
			    {
			      if(posvel>=4 && posvel<=11)
				{
				  if(type==0)
				    my_fwrite(&blocku[0],sizeof(float), pc ,fd); 
				}
			      else
				{
				  if(type==0 || type==4)
				    {
				      if(posvel==12)
					my_fwrite(&blocku[0],sizeof(float), pc ,fd); 
				    }
				}
			    }
			}
		    }
		}
	      else
		{
		  MPI_Send(&pc, 1, MPI_INT, masterTask , TAG_N, MPI_COMM_WORLD);
		  MPI_Recv(&flag, 1, MPI_INT, masterTask , TAG_FLAG, MPI_COMM_WORLD,&status); 
		  /* the receive statements protects against overflows of the send buffer */
		  if(posvel<=1)
		    MPI_Send(&block[0][0], 3*pc, MPI_FLOAT, masterTask , TAG_ANY, MPI_COMM_WORLD);
		  else
		    {
		      if(posvel==2)
			MPI_Send(&blockid[0], 4*pc, MPI_BYTE, masterTask , TAG_ANY, MPI_COMM_WORLD);
		      else
			{
			  if(posvel==3)
			    MPI_Send(&blockmass[0], pc, MPI_FLOAT, masterTask , TAG_ANY, MPI_COMM_WORLD);
			  else
			    {
			      if(posvel>=4 && posvel<=11)
				{
				  if(type==0)
				    MPI_Send(&blocku[0], pc, MPI_FLOAT, masterTask , TAG_ANY, MPI_COMM_WORLD);
				}
			      else
				{
				  if(type==0 || type==4)
				    {
				      if(posvel==12)
					MPI_Send(&blocku[0], pc, MPI_FLOAT, masterTask , TAG_ANY, MPI_COMM_WORLD);
				    }
				}
			    }
					    
			}
		    }
		}
	      pc=0;
	    }

	  if(ThisTask!=masterTask)
	    MPI_Send(&pc, 1, MPI_INT, masterTask , TAG_N, MPI_COMM_WORLD);
	

	  
	  if(ThisTask==masterTask)
	    {
	      for(i=1;i<nprocgroup;i++)
		{
		  do
		    {
		      MPI_Recv(&pc, 1, MPI_INT, masterTask+i , TAG_N, MPI_COMM_WORLD,&status);
		      if(pc>0)
			{
			  MPI_Send(&flag, 1, MPI_INT, masterTask+i , TAG_FLAG, MPI_COMM_WORLD);
			  if(posvel<=1)
			    {
			      MPI_Recv(&block[0][0], 3*pc, MPI_FLOAT, masterTask+i , TAG_ANY, MPI_COMM_WORLD,&status);
			      my_fwrite(&block[0][0],sizeof(float), 3*pc ,fd); 
			    }
			  else
			    {
			      if(posvel==2)
				{
				  MPI_Recv(&blockid[0], 4*pc, MPI_BYTE, masterTask+i , TAG_ANY, MPI_COMM_WORLD,&status);
				  my_fwrite(&blockid[0],sizeof(int4byte), pc ,fd); 
				}
			      else
				{
				  if(posvel==3)
				    {
				      MPI_Recv(&blockmass[0], pc, MPI_FLOAT, masterTask+i , TAG_ANY, MPI_COMM_WORLD,&status);
				      if(All.MassTable[type]==0)
					my_fwrite(&blockmass[0],sizeof(int4byte), pc ,fd); 
				    }
				  else
				    {
				      if(posvel>=4 && posvel<=11)
					{
					  if(type==0)
					    {
					      MPI_Recv(&blocku[0], pc, MPI_FLOAT, masterTask+i , TAG_ANY, MPI_COMM_WORLD,&status);
					      my_fwrite(&blocku[0],sizeof(int4byte), pc ,fd); 
					    }
					}
				      else  /* 12 */
					{
					  if(type==0 || type==4)
					    {
					      if(posvel==12)
						{
						  MPI_Recv(&blocku[0], pc, MPI_FLOAT, masterTask+i , TAG_ANY, MPI_COMM_WORLD,&status);
						  my_fwrite(&blocku[0],sizeof(int4byte), pc ,fd); 
						}
					    }
					}
				    }
				}
			    }
			}
		    }
		  while(pc>0);
		}
	    }
	  
	}  
	  
      if(ThisTask==masterTask)
	{
	  if(dummy>0)
	    my_fwrite(&dummy, sizeof(dummy), 1, fd);
	}
    }

  
  if(ThisTask==masterTask)
    {
      fclose(fd); 
    }

#undef BLOCKSIZE 
}

 
/* This catches I/O errors occuring for my_fwrite(). In this case we better stop.
 */
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
  size_t nwritten;

  if((nwritten=fwrite(ptr, size, nmemb, stream))!=nmemb)
    {
      printf("I/O error (fwrite) on task=%d has occured.\n", ThisTask);
      fflush(stdout);
      endrun(777);
    }
  return nwritten;
}


/* This catches I/O errors occuring for fread(). In this case we better stop.
 */
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
  size_t nread;

  if((nread=fread(ptr, size, nmemb, stream))!=nmemb)
    {
      printf("I/O error (fread) on task=%d has occured.\n", ThisTask);
      fflush(stdout);
      endrun(778);
    }
  return nread;
}









