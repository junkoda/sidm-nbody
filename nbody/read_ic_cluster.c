#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"


/* This is a special-purpose read-in routine for initial conditions
 * produced with Bepi Tormen's initial conditions generator ZIC.
 * Using this routine requires that you know pretty well what you're
 * doing...
 * Be aware that there are unit conversion factor `massfac', `posfac',
 * and `velfac', that have to be set appropriately.
 * Also note that there is a boundary for the intermediate resolution
 * zone set by hand below (to a value of 24000.0 in this example).
 */
void read_ic_cluster(char *fname)
{
#define BLOCKSIZE 10000
  FILE   *fd=0;
  int    i,j,k,n,files;
  double sqr_a;
  float  *fp,dummy[BLOCKSIZE][3];
  int    pc,id,pc_here;
  int    type,pr,left,groupsize;
  MPI_Status status;
  int4byte   npart;
  float  a0;
  char   buf[100];
  int    n_for_this_task, n_in_file;
  double massfac,posfac,velfac;
  double r2;
  int    blocks,subblock;
  int4byte nhr,nlr;
  int    counttype3,counttype2;
  int    nhr_blocks, nlr_blocks;
  float  pmhr;
#ifdef T3E
  short int dummy4byte;   /* Note: int has 8 Bytes on the T3E ! */
#else
  int dummy4byte;
#endif

#define SKIP my_fread(&dummy4byte, sizeof(dummy4byte), 1, fd);

  massfac= 0.3*3*0.1*0.1/ (8*PI*All.G) * pow(141300.0/760, 3);
  posfac=  141300.0;
  velfac=  14130.0;


  /* Below, Bepi's new format is assumed !!!!!! */
  /* for the old one, the HR particle mass 'pmhr' has to be set
   * by hand! 
   */

  if(ThisTask==0) 
    {
      for(i=0;i<5;i++)
	{
	  All.MassTable[i]= 0;
	}


      if(!(fd=fopen(fname,"r")))
	{
	  printf("can't open file `%s'.\n",fname);
	  endrun(123);
	}

      fprintf(stdout,"READING FILE '%s'....\n",fname); fflush(stdout);
      
      SKIP;
      my_fread(&nhr,sizeof(int4byte),1,fd);
      my_fread(&nlr,sizeof(int4byte),1,fd);
      my_fread(&a0,sizeof(float),1,fd);
      if(dummy4byte==16)
	my_fread(&pmhr,sizeof(float),1,fd);
      else
	{
	  pmhr=  1.0000;    /* here set by hand, if necessary */
	}
      SKIP;


      All.MassTable[1]=  pmhr * massfac;  /* high-res particles */


      printf("All.MassTable[1]=%g\n", All.MassTable[1]);
      
      printf("file contains %d HR and %d LR particles.\n",
	      nhr,nlr);
      
      All.TotN_gas  = 0; 
      All.TotN_halo = nhr;
      All.TotN_disk = nlr;
      All.TotN_bulge= 0;
      All.TotN_stars= 0;

     
      printf("\nN_sph: %d\nN_halo: %d\nN_disk: %d\nN_bulge: %d\nN_stars: %d\n",
	     All.TotN_gas, All.TotN_halo, All.TotN_disk, All.TotN_bulge, All.TotN_stars);
      
      
      All.TotNumPart =    All.TotN_gas  + All.TotN_halo 
	+ All.TotN_disk + All.TotN_bulge + All.TotN_stars; 

      fclose(fd);
    }


  MPI_Bcast(&All, sizeof(struct global_data_all_processes), MPI_BYTE, 0, MPI_COMM_WORLD);

  All.MaxPart=    All.PartAllocFactor * (All.TotNumPart/NTask);    /* sets the maximum number of particles that may 
                                 	                              reside on a processor */

  allocate_memory();


  pc=1;
  id=1;


  NumPart = 0;



  for(files=0;files<=0;files++)   /* only one file here */
    {
      if(ThisTask==0)
	{
	  sprintf(buf,"%s",fname);
	  
	  if(!(fd=fopen(buf,"r")))
	    {
	      printf("can't open file `%s'.\n",buf);
	      endrun(123);
	    }

	  SKIP;
	  my_fread(&nhr,sizeof(int4byte),1,fd);
	  my_fread(&nlr,sizeof(int4byte),1,fd);
	  my_fread(&a0,sizeof(float),1,fd);
	  if(dummy4byte==16)
	    my_fread(&pmhr,sizeof(float),1,fd);
	  SKIP;

	  nhr_blocks=nhr/1000000 + 1;
	  nlr_blocks=nlr/1000000 + 1;
       	}


      MPI_Bcast(&nhr_blocks, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&nlr_blocks, 1, MPI_INT, 0, MPI_COMM_WORLD);
     

      for(blocks=0; blocks<(nhr_blocks+nlr_blocks); blocks++)
	{
	  if(ThisTask==0)
	    {
	      SKIP;
	      my_fread(&npart,sizeof(int4byte),1,fd);
	      SKIP;
	      n_in_file=npart;

	      if(blocks<nhr_blocks)
		type=1;
	      else
		type=2;
	    }

	  MPI_Bcast(&n_in_file, 1, MPI_INT, 0, MPI_COMM_WORLD);    /* receive type of particle and total number */
	  MPI_Bcast(&type,      1, MPI_INT, 0, MPI_COMM_WORLD);

	  n_for_this_task = n_in_file / NTask;
	  if(ThisTask <  (n_in_file % NTask))
	    n_for_this_task++;



	  for(subblock=0; subblock<3; subblock++)
	    {
	      if(type==1 && subblock==2)
		continue;  /* HR part's have no mass array */

	      if(ThisTask==0)
		{
		  SKIP;
		}

	      pc_here=pc;

	      for(pr=0;pr<NTask;pr++) /* go through all processes, note: pr is the receiving process */
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
			      if(subblock<2)
				my_fread(&dummy[0][0],sizeof(float), 3*groupsize ,fd);
			      else
				my_fread(&dummy[0][0],sizeof(float),  groupsize ,fd);
			    }
		      
			  if(ThisTask==0 && pr!=0)
			    {
			      if(subblock<2)
				MPI_Send(&dummy[0][0], 3*groupsize, MPI_FLOAT, pr, TAG_ANY, MPI_COMM_WORLD);
			      else
				MPI_Send(&dummy[0][0], groupsize, MPI_FLOAT, pr, TAG_ANY, MPI_COMM_WORLD);
			    }
		      
			  if(ThisTask!=0 && pr!=0)
			    {
			      if(subblock<2)
				MPI_Recv(&dummy[0][0], 3*groupsize, MPI_FLOAT, 0 , TAG_ANY, MPI_COMM_WORLD, &status);
			      else
				MPI_Recv(&dummy[0][0], groupsize, MPI_FLOAT, 0 , TAG_ANY, MPI_COMM_WORLD, &status);
			    }
		      
			  if(ThisTask==pr)
			    {
			      for(i=0,fp=&dummy[0][0];i<groupsize;i++)
				{
				  if(subblock==0)
				    {
				      P[pc_here].Type=type;
				      P[pc_here].ID  =id;    /* now set ID */
				      for(k=0;k<3;k++)
					P[pc_here].Pos[k]=dummy[i][k]; 
				      pc_here++;
				      id++;
				    }

				  if(subblock==1)
				    {
				      for(k=0;k<3;k++)
					P[pc_here].Vel[k]=dummy[i][k]; 
				      pc_here++;
				    }
				  if(subblock==2)
				    {
				      P[pc_here].Mass=fp[i];  
				      pc_here++;
				    }
				}

			    }		      
			  
			  left-=groupsize;
			}
		    }	      
		  
			      
	      
		  MPI_Barrier(MPI_COMM_WORLD);
	  
		  MPI_Bcast(&id, 1, MPI_INT, pr, MPI_COMM_WORLD);
		}

	      if(ThisTask==0)
		SKIP;
	    }

	  pc += n_for_this_task;
	  NumPart  += n_for_this_task;
				  
	}
  
    
      if(ThisTask==0)
	{
	  fclose(fd);
	}
      
    }
  

  MPI_Barrier(MPI_COMM_WORLD);
  if(ThisTask==0)
    {
      fprintf(stdout,"\nreading done.\n\n");
      fflush(stdout);
    }








  /* now convert the units */

  sqr_a=sqrt(All.Time);
      
  counttype2=counttype3=0;

  for(i=1;i<=NumPart;i++)
    {
      for(j=0;j<3;j++)
	{
	  P[i].Pos[j] = P[i].Pos[j] * posfac; /* here in units of kpc/h */

	  P[i].Vel[j] = P[i].Vel[j] * velfac; /* comoving velocity xdot on km/sec */

	  P[i].Vel[j] *= sqr_a;  /* transform to velocity variable u */
	}

      if(P[i].Type==1)
	{
	  P[i].Mass=All.MassTable[1];
	}
      else
	{
	  P[i].Mass *= massfac;
	  
	  r2=P[i].Pos[0]*P[i].Pos[0] + P[i].Pos[1]*P[i].Pos[1] + P[i].Pos[2]*P[i].Pos[2];
	  
	  if(sqrt(r2) > 24000.0)   /* boundary of inner LR zone */ 
	    {                     
	      P[i].Type=3;
	      counttype3++;
	    }
	  else
	    counttype2++;
	    
	}
    }


  printf("Task: %d has %d particles, and %d of type 3\n",ThisTask, NumPart, counttype3);

  MPI_Reduce(&counttype2, &All.TotN_disk, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&counttype3, &All.TotN_bulge, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Bcast(&All.TotN_disk, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&All.TotN_bulge, 1, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Reduce(&NumPart, &All.TotNumPart, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Bcast(&All.TotNumPart, 1, MPI_INT, 0, MPI_COMM_WORLD);

  fflush(stdout);

  if(ThisTask==0)
    {
      fprintf(stdout,"particles loaded: %d \n\n",All.TotNumPart );
      fprintf(stdout,"particles of type 2: %d \n",All.TotN_disk);
      fprintf(stdout,"particles of type 3: %d  (together %d)  \n\n",All.TotN_bulge,All.TotN_bulge+All.TotN_disk);

      fprintf(stdout,"\n");
      fprintf(stdout,"Collisionless particles   :  %d\n", All.TotNumPart-All.TotN_gas);
      fprintf(stdout,"Baryonic particles        :  %d\n", All.TotN_gas);
      fprintf(stdout,"                             ---------\n");
      fprintf(stdout,"Total number of particles :  %d\n\n", All.TotNumPart);
    }

#undef BLOCKSIZE 
}






