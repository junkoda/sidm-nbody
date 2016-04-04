#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/file.h>
#include <unistd.h>

#include "allvars.h"
#include "proto.h"


static int fd;

static int modus;   /* modus>0 read, modus==0 write */

static void in(int *x);
static void byten(void *x,int n);



/* This function reads or writes the restart files.
 * Each processor writes its own restart file, with the
 * I/O being done in parallel. To avoid congestion of the disks
 * you can tell the program to restrict the number of files
 * that are simultaneously written to NumFilesWrittenInParallel.
 * On T3E systems, you can also enable preallocation on
 * different partitions to ensure high I/O bandwidth.
 *
 * If modus>0  the restart()-routine reads, 
 * if modus==0 it writes a restart file. 
 */
void restart(int mod)
{
  char buf[200],buf_bak[200],buf_mv[500];
  double save_PartAllocFactor;
  MPI_Status status;
  int nprocgroup, masterTask, groupTask, otherTask, dummy;
  struct global_data_all_processes all_local;


  sprintf(buf,"%s%s.%d",All.OutputDir,All.RestartFile,ThisTask); 
  sprintf(buf_bak,"%s%s.%d.bak",All.OutputDir,All.RestartFile,ThisTask); 
  sprintf(buf_mv,"mv %s %s", buf, buf_bak);

  modus=mod; 

  nprocgroup=NTask/All.NumFilesWrittenInParallel;
  
  if((NTask % nprocgroup))
    {
      printf("Fatal error.\nNumber of processors must be a multiple of `NumFilesWrittenInParallel'.\n");
      endrun(213);
    }

  masterTask=(ThisTask/nprocgroup)*nprocgroup;

  for(groupTask=0; groupTask< nprocgroup; groupTask++)
    {
      if(ThisTask == (masterTask + groupTask))  /* ok, it's this processor's turn */
	{
	  if(modus)  
	    { 
	      if((fd=open(buf,O_RDONLY))<0) 
		{ 
		  fprintf(stdout,"Restart file '%s' not found.\n",buf); 
		  endrun(7870); 
		} 
	    }
	  else
	    {
	      //system(buf_mv); jk Sep12 GadgetSidm26 and later. lonestar have a problem with system();
               /* move old restart files to .bak files */
	      
	      if((fd=open(buf, O_WRONLY|O_CREAT|O_TRUNC, S_IRUSR|S_IWUSR))<0) 
		{ 
		  fprintf(stdout,"Restart file '%s' cannot be opened.\n",buf); 
		  endrun(7878); 
		}
	    }

	  
	  save_PartAllocFactor= All.PartAllocFactor;
	      
	  byten(&All,sizeof(struct global_data_all_processes));   /* common data  */
	      
	  if(modus) /* read */
	    {
	      All.PartAllocFactor= save_PartAllocFactor;
	      All.MaxPart =   All.PartAllocFactor * (All.TotNumPart/NTask); 
	      All.MaxPartSph= All.PartAllocFactor * (All.TotN_gas/NTask);   
	    }

	  if(modus>0 && groupTask==0)  /* read */
	    {
	      all_local= All;
	      MPI_Bcast( &All, sizeof(struct global_data_all_processes), MPI_BYTE, 0, MPI_COMM_WORLD);

	      if(all_local.Time != All.Time)
		{
		  printf("The restart file on task=%d is not consistent with the one on task=0\n", ThisTask);
		  fflush(stdout);
		  endrun(16);
		}

	      allocate_memory();
	    }


	  in(&NumPart); 

	  if(NumPart>All.MaxPart)
	    {
	      printf("it seems you have reduced(!) 'PartAllocFactor' below the value of %g needed to load the restart file.\n",
		     ((double)(NumPart/(All.TotNumPart/NTask))));
	      printf("fatal error\n");
	      endrun(22);
	    }
	  
	  byten(&P[1], NumPart*sizeof(struct particle_data));  /* Particle data  */

	  in(&N_gas);

          if(N_gas>0)
	    byten(&SphP[1], N_gas*sizeof(struct sph_particle_data));  /* Sph-Particle data  */


	  close(fd);


	  /* send a continuation message to the other group members */
	  
	  for(otherTask=0; otherTask< nprocgroup; otherTask++)
	    {
	      if(otherTask != groupTask)
		MPI_Ssend(&otherTask, 1, MPI_INT, masterTask + otherTask,  masterTask+ groupTask, MPI_COMM_WORLD);
	    }
	}
      else  /* wait inside the group */
	{
	  if(modus>0 && groupTask==0)  /* read */
	    {
	      MPI_Bcast( &All, sizeof(struct global_data_all_processes), MPI_BYTE, 0, MPI_COMM_WORLD);
	      allocate_memory();
	    }

	  MPI_Recv(&dummy,   1, MPI_INT, masterTask+ groupTask, masterTask+ groupTask, MPI_COMM_WORLD, &status);
	}
    }
}



/* reads/writes n bytes 
 */
void byten(void *x,int n)
{
  if(modus)
    {
      if(read(fd,x,n*sizeof(char)) != n*sizeof(char))
	{
	  printf("read error on task=%d (restart file appears to be truncated)\n", ThisTask);
	  fflush(stdout);
	  endrun(7);
	}
    }
  else
    {
      if(write(fd,x,n*sizeof(char)) != n*sizeof(char))
	{
	  printf("write error on task=%d upon writing restart file\n", ThisTask);
	  fflush(stdout);
	  endrun(8);
	}
    }
}


/* reads/writes one int 
 */
void in(int *x)
{
  if(modus)
    {
      if(read(fd, x,1*sizeof(int)) != sizeof(int))
	{
	  printf("read error on task=%d (restart file appears to be truncated)\n", ThisTask);
	  fflush(stdout);
	  endrun(7);
	}
    }
  else
    {
      if(write(fd, x,1*sizeof(int)) != sizeof(int))
	{
	  printf("write error on task=%d upon writing restart file\n", ThisTask);
  	  fflush(stdout);
	  endrun(8);
	}
    }
}



















