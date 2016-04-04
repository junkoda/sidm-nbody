#ifdef PERIODIC

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"

#define EN  64

#define ED (EN/2)

static float fcorrx[ED+1][ED+1][ED+1];
static float fcorry[ED+1][ED+1][ED+1];
static float fcorrz[ED+1][ED+1][ED+1];
static float potcorr[ED+1][ED+1][ED+1];
static double fac_intp;


/* set-up tables with the correction force and the correction
 * potential due to the periodic images of a point mass
 * located at the origin. These corrections are obtained by
 * Ewald summation. (see Hernquist, Bouchet, Suto, ApJS, 1991, 75, 231)
 * These corrections fields are used to obtain
 * the full periodic force during the tree walk.
 *
 * The fields are stored on disk once they are computed. If a corresponding
 * file is found, they are loaded from disk to speed up the initialization.
 * The Ewald summation is done in parallel, i.e. the processors share the
 * work to set-up the tables evenly.
 */
void ewald_init(void)
{
  int i, j, k, beg, len, size, n, task, count;
  double x[3], force[3];
  char   buf[200];
  FILE  *fd;

  if(ThisTask==0)
    {
      printf("initialize Ewald correction...\n");
      fflush(stdout);
    }

  sprintf(buf, "ewald_table_%d.dat", EN);

  if((fd=fopen(buf, "r")))
    {
      if(ThisTask==0)
	{
	  printf("\nreading Ewald tables from file `%s'\n", buf);
	  fflush(stdout);
	}

      my_fread(&fcorrx[0][0][0], sizeof(float), (ED+1)*(ED+1)*(ED+1), fd);
      my_fread(&fcorry[0][0][0], sizeof(float), (ED+1)*(ED+1)*(ED+1), fd);
      my_fread(&fcorrz[0][0][0], sizeof(float), (ED+1)*(ED+1)*(ED+1), fd);
      my_fread(&potcorr[0][0][0], sizeof(float), (ED+1)*(ED+1)*(ED+1), fd);

      fclose(fd);
    }
  else
    {
      if(ThisTask==0)
	{
	  printf("\nNo Ewald tables in file `%s' found.\nRecomputing them...\n", buf);
	  fflush(stdout);
	}

      /* ok, let's recomoute things. Actually, we can do that in parallel. */

      size= (ED+1)*(ED+1)*(ED+1)/NTask;


      beg=  ThisTask*size;
      len=  size;
      if(ThisTask == (NTask-1))
	len= (ED+1)*(ED+1)*(ED+1)-beg;
      
      for(i=0, count=0; i<=ED; i++)
	for(j=0; j<=ED; j++)
	  for(k=0; k<=ED; k++)
	    {
	      n=(i*(ED+1) + j)*(ED+1) + k;
	      if(n>=beg && n<(beg+len))
		{
		  if(ThisTask==0)
		    {
		      if((count % (len/20)) == 0)
			printf("%4.1f percent done\n", count/(len/100.0) );
		    }

		  x[0]= ((double)i)/EN;
		  x[1]= ((double)j)/EN;
		  x[2]= ((double)k)/EN;
		  
		  ewald_force(x, force);
		  
		  if(i+j+k==0)
		    potcorr[i][j][k]= 2.8372975;
		  else
		    potcorr[i][j][k]= ewald_psi(x);


		  fcorrx[i][j][k]= force[0];
		  fcorry[i][j][k]= force[1];
		  fcorrz[i][j][k]= force[2];

		  count++;
		}
	    }

      for(task=0; task<NTask; task++)
	{
	  beg=  task*size;
	  len=  size;
	  if(task == (NTask-1))
	    len= (ED+1)*(ED+1)*(ED+1)- beg;

	  MPI_Bcast(&fcorrx[0][0][beg],  len, MPI_FLOAT, task, MPI_COMM_WORLD);
	  MPI_Bcast(&fcorry[0][0][beg],  len, MPI_FLOAT, task, MPI_COMM_WORLD);
	  MPI_Bcast(&fcorrz[0][0][beg],  len, MPI_FLOAT, task, MPI_COMM_WORLD);
	  MPI_Bcast(&potcorr[0][0][beg], len, MPI_FLOAT, task, MPI_COMM_WORLD);
	}

      if(ThisTask==0)
	{ 
	  printf("\nwriting Ewald tables to file `%s'\n", buf);
	  fflush(stdout);

	  if((fd=fopen(buf, "w")))
	    {
	      my_fwrite(&fcorrx[0][0][0], sizeof(float), (ED+1)*(ED+1)*(ED+1), fd);
	      my_fwrite(&fcorry[0][0][0], sizeof(float), (ED+1)*(ED+1)*(ED+1), fd);
	      my_fwrite(&fcorrz[0][0][0], sizeof(float), (ED+1)*(ED+1)*(ED+1), fd);
	      my_fwrite(&potcorr[0][0][0], sizeof(float), (ED+1)*(ED+1)*(ED+1), fd);
	      fclose(fd);
	    }
	}
    }

  fac_intp= EN/All.BoxSize;

  for(i=0;i<=ED;i++)
    for(j=0;j<=ED;j++)
      for(k=0;k<=ED;k++)
	{
	  potcorr[i][j][k]/= All.BoxSize;
	  fcorrx[i][j][k]/= All.BoxSize*All.BoxSize;
	  fcorry[i][j][k]/= All.BoxSize*All.BoxSize;
	  fcorrz[i][j][k]/= All.BoxSize*All.BoxSize;
	}

  if(ThisTask==0)
    {
      printf("initialization of periodic boundaries finished.\n");
      fflush(stdout);
    }
}

/* This function computes the correction force due
 * to the infinte number of periodic particle/node images.
 * This correction force is usually computed by the
 * Ewald summation method, but it can also be derived
 * by an FFT techniqe.
 */
// originally INLINE_FUNC
void ewald_corr(double dx, double dy, double dz, double *fper)
{
  int signx, signy, signz;
  int i,j,k;
  double u,v,w;
  double f1,f2,f3,f4,f5,f6,f7,f8;
  
  if(dx<0)
    {
      dx=-dx; signx=+1;
    }
  else
    signx=-1;

  if(dy<0)
    {
      dy=-dy; signy=+1;
    }
  else
    signy=-1;

  if(dz<0)
    {
      dz=-dz; signz=+1;
    }
  else
    signz=-1;
 
  u=dx*fac_intp; i=(int)u; if(i>=ED) i=ED-1; u-=i;
  v=dy*fac_intp; j=(int)v; if(j>=ED) j=ED-1; v-=j;
  w=dz*fac_intp; k=(int)w; if(k>=ED) k=ED-1; w-=k;

  f1= (1-u)*(1-v)*(1-w);
  f2= (1-u)*(1-v)*(w);
  f3= (1-u)*(v)*(1-w);
  f4= (1-u)*(v)*(w);
  f5= (u)*(1-v)*(1-w);
  f6= (u)*(1-v)*(w);
  f7= (u)*(v)*(1-w);
  f8= (u)*(v)*(w);

  fper[0]= signx * (fcorrx[i][j][k]*f1+
		    fcorrx[i][j][k+1]*f2+
		    fcorrx[i][j+1][k]*f3+
		    fcorrx[i][j+1][k+1]*f4+
		    fcorrx[i+1][j][k]*f5+
		    fcorrx[i+1][j][k+1]*f6+
		    fcorrx[i+1][j+1][k]*f7+
		    fcorrx[i+1][j+1][k+1]*f8);

  fper[1]= signy * (fcorry[i][j][k]*f1+
		    fcorry[i][j][k+1]*f2+
		    fcorry[i][j+1][k]*f3+
		    fcorry[i][j+1][k+1]*f4+
		    fcorry[i+1][j][k]*f5+
		    fcorry[i+1][j][k+1]*f6+
		    fcorry[i+1][j+1][k]*f7+
		    fcorry[i+1][j+1][k+1]*f8);

  fper[2]= signz * (fcorrz[i][j][k]*f1+
		    fcorrz[i][j][k+1]*f2+
		    fcorrz[i][j+1][k]*f3+
		    fcorrz[i][j+1][k+1]*f4+
		    fcorrz[i+1][j][k]*f5+
		    fcorrz[i+1][j][k+1]*f6+
		    fcorrz[i+1][j+1][k]*f7+
		    fcorrz[i+1][j+1][k+1]*f8);
}


/* This function computes the correction to the potential due
 * to the infinte number of periodic particle/node images.
 * This correction potential is usually computed by the
 * Ewald summation method, but it can also be derived
 * by an FFT techniqe.
 */
// originally INLINE_FUNC
double ewald_pot_corr(double dx, double dy, double dz)
{
  int    i,j,k;
  double u,v,w;
  double f1,f2,f3,f4,f5,f6,f7,f8;

  if(dx<0)
    dx=-dx; 
    
  if(dy<0)
    dy=-dy; 

  if(dz<0)
    dz=-dz; 
 
  u=dx*fac_intp; i=(int)u; if(i>=ED) i=ED-1; u-=i;
  v=dy*fac_intp; j=(int)v; if(j>=ED) j=ED-1; v-=j;
  w=dz*fac_intp; k=(int)w; if(k>=ED) k=ED-1; w-=k;

  f1= (1-u)*(1-v)*(1-w);
  f2= (1-u)*(1-v)*(w);
  f3= (1-u)*(v)*(1-w);
  f4= (1-u)*(v)*(w);
  f5= (u)*(1-v)*(1-w);
  f6= (u)*(1-v)*(w);
  f7= (u)*(v)*(1-w);
  f8= (u)*(v)*(w);

  return          potcorr[i][j][k]*f1+
	          potcorr[i][j][k+1]*f2+
		  potcorr[i][j+1][k]*f3+
		  potcorr[i][j+1][k+1]*f4+
		  potcorr[i+1][j][k]*f5+
		  potcorr[i+1][j][k+1]*f6+
		  potcorr[i+1][j+1][k]*f7+
		  potcorr[i+1][j+1][k+1]*f8;
}



/* This function computes the potential correction term
 * by means of Ewald summation
 */
double ewald_psi(double x[3])
{
  double alpha, psi;
  double r, sum1, sum2, hdotx;
  double dx[3];
  int    i, n[3], h[3], h2;

  alpha = 2.0;
 
  for(n[0]=-4, sum1=0; n[0]<=4; n[0]++)
    for(n[1]=-4; n[1]<=4; n[1]++)
      for(n[2]=-4; n[2]<=4; n[2]++)
	{
	  for(i=0; i<3; i++)
	    dx[i]= x[i]-n[i];

	  r=  sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
	  sum1+= erfc(alpha*r)/r;
	}
  
  for(h[0]=-4, sum2=0; h[0]<=4; h[0]++)
    for(h[1]=-4; h[1]<=4; h[1]++)
      for(h[2]=-4; h[2]<=4; h[2]++)
	{
	  hdotx= x[0]*h[0] + x[1]*h[1] + x[2]*h[2];
	  h2=    h[0]*h[0] + h[1]*h[1] + h[2]*h[2];
	  if(h2>0)
	    sum2+= 1/(PI*h2)*exp(-PI*PI*h2/(alpha*alpha))*cos(2*PI*hdotx);
	}
  
  r=  sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

  psi= PI/(alpha*alpha) - sum1 - sum2 + 1/r;

  return psi;
}


/* This function computes the force correction term
 * by means of Ewald summation
 */
void ewald_force(double x[3], double force[3])
{
  double alpha;
  double r, r2, val, hdotx, dx[3];
  int    i, h[3], n[3], h2;

  alpha = 2.0;

  for(i=0; i<3; i++)
    force[i]=0;
  
  r2= x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
  if(r2==0)
    return;
  
  for(i=0; i<3; i++)
    force[i]+= x[i]/(r2*sqrt(r2));
  
  for(n[0]=-4; n[0]<=4; n[0]++)
    for(n[1]=-4; n[1]<=4; n[1]++)
      for(n[2]=-4; n[2]<=4; n[2]++)
	{
	  for(i=0; i<3; i++)
	    dx[i]= x[i]-n[i];
	  
	  r=  sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
	  
	  val= erfc(alpha*r) + 2*alpha*r/sqrt(PI)*exp(-alpha*alpha*r*r);
	  
	  for(i=0; i<3; i++)
	    force[i]-= dx[i]/(r*r*r)*val;
	}
  
  for(h[0]=-4; h[0]<=4; h[0]++)
    for(h[1]=-4; h[1]<=4; h[1]++)
      for(h[2]=-4; h[2]<=4; h[2]++)
	{
	  hdotx= x[0]*h[0] + x[1]*h[1] + x[2]*h[2];
	  h2=    h[0]*h[0] + h[1]*h[1] + h[2]*h[2];
	  
	  if(h2>0)
	    {
	      val= 2.0/((double)h2)*exp(-PI*PI*h2/(alpha*alpha))*
		   sin(2*PI*hdotx);
	      
	      for(i=0; i<3; i++)
		force[i]-= h[i]*val;
	    }
	}
}


#endif




