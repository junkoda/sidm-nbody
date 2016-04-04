#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"



/* Allocates communication buffers, plus the list of the interior
 * domains (for neighbour search and hydro). 
 */
void allocate_commbuffers(void) 
{
  int bytes;

  if(!(InteriorMin= malloc(bytes= 3*NTask*sizeof(float))))
    {
      printf("failed to allocate memory for `InteriorGasMin' (%d bytes).\n",bytes);
      endrun(2);
    }

  if(!(InteriorMax= malloc(bytes= 3*NTask*sizeof(float))))
    {
      printf("failed to allocate memory for `InteriorGasMax' (%d bytes).\n",bytes);
      endrun(2);
    }

  if(!(CommBuffer= malloc(bytes= All.BufferSize*1024*1024)))
    {
      printf("failed to allocate memory for `CommBuffer' (%d bytes).\n",bytes);
      endrun(2);
    }
  
  All.BunchSizeForce= (All.BufferSize*1024*1024)/(sizeof(struct gravdata_in) + sizeof(struct gravdata_out) + sizeof(struct gravdata_out) + 2*sizeof(double));

  if(All.BunchSizeForce&1)
    All.BunchSizeForce-= 1;  /* make sure that All.BunchSizeForce is even 
                                --> 8-byte alignment of GravDataResult for 64bit processors */
  
  GravDataIn= (struct gravdata_in *)CommBuffer;
  GravDataResult= (struct gravdata_out *)(GravDataIn + All.BunchSizeForce);
  GravDataPartialResult= GravDataResult + All.BunchSizeForce;

  GravDataPotential= (double *)(GravDataPartialResult + All.BunchSizeForce);
  GravDataPartialPotential= GravDataPotential + All.BunchSizeForce;

  All.BunchSizeDensity= (All.BufferSize*1024*1024)/(sizeof(struct densdata_in) + 2*sizeof(struct densdata_out));

  DensDataIn= (struct densdata_in *)CommBuffer;
  DensDataResult= (struct densdata_out *)(DensDataIn + All.BunchSizeDensity);
  DensDataPartialResult= DensDataResult +  All.BunchSizeDensity;

#ifdef VELDISP
  All.BunchSizeVelDisp= (All.BufferSize*1024*1024)/(sizeof(struct veldispdata_in) + 2*sizeof(struct veldispdata_out));

  VelDispDataIn= (struct veldispdata_in *)CommBuffer;
  VelDispDataResult= (struct veldispdata_out *)(VelDispDataIn + All.BunchSizeVelDisp);
  VelDispDataPartialResult= VelDispDataResult +  All.BunchSizeVelDisp;
#endif

#ifdef SIDM
  All.BunchSizeSidm= (All.BufferSize*1024*1024)/(sizeof(struct sidmdata_in) + 2*sizeof(struct sidmdata_out) + 
                                                 sizeof(struct sidmdata_confirm) + sizeof(int4byte));

  SidmDataIn= (struct sidmdata_in *)CommBuffer;
  SidmDataResult= (struct sidmdata_out *)(SidmDataIn + All.BunchSizeSidm);
  SidmDataPartialResult= SidmDataResult +  All.BunchSizeSidm;
  SidmDataConfirm= (struct sidmdata_confirm *)(SidmDataPartialResult + All.BunchSizeSidm);
  SidmTarget= (int4byte*)(SidmDataConfirm + All.BunchSizeSidm);
#endif

  All.BunchSizeHydro= (All.BufferSize*1024*1024)/(sizeof(struct hydrodata_in) + 2*sizeof(struct hydrodata_out));

  HydroDataIn= (struct hydrodata_in *)CommBuffer;
  HydroDataResult= (struct hydrodata_out *)(HydroDataIn + All.BunchSizeHydro);
  HydroDataPartialResult= HydroDataResult +  All.BunchSizeHydro;

  All.BunchSizeDomain= (All.BufferSize*1024*1024)/(sizeof(struct particle_data) + sizeof(struct sph_particle_data));

  DomainPartBuf=(struct particle_data *)CommBuffer;
  DomainSphBuf= (struct sph_particle_data *)(DomainPartBuf + All.BunchSizeDomain);

#ifdef SFR
  All.BunchSizeWeight= (All.BufferSize*1024*1024)/(sizeof(struct weightdata_in) + 2*sizeof(struct weightdata_out));

  WeightDataIn= (struct weightdata_in *)CommBuffer;
  WeightDataResult= (struct weightdata_out *)(WeightDataIn + All.BunchSizeWeight);
  WeightDataPartialResult= WeightDataResult +  All.BunchSizeWeight;


  All.BunchSizeDissolve= (All.BufferSize*1024*1024)/(sizeof(struct dissolvedata_in) + 2*sizeof(struct dissolvedata_out));

  DissolveDataIn= (struct dissolvedata_in *)CommBuffer;
  DissolveDataResult= (struct dissolvedata_out *)(DissolveDataIn + All.BunchSizeDissolve);
  DissolveDataPartialResult= DissolveDataResult +  All.BunchSizeDissolve;
#endif

  if(ThisTask==0)
    {
      printf("\nAllocated %d MByte communication buffer per processor.\n\n", All.BufferSize);
      printf("Communication buffer has room for %d particles in gravity computation\n", All.BunchSizeForce);
      printf("Communication buffer has room for %d particles in density computation\n", All.BunchSizeDensity);
      printf("Communication buffer has room for %d particles in hydro computation\n", All.BunchSizeHydro);
      printf("Communication buffer has room for %d particles in domain decomposition\n", All.BunchSizeDomain);
      printf("\n");
    }
}





/* This routine allocates memory for 
 * particle storage, both the collisionless and the SPH particles.
 * The memory for the ordered binary tree of the timeline
 * is also allocated.
 */
void allocate_memory(void)
{
  int bytes,bytes_tot=0;

  if(All.MaxPart>0)
    {
      if(!(P_data=malloc(bytes=All.MaxPart*sizeof(struct particle_data))))
	{
	  printf("failed to allocate memory for `P_data' (%d bytes).\n",bytes);
	  endrun(1);
	}
      bytes_tot+=bytes;


      if(!(PTimeTree=malloc(bytes=All.MaxPart*sizeof(struct timetree_data))))
	{
	  printf("failed to allocate memory for `PTimeTree' (%d bytes).\n",bytes);
	  endrun(1);
	}
      bytes_tot+=bytes;


      P= P_data-1;   /* start with offset 1 */
      PTimeTree--;

      if(ThisTask==0)
	printf("\nAllocated %g MByte for particle storage.\n\n",bytes_tot/(1024.0*1024.0));
    }

  if(All.MaxPartSph>0)
    {
      bytes_tot=0;

      if(!(SphP_data=malloc(bytes=All.MaxPartSph*sizeof(struct  sph_particle_data))))
	{
	  printf("failed to allocate memory for `SphP_data' (%d bytes).\n",bytes);
	  endrun(1);
	}
      bytes_tot+=bytes;

      SphP= SphP_data-1;   /* start with offset 1 */

      if(ThisTask==0)
	printf("Allocated %g MByte for storage of SPH data.\n\n",bytes_tot/(1024.0*1024.0));
    }
}




/* This routine frees the memory for the particle storage
 * and for the communication buffers,
 * but we don't actually call it in the code. 
 * When the program terminats, the memory will be automatically
 * freed by the operating system.
 */
void free_memory(void)
{
  if(All.MaxPart>0)
    free(P_data);

  if(All.MaxPartSph>0)
    free(SphP_data);

  free(CommBuffer);
}


