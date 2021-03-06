#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"



/* This routine computes various global properties of the particle
 * distribution and stores the result in the struct `SysState'.
 * Currently, not all the information that's computed here is 
 * actually used (e.g. momentum is not really used anywhere),
 * just the energies are written to a log-file every once in a while.
 */
void compute_global_quantities_of_system(void)
{
  int    i, j, n;
  struct state_of_system sys;

  for(n=0; n<5; n++)   
    { 
       sys.MassComp[n]= sys.EnergyKinComp[n]= sys.EnergyPotComp[n]= sys.EnergyIntComp[n]= 0; 
       
       for(j=0;j<4;j++) 
	 sys.CenterOfMassComp[n][j]= sys.MomentumComp[n][j] = sys.AngMomentumComp[n][j] = 0; 
    } 

  for(i=1; i<=NumPart; i++) 
    { 
      sys.MassComp[P[i].Type] += P[i].Mass; 

      sys.EnergyPotComp[P[i].Type] += 0.5*P[i].Mass*P[i].Potential; 
      
      if(P[i].Type==0)
	sys.EnergyIntComp[0] += P[i].Mass * SphP[i].EgySpecPred; 
	
      sys.EnergyKinComp[P[i].Type] += 0.5*P[i].Mass*(P[i].VelPred[0]*P[i].VelPred[0] +
						       P[i].VelPred[1]*P[i].VelPred[1] + 
						       P[i].VelPred[2]*P[i].VelPred[2]); 
      
      for(j=0; j<3; j++) 
	{ 
	  sys.MomentumComp[P[i].Type][j] += P[i].Mass* P[i].VelPred[j]; 
	  sys.CenterOfMassComp[P[i].Type][j] += P[i].Mass*P[i].PosPred[j]; 
	} 
      
      sys.AngMomentumComp[P[i].Type][0] += P[i].Mass* ( P[i].PosPred[1]*P[i].VelPred[2] -
							  P[i].PosPred[2]*P[i].VelPred[1] ); 
      sys.AngMomentumComp[P[i].Type][1] += P[i].Mass* ( P[i].PosPred[2]*P[i].VelPred[0] -
							  P[i].PosPred[0]*P[i].VelPred[2] ); 
      sys.AngMomentumComp[P[i].Type][2] += P[i].Mass* ( P[i].PosPred[0]*P[i].VelPred[1] -
							  P[i].PosPred[1]*P[i].VelPred[0] ); 
    } 


  /* some the stuff over all processors */
  MPI_Reduce(&sys.MassComp[0], &SysState.MassComp[0], 5, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sys.EnergyPotComp[0], &SysState.EnergyPotComp[0], 5, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sys.EnergyIntComp[0], &SysState.EnergyIntComp[0], 5, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sys.EnergyKinComp[0], &SysState.EnergyKinComp[0], 5, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sys.MomentumComp[0][0], &SysState.MomentumComp[0][0], 5*4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sys.AngMomentumComp[0][0], &SysState.AngMomentumComp[0][0], 5*4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sys.CenterOfMassComp[0][0], &SysState.CenterOfMassComp[0][0], 5*4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


  if(ThisTask==0)
    {
      for(i=0;i<5;i++) 
	  SysState.EnergyTotComp[i] = SysState.EnergyKinComp[i] + 
	                              SysState.EnergyPotComp[i] + 
                                      SysState.EnergyIntComp[i];
   
      SysState.Mass= SysState.EnergyKin= SysState.EnergyPot= SysState.EnergyInt= SysState.EnergyTot=0; 

      for(j=0; j<3; j++) 
	SysState.Momentum[j]= SysState.AngMomentum[j]= SysState.CenterOfMass[j]= 0; 
 
      for(i=0; i<5; i++) 
	{ 
	  SysState.Mass += SysState.MassComp[i]; 
	  SysState.EnergyKin += SysState.EnergyKinComp[i]; 
	  SysState.EnergyPot += SysState.EnergyPotComp[i]; 
	  SysState.EnergyInt += SysState.EnergyIntComp[i]; 
	  SysState.EnergyTot += SysState.EnergyTotComp[i]; 

	  for(j=0; j<3; j++) 
	    { 
	      SysState.Momentum[j] += SysState.MomentumComp[i][j]; 
	      SysState.AngMomentum[j] += SysState.AngMomentumComp[i][j]; 
	      SysState.CenterOfMass[j] += SysState.CenterOfMassComp[i][j]; 
	    } 
	} 
      
      for(i=0; i<5; i++) 
	for(j=0; j<3; j++) 
	  if(SysState.MassComp[i] > 0) 
	    SysState.CenterOfMassComp[i][j] /= SysState.MassComp[i]; 
      
      for(j=0; j<3; j++) 
	if(SysState.Mass > 0) 
	  SysState.CenterOfMass[j] /= SysState.Mass; 

      for(i=0; i<5; i++) 
	{ 
	  SysState.CenterOfMassComp[i][3]= SysState.MomentumComp[i][3]= SysState.AngMomentumComp[i][3]=0; 
	  for(j=0; j<3; j++) 
	    { 
	      SysState.CenterOfMassComp[i][3] += SysState.CenterOfMassComp[i][j]*SysState.CenterOfMassComp[i][j]; 
	      SysState.MomentumComp[i][3] += SysState.MomentumComp[i][j]*SysState.MomentumComp[i][j]; 
	      SysState.AngMomentumComp[i][3] += SysState.AngMomentumComp[i][j]*SysState.AngMomentumComp[i][j]; 
	    } 
	  SysState.CenterOfMassComp[i][3]= sqrt(SysState.CenterOfMassComp[i][3]); 
	  SysState.MomentumComp[i][3]= sqrt(SysState.MomentumComp[i][3]); 
	  SysState.AngMomentumComp[i][3]= sqrt(SysState.AngMomentumComp[i][3]); 
	} 
      
      SysState.CenterOfMass[3]= SysState.Momentum[3]= SysState.AngMomentum[3]= 0; 
      
      for(j=0; j<3; j++) 
	{ 
	  SysState.CenterOfMass[3] += SysState.CenterOfMass[j]*SysState.CenterOfMass[j]; 
	  SysState.Momentum[3] += SysState.Momentum[j]*SysState.Momentum[j]; 
	  SysState.AngMomentum[3] += SysState.AngMomentum[j]*SysState.AngMomentum[j]; 
	} 
      
      SysState.CenterOfMass[3]= sqrt(SysState.CenterOfMass[3]); 
      SysState.Momentum[3]= sqrt(SysState.Momentum[3]); 
      SysState.AngMomentum[3]= sqrt(SysState.AngMomentum[3]); 
    }

  /* give everyone the result, maybe the want to do something with it */
  MPI_Bcast(&SysState, sizeof(struct state_of_system), MPI_BYTE, 0,  MPI_COMM_WORLD);
}















