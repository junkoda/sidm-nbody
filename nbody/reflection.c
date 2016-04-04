// reflection.c
// spherical reflection boundary

#include "allvars.h"

#ifdef  REFLECTIONBOUNDARY
void reflect(void)
{
  const float r_ref2= (float)(All.ReflectionRadius*All.ReflectionRadius);

  float r2, r2inv2, rv;
  float x, y, z;
  int i, count;
  for(i=IndFirstUpdate, count=0; count<NumForceUpdate;
      i=P[i].ForceFlag, count++){
    
    r2= P[i].Pos[0]*P[i].Pos[0] + P[i].Pos[1]*P[i].Pos[1]
      + P[i].Pos[2]*P[i].Pos[2];
    if(r2>r_ref2){
      x = P[i].Pos[0];
      y = P[i].Pos[1];
      z = P[i].Pos[2];
      rv= x*P[i].Vel[0] + y*P[i].Vel[1] + z*P[i].Vel[2];
      if(rv > 0){
	r2inv2= 2.0f/r2;
	P[i].Vel[0] -= rv*x*r2inv2;
	P[i].Vel[1] -= rv*y*r2inv2;
	P[i].Vel[2] -= rv*z*r2inv2;
      }
    }
  }
}

#endif
