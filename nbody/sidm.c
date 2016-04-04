#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>
#include "allvars.h"
#include "proto.h"
#include "sidm.h"
#include "sidm_rand.h"

#ifdef SIDM

#ifdef VELDISP
#error VELDISP cannot be used with SIDM
#endif

/* Self interaction +
    neighbor counting */

/* version 0.1.0 */

#if (CROSS_SECTION_TYPE == 4)
static inline double norm(const double x[])
{
  return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}
    
static inline void perp(const double x[], const double y[], double out[])
{
  // out is a unit vector perpendicular to both x and y
  double oo;
  out[0]= x[1]*y[2] - x[2]*y[1];
  out[1]= x[2]*y[0] - x[0]*y[2];
  out[2]= x[0]*y[1] - x[1]*y[0];
  
  oo= sqrt(out[0]*out[0] + out[1]*out[1] + out[2]*out[2]);

  if(oo == 0.0) {
    fprintf(stderr, "x= %e %e %e\n", x[0], x[1], x[2]);
    fprintf(stderr, "y= %e %e %e\n", y[0], y[1], y[2]);
    endrun(4003);
  }

  out[0] /= oo;
  out[1] /= oo;
  out[2] /= oo;

  /*
  if(out[0] != out[0] || out[1] != out[1] || out[2] != out[2])
    endrun(4004);
  */
}
#endif


void sidm(void)
{
  int     i,ii,j,k,n, numngb;
  double  h, hinv, hinv3, r, u, wk=0.0;      
  float   *r2list;
  int     *ngblist;
  double  tstart,tend;
  int     ntot,ntotleft,npleft,nthis;     
  int     nstart,nbuffer,ncount,nchunk;
  int     timelinecounter,startcounter;  
  int     *nrecv,*noffset;
  int     level,sendTask,recvTask;
  int     place, nexport;
  int     num_collisionless;
  MPI_Status status;
  
  double  dt_h0;
  double  rand, C_Pmax, P_max, Prob;
  double  rv, rvx, rvy, rvz;
  double  rmass;
  double  nx[3];
  double  dvx, dvy, dvz;
  double  s_a, s_a_inverse;
  double  CrossSectionCo;
  int    nscat[3], nscatTot[3]; 

#if (CROSS_SECTION_TYPE == 2)
  double beta, v_dep;
#elif (CROSS_SECTION_TYPE == 3)
  double n_cross_section; // sigma = sigma0*(dv/v_scale)^n
  double v_scale;
#elif (CROSS_SECTION_TYPE == 4)
  double beta, cosO, sin22, sinO, denom;
  double rvv[3], nperp[3];
#endif

  
  
  /* file for scattering log */
#ifdef SCATTERLOG
  FILE   *fp;
  char   filename[128];
  struct scatlog ScatLog;
  ScatLog.time= All.Time;
  sprintf(filename, "sct_%03d.%d", All.SnapshotFileCount, ThisTask);
  //sprintf(filename, "sct.%d", ThisTask);
  fp = fopen(filename, "ab");
#endif

  for(k=0; k<3; k++)
    nscat[k]= nscatTot[k]= 0;
    
  num_collisionless= NumForceUpdate-NumSphUpdate;

  MPI_Allreduce(&num_collisionless, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  /* compute maximum size of bunch that may come from this task */
  if(num_collisionless==0)
    nthis=1;
  else
    nthis = ((double)(num_collisionless)*(All.BunchSizeSidm-NTask))/ntot + 1; 


  if(num_collisionless>0) {
    nchunk= num_collisionless/nthis; /* number of bunches needed */  
    if((num_collisionless % nthis)>0) nchunk+=1;
  }
  else
    nchunk= 0;


  nrecv= malloc(sizeof(int)*NTask);  /* list of particle numbers that constituate current bunch */
  noffset= malloc(sizeof(int)*NTask);  /* offsets of bunches in common list */

  ntotleft= ntot;              /* particles left for all tasks together */
  npleft= num_collisionless;   /* particles left for this task */

  nstart= IndFirstUpdate;      /* first particle for this task */
  startcounter=0;

  while(ntotleft > 0){
    if(nthis > npleft)
      nthis= npleft;

    for(i=nstart, ncount=0, nexport=0, timelinecounter=startcounter; 
        timelinecounter<NumForceUpdate && ncount<nthis; 
        i=P[i].ForceFlag, timelinecounter++) {
      if(P[i].Type>0) {
        ncount++;
        
        for(j=0; j<3; j++) {
          if(P[i].PosPred[j] < (DomainMin[P[i].Type][j] + P[i].HsmlVelDisp))
            break;
          if(P[i].PosPred[j] > (DomainMax[P[i].Type][j] - P[i].HsmlVelDisp))
            break;
        }
        
        if(j != 3) {
	  /* particle lies NOT completely inside. 
	     needs to be sent to other processors */
          nexport++;
          P[i].Type |= 8;
        }
      }
    }

    MPI_Allgather(&nexport, 1, MPI_INT, nrecv, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allreduce(&ncount, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      
    for(i=nbuffer=0; i<NTask; i++)  /* compute length of common list */
      nbuffer += nrecv[i];
      
    for(i=1, noffset[0]=0; i<NTask; i++) /* set-up offset-table */
      noffset[i]= noffset[i-1] + nrecv[i-1];

    /* fill in the own particles at the right place */
    for(i=nstart, ncount=0, nexport=0, timelinecounter=startcounter; 
        timelinecounter < NumForceUpdate && ncount < nthis; 
        i=P[i].ForceFlag, timelinecounter++) {
      if(P[i].Type > 0) {
        if(P[i].Type & 8) {
          place= noffset[ThisTask] + nexport; 
          nexport++;
        }
        else
          place= nbuffer + (ncount-nexport);
        
        for(k=0; k<3; k++) {
          SidmDataIn[place].Pos[k]= P[i].PosPred[k];
          SidmDataIn[place].Vel[k]= P[i].Vel[k];
        }
        SidmDataIn[place].Hsml= P[i].HsmlVelDisp;
        SidmDataIn[place].Type= P[i].Type & 7;
        SidmDataIn[place].Mass= P[i].Mass;
	if(P[i].dVel[0] != 0.0f)
	  SidmDataIn[place].ID= 0; // Oct 9, 2005
	else
	  SidmDataIn[place].ID= P[i].ID;
        SidmDataIn[place].dt= 2*(All.Time - P[i].CurrentTime);
        
        SidmDataConfirm[place].Task= -1;
        ncount++;
      }
    }

    /* 1st big communication */
    tstart= second();      
    for(level=1; level<NTask; level++) {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level; 
      
      MPI_Sendrecv(&SidmDataIn[noffset[sendTask]], 
		   nrecv[sendTask]*sizeof(struct sidmdata_in), 
		   MPI_BYTE, recvTask, TAG_ANY, 
                   &SidmDataIn[noffset[recvTask]], 
		   nrecv[recvTask]*sizeof(struct sidmdata_in), MPI_BYTE, 
		   recvTask, TAG_ANY, MPI_COMM_WORLD, &status);
    }
    tend=second();
    All.CPU_CommSum += timediff(tstart, tend);

    /* ok, all preparations are done. */
    /* SIDM: MAIN PART */

    /* Common Factors */
    /* s_a= a^{3/2} H
       ==> s_a_inverse da = a^{-1/2} dt
       ==> s_a_inverse da v_internal= a^{-1} v_phys dt = dx_comoving  
    */
    if(All.ComovingIntegrationOn) {
      s_a= All.Hubble*sqrt(All.Omega0 + All.Time*
           (1-All.Omega0-All.OmegaLambda)+pow(All.Time,3)*All.OmegaLambda);
      s_a_inverse= 1/s_a;
#if (CROSS_SECTION_TYPE == 0)
      CrossSectionCo = All.CrossSectionInternal/pow(All.Time,2);  
      C_Pmax= SAFEFACTOR * 
              BALLINVERSE*(All.DesNumNgb+All.MaxNumNgbDeviation)*2*vmax*
	      CrossSectionCo;
#elif (CROSS_SECTION_TYPE == 1)
      CrossSectionCo = All.CrossSectionInternal/pow(All.Time,2.5);
      // a^(1/2) factor considers that in s_a_inverse
      C_Pmax= SAFEFACTOR * 
              BALLINVERSE*(All.DesNumNgb+All.MaxNumNgbDeviation)*
              CrossSectionCo;
#elif (CROSS_SECTION_TYPE == 2)
      CrossSectionCo = All.CrossSectionInternal/pow(All.Time,2);  
      vc= All.YukawaVelocity/sqrt(All.Time); /* In internal unit */
      if(2.0*vmax < vc/sqrt(3.0)) {
	beta= 2.0*vmax/vc;
	v_dep= 1.0/(1.0 + beta*beta);
	C_Pmax= SAFEFACTOR * 
	        BALLINVERSE*(All.DesNumNgb+All.MaxNumNgbDeviation)*
	        2.0*vmax*v_dep*v_dep*
                CrossSectionCo;
      }
      else {
	C_Pmax= SAFEFACTOR * 
                BALLINVERSE*(All.DesNumNgb+All.MaxNumNgbDeviation)*
	        (3.0*sqrt(3.0)/16.0)*vc*
                CrossSectionCo;
      }
#elif (CROSS_SECTION_TYPE == 3)
      CrossSectionCo = All.CrossSectionInternal/pow(All.Time, 2);
      n_cross_section= All.CrossSectionPowLaw;
      v_scale= All.CrossSectionVelScale;
      C_Pmax= SAFEFACTOR * 
	      BALLINVERSE*(All.DesNumNgb+All.MaxNumNgbDeviation)*2*v_scale*
	      CrossSectionCo;
#elif (CROSS_SECTION_TYPE == 4)
      vc= All.YukawaVelocity/sqrt(All.Time);
      CrossSectionCo = All.CrossSectionInternal/pow(All.Time,2);  
      C_Pmax= SAFEFACTOR * 
              BALLINVERSE*(All.DesNumNgb+All.MaxNumNgbDeviation)*2*vmax*
	      CrossSectionCo;

#endif
    }

    else{
      s_a_inverse= 1.0;
#if (CROSS_SECTION_TYPE == 0)
      C_Pmax= SAFEFACTOR * 
	      BALLINVERSE*(All.DesNumNgb+All.MaxNumNgbDeviation)*2*vmax*
              All.CrossSectionInternal;
#elif (CROSS_SECTION_TYPE == 1)
      C_Pmax= SAFEFACTOR * 
	      BALLINVERSE*(All.DesNumNgb+All.MaxNumNgbDeviation)*
              All.CrossSectionInternal;
#elif (CROSS_SECTION_TYPE == 2)
      vc= All.YukawaVelocity;
      if(2.0*vmax < vc/sqrt(3.0)) {
	beta= 2.0*vmax/vc;
	v_dep= 1.0/(1.0 + beta*beta);
	C_Pmax= SAFEFACTOR * 
	        BALLINVERSE*(All.DesNumNgb+All.MaxNumNgbDeviation)*
	        2.0*vmax*v_dep*v_dep*
	        All.CrossSectionInternal;
      }
      else {
	C_Pmax= SAFEFACTOR * 
                BALLINVERSE*(All.DesNumNgb+All.MaxNumNgbDeviation)*
                (3.0*sqrt(3.0)/16.0)*vc*
	        All.CrossSectionInternal;
      }
#elif (CROSS_SECTION_TYPE == 3)
      n_cross_section= All.CrossSectionPowLaw;
      v_scale= All.CrossSectionVelScale;

      C_Pmax= SAFEFACTOR * 
	      BALLINVERSE*(All.DesNumNgb+All.MaxNumNgbDeviation)*2*v_scale*
              All.CrossSectionInternal;
#elif (CROSS_SECTION_TYPE == 4)
      vc= All.YukawaVelocity;
      C_Pmax= SAFEFACTOR * 
	      BALLINVERSE*(All.DesNumNgb+All.MaxNumNgbDeviation)*2*vmax*
              All.CrossSectionInternal;

#endif
      CrossSectionCo = All.CrossSectionInternal;
    }

    /* For each active particles */
    for(i=0; i<(nbuffer+ncount-nexport); i++) {
      numngb= ngb_treefind_variable(&SidmDataIn[i].Pos[0], 
                    SidmDataIn[i].Hsml, 
                    SidmDataIn[i].Type, 
                    &ngblist, &r2list);   
      SidmDataResult[i].Ngb= numngb;

      dt_h0= SidmDataIn[i].dt * s_a_inverse;

      h     = SCATKERNELFACTOR * SidmDataIn[i].Hsml; // Oct 31, 2005
      hinv  = 1.0/h;
      hinv3 = hinv*hinv*hinv;
      
      /* Initialize Variables */
      SidmDataResult[i].dv[0] = 0.0;
      SidmDataResult[i].dv[1] = 0.0;
      SidmDataResult[i].dv[2] = 0.0;
      SidmTarget[i] = 0; /* for security */

      P_max = C_Pmax*SidmDataIn[i].Mass*hinv3*dt_h0;
      // This P_max is not correct if Masses are different. 

      rand = ran2(&iseed);
      
      if(P_max < rand)
        continue;
      else if(SidmDataIn[i].ID == 0)
	continue; // Already Scattered

      nscat[0]++; /* Passed first approximation */
      
      Prob= 0.0;
      /* for each neighbors */
      for(n=0; n<numngb; n++) {
        j  = ngblist[n]+1; 
        r  = sqrt(r2list[n]); /* There was a BUG! */

	if(P[j].dVel[0] != 0.0f)
	  continue; // No double scattering. Oct 9, 2005
        
        if(r < h) {
          u = r*hinv;
          ii = (int)(u*KERNEL_TABLE);
          wk =hinv3*( Kernel[ii]  + (Kernel[ii+1]-Kernel[ii])*(u-KernelRad[ii])*KERNEL_TABLE);
        }
	
        /* Relative Velocities */
        rvx = SidmDataIn[i].Vel[0] - P[j].Vel[0];
        rvy = SidmDataIn[i].Vel[1] - P[j].Vel[1];
        rvz = SidmDataIn[i].Vel[2] - P[j].Vel[2];
        rv  = sqrt(rvx*rvx + rvy*rvy + rvz*rvz);

#if (CROSS_SECTION_TYPE == 0)
        Prob += 0.5*P[j].Mass*wk*rv*CrossSectionCo*dt_h0;
#elif (CROSS_SECTION_TYPE == 1)
        Prob += 0.5*P[j].Mass*wk*   CrossSectionCo*dt_h0;
#elif (CROSS_SECTION_TYPE == 2)
	beta= rv/vc;
	v_dep= 1.0/(1.0+beta*beta);
	Prob += 0.5*P[j].Mass*wk*rv*v_dep*v_dep*CrossSectionCo*dt_h0;
#elif (CROSS_SECTION_TYPE == 3)
        Prob += 0.5*P[j].Mass*wk*rv*pow(rv/v_scale, n_cross_section)*CrossSectionCo*dt_h0;
#elif (CROSS_SECTION_TYPE == 4)
        Prob += 0.5*P[j].Mass*wk*rv*CrossSectionCo*dt_h0;
#endif

        if(Prob < rand)
          continue;

	SidmTarget[i] = j;        
        rmass = P[j].Mass / (SidmDataIn[i].Mass + P[j].Mass);

#if (CROSS_SECTION_TYPE == 4)
	// Furture selection of angle
	beta= rv/vc;
	rand = ran2(&iseed);
	cosO = 2.0*ran2(&iseed) - 1.0;    // cos(theta) of scatter direction
	sin22 = 0.5*(1.0 - cosO);         // sin^2(theta/2)
	denom = 1.0 + beta*beta*sin22*sin22;
	if(rand >= 1/(denom*denom))
	  continue;

	if(rv == 0.0) { // debug
	  //printf("Error: rv == 0. %e %e\n", rv, 0.5*P[j].Mass*wk*rv*CrossSectionCo*dt_h0);
	  endrun(4006);
	}

	rvv[0]= rvx; rvv[1]= rvy; rvv[2]= rvz;
	random_direction(nx);


	perp(rvv, nx, nperp);

	//if(nperp[0] != nperp[0] || nperp[1] != nperp[1] || nperp[2] != nperp[2])
	//  endrun(4002);

	//assert(1.0 - cosO*cosO >= 0.0);

	sinO = 1.0 - cosO*cosO > 0.0 ? sqrt(1.0 - cosO*cosO) : 0.0;

	//if(sinO != sinO) endrun(4001);

	dvx= rmass*(-rvx + cosO*rvv[0] + sinO*rv*nperp[0]);
	dvy= rmass*(-rvy + cosO*rvv[1] + sinO*rv*nperp[1]);
	dvz= rmass*(-rvz + cosO*rvv[2] + sinO*rv*nperp[2]);

	if(dvx != dvx || dvy != dvy || dvz != dvz) {
	  endrun(4000);
	}

	//debug
	/*
	rvv[0]= cosO*rvv[0] + sinO*rv*nperp[0];
	rvv[1]= cosO*rvv[1] + sinO*rv*nperp[1];
	rvv[2]= cosO*rvv[2] + sinO*rv*nperp[2];
	rvv[0]= SidmDataIn[i].Vel[0] + dvx;
	rvv[1]= SidmDataIn[i].Vel[1] + dvy;
	rvv[2]= SidmDataIn[i].Vel[2] + dvz;
	fprintf(stderr, "%e %e\n", rv, norm(rvv));
	*/
#endif

	
        /* Scatter. Calc velocity change */

	// Better way of producing spherical random numbers
	// Sep 19, 2005
#if (CROSS_SECTION_TYPE < 4)
	// isotropic scattering
	random_direction(nx);
	dvx= rmass*(-rvx+rv*nx[0]);
	dvy= rmass*(-rvy+rv*nx[1]);
	dvz= rmass*(-rvz+rv*nx[2]);
#endif
	
        SidmDataResult[i].dv[0]= dvx;
        SidmDataResult[i].dv[1]= dvy;
        SidmDataResult[i].dv[2]= dvz;
        SidmDataConfirm[i].Task= ThisTask;
        
        break;
      }
    }  

    /* 2nd communicate contributions, and sum up */

    tstart=second();
    for(level=1; level<NTask; level++) {
      sendTask= ThisTask;
      recvTask= ThisTask ^ level; 
    
      MPI_Sendrecv(&SidmDataResult[noffset[recvTask]], 
		   nrecv[recvTask]*sizeof(struct sidmdata_out), MPI_BYTE,
		   recvTask, TAG_ANY, &SidmDataPartialResult[0],
		   nrecv[sendTask]*sizeof(struct sidmdata_out), MPI_BYTE, 
		   recvTask, TAG_ANY, MPI_COMM_WORLD, &status);

      for(i=0; i<nrecv[ThisTask]; i++) {
        SidmDataResult[noffset[ThisTask] + i].Ngb += 
	  SidmDataPartialResult[i].Ngb;
        
        // Keep First collision only
        if(SidmDataResult[noffset[ThisTask] + i].dv[0] == 0.0 &&
           SidmDataPartialResult[i].dv[0] != 0.0 ) {
          for(k=0; k<3; k++) {
            SidmDataResult[noffset[ThisTask] + i].dv[k] 
              = SidmDataPartialResult[i].dv[k];
          }
          SidmDataConfirm[noffset[ThisTask] + i].Task
            = recvTask;
        }
      }
    }
    tend=second();
    All.CPU_CommSum += timediff(tstart,tend);

    /* transfer the result to the particles */
      
    for(i=nstart, ncount=0, nexport=0, timelinecounter=startcounter; 
        timelinecounter<NumForceUpdate && ncount<nthis; 
        i=P[i].ForceFlag, timelinecounter++)
      if(P[i].Type > 0) {
        if(P[i].Type & 8) {
          place= noffset[ThisTask] + nexport; 
          nexport++;
          P[i].Type&= 7;
        }
        else
          place= nbuffer + (ncount-nexport);
        
        P[i].NgbVelDisp = SidmDataResult[place].Ngb;
        
        // Check # of neighbours
        if(P[i].NgbVelDisp < (All.DesNumNgb-All.MaxNumNgbDeviation) || 
	   (P[i].NgbVelDisp > (All.DesNumNgb+All.MaxNumNgbDeviation))) {
          // Scattering is invalid
          // Redo sidm in ensure_neighbours() function
          SidmDataConfirm[place].Task= -1;
          if(SidmDataResult[place].dv[0] != 0.0)
            nscat[2]++; // rejected;
        }
        else if(SidmDataResult[place].dv[0] != 0.0) {
          // OK. Update Velocity.
	  /* Oct 09, 2005
          P[i].Vel[0] += SidmDataResult[place].dv[0];
          P[i].Vel[1] += SidmDataResult[place].dv[1];
          P[i].Vel[2] += SidmDataResult[place].dv[2];
	  */
	  P[i].dVel[0] = SidmDataResult[place].dv[0];
          P[i].dVel[1] = SidmDataResult[place].dv[1];
          P[i].dVel[2] = SidmDataResult[place].dv[2];
          nscat[1]++; // scattered;
        }
	else {
	  // Feb 15,2006 Double scattering rejection (very rare)
          SidmDataConfirm[place].Task= -1;
	}
        ncount++;
      }
    
    /* Scattering Confirmation */
    /* Communication; 3rd time */

    tstart= second();      
    for(level=1; level<NTask; level++){
      sendTask = ThisTask;
      recvTask = ThisTask ^ level; 
      
      MPI_Sendrecv(&SidmDataConfirm[noffset[sendTask]], 
		   nrecv[sendTask]*sizeof(struct sidmdata_confirm), MPI_BYTE, 
		   recvTask, TAG_ANY, 
		   &SidmDataConfirm[noffset[recvTask]], 
		   nrecv[recvTask]*sizeof(struct sidmdata_confirm), MPI_BYTE, 
		   recvTask, TAG_ANY, MPI_COMM_WORLD, &status);
    }
    tend= second();
    All.CPU_CommSum += timediff(tstart, tend);
    
    /* Update the targets of active particles if confirmed */

    for(j=0; j<(nbuffer+ncount-nexport); j++) {
      if(SidmDataConfirm[j].Task == ThisTask) {
        if(SidmTarget[j] <= 0 || SidmTarget[j] > NumPart) {
          // for debug only
          printf("SIDM ERROR: SidmTarget is wrong\n");
          endrun(0);
        }
	// Oct 9,2005 Vel[] => dVel[] -=-> =-
        P[SidmTarget[j]].dVel[0] = -SidmDataResult[j].dv[0];
        P[SidmTarget[j]].dVel[1] = -SidmDataResult[j].dv[1];
        P[SidmTarget[j]].dVel[2] = -SidmDataResult[j].dv[2];
        
#ifdef SCATTERLOG
        /* Scatter Log */
	ScatLog.id1= SidmDataIn[j].ID;
	ScatLog.id2= P[SidmTarget[j]].ID;

	ScatLog.Hsml1= SidmDataIn[j].Hsml;
	ScatLog.Hsml2= P[SidmTarget[j]].HsmlVelDisp;

	ScatLog.x1[0]= SidmDataIn[j].Pos[0];
	ScatLog.x1[1]= SidmDataIn[j].Pos[1];
	ScatLog.x1[2]= SidmDataIn[j].Pos[2];

	ScatLog.x2[0]= P[SidmTarget[j]].PosPred[0];
	ScatLog.x2[1]= P[SidmTarget[j]].PosPred[1];
	ScatLog.x2[2]= P[SidmTarget[j]].PosPred[2];


	ScatLog.v1[0]= SidmDataIn[j].Vel[0];
	ScatLog.v1[1]= SidmDataIn[j].Vel[1];
	ScatLog.v1[2]= SidmDataIn[j].Vel[2];

	ScatLog.v2[0]= P[SidmTarget[j]].Vel[0];
	ScatLog.v2[1]= P[SidmTarget[j]].Vel[1];
	ScatLog.v2[2]= P[SidmTarget[j]].Vel[2];

	ScatLog.dv[0]= SidmDataResult[j].dv[0];
	ScatLog.dv[1]= SidmDataResult[j].dv[1];
	ScatLog.dv[2]= SidmDataResult[j].dv[2];

	fwrite(&ScatLog, sizeof(ScatLog), 1, fp);
#endif
      }
    }
    
    nstart= i; /* continue at this point in the timeline in next iteration */
    startcounter=timelinecounter;
      
    npleft-= ncount;
    ntotleft-= ntot;
  }

  // Summerize the scattering statistics
  MPI_Reduce(nscat, nscatTot, 3, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
#ifdef FINDNBRLOG
  if(ThisTask == 0){
    //printf("SIDM SCT tot= %d 1st= %d scattered= %d rejected= %d\n",   
    fprintf(stdout, "SCT %d %d %d %d\n",     
     ntot, nscatTot[0], nscatTot[1], nscatTot[2]);
  }
#endif

#ifdef SCATTERLOG
  fclose(fp);
#endif
  free(nrecv);
  free(noffset);
}


void setup_nbr_sidm(void)
{
  int     i,j,k, numngb;
  float   *r2list;
  int     *ngblist;
  double  tstart,tend;
  int     ntot,ntotleft,npleft,nthis;     
  int     nstart,nbuffer,ncount,nchunk;
  int     timelinecounter,startcounter;   
  int     *nrecv,*noffset;
  int     level,sendTask,recvTask;
  int     place, nexport;
  int     num_collisionless;
  MPI_Status status;
  
  num_collisionless= NumForceUpdate-NumSphUpdate;

  MPI_Allreduce(&num_collisionless, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  /* compute maximum size of bunch that may come from this task */
  if(num_collisionless==0)
    nthis=1;
  else
    nthis = ((double)(num_collisionless)*(All.BunchSizeSidm-NTask))/ntot + 1; 


  if(num_collisionless>0) {
    nchunk= num_collisionless/nthis;   /* number of bunches needed */  
    if((num_collisionless%nthis)>0) nchunk+=1;
  }
  else
    nchunk=0;

  nrecv=   malloc(sizeof(int)*NTask);  /* list of particle numbers that constituate current bunch */
  noffset= malloc(sizeof(int)*NTask);  /* offsets of bunches in common list */

  ntotleft= ntot;              /* particles left for all tasks together */
  npleft= num_collisionless;   /* particles left for this task */

  nstart= IndFirstUpdate;      /* first particle for this task */
  startcounter=0;

  while(ntotleft > 0) {
    if(nthis > npleft)
      nthis=npleft;

    for(i=nstart, ncount=0, nexport=0, timelinecounter=startcounter; 
        timelinecounter < NumForceUpdate && ncount<nthis; 
        i=P[i].ForceFlag, timelinecounter++) {
      if(P[i].Type > 0) {
        ncount++;
        
        for(j=0; j<3; j++) {
          if(P[i].PosPred[j] < (DomainMin[P[i].Type][j] + P[i].HsmlVelDisp))
            break;
          if(P[i].PosPred[j] > (DomainMax[P[i].Type][j] - P[i].HsmlVelDisp))
            break;
        }
        
        if(j != 3) {
	  /* particle lies NOT completely inside.
	     needs to be sent to other processors */
          nexport++;
          P[i].Type |= 8;
        }
      }
    }

    MPI_Allgather(&nexport, 1, MPI_INT, nrecv, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allreduce(&ncount, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      
    for(i=nbuffer=0; i<NTask; i++)  /* compute length of common list */
      nbuffer += nrecv[i];
      
    for(i=1, noffset[0]=0; i<NTask; i++) /* set-up offset-table */
      noffset[i]= noffset[i-1] + nrecv[i-1];

    /* fill in the own particles at the right place */
    for(i=nstart, ncount=0, nexport=0, timelinecounter=startcounter; 
        timelinecounter<NumForceUpdate && ncount<nthis; 
        i=P[i].ForceFlag, timelinecounter++) {
      if(P[i].Type>0) {
        if(P[i].Type & 8) {
          place= noffset[ThisTask] + nexport; 
          nexport++;
        }
        else
          place= nbuffer + (ncount-nexport);
      
        for(k=0; k<3; k++) {
          SidmDataIn[place].Pos[k]= P[i].PosPred[k];
        }
        SidmDataIn[place].Hsml  = P[i].HsmlVelDisp;
        SidmDataIn[place].Type  = P[i].Type&7;
        ncount++;
      }
    }

    /* 4th big communication */

    tstart= second();      
    for(level=1; level<NTask; level++) {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level; 
    
      MPI_Sendrecv(&SidmDataIn[noffset[sendTask]], 
		   nrecv[sendTask]*sizeof(struct sidmdata_in), MPI_BYTE, 
		   recvTask, TAG_ANY, &SidmDataIn[noffset[recvTask]], 
		   nrecv[recvTask]*sizeof(struct sidmdata_in), MPI_BYTE, 
		   recvTask, TAG_ANY, MPI_COMM_WORLD, &status);
    }
    tend=second();
    All.CPU_CommSum+= timediff(tstart, tend);

    /* ok, all preparations are done. 
       Now density evaluation and velocity stuff */
    
    for(i=0; i<(nbuffer+ncount-nexport); i++) {
      numngb= ngb_treefind_variable(&SidmDataIn[i].Pos[0], 
                    SidmDataIn[i].Hsml, 
                    SidmDataIn[i].Type, 
                    &ngblist, &r2list);
      
      SidmDataResult[i].Ngb= numngb;
    }  

    /* now communicate contributions, and sum up */

    tstart=second();
    for(level=1; level<NTask; level++) {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level; 
    
      MPI_Sendrecv(&SidmDataResult[noffset[recvTask]], 
		   nrecv[recvTask]*sizeof(struct sidmdata_out), MPI_BYTE, 
		   recvTask, TAG_ANY, &SidmDataPartialResult[0], 
		   nrecv[sendTask]*sizeof(struct sidmdata_out), MPI_BYTE,
		   recvTask, TAG_ANY, MPI_COMM_WORLD, &status);

      for(i=0; i<nrecv[ThisTask]; i++) {
        SidmDataResult[noffset[ThisTask] + i].Ngb += 
	  SidmDataPartialResult[i].Ngb;          
      }
    }
    tend=second();
    All.CPU_CommSum+= timediff(tstart,tend);

    /* transfer the result to the particles */
    
    for(i=nstart, ncount=0, nexport=0, timelinecounter=startcounter; 
        timelinecounter<NumForceUpdate && ncount<nthis; 
        i=P[i].ForceFlag, timelinecounter++) {
      if(P[i].Type>0) {
        if(P[i].Type & 8) {
          place= noffset[ThisTask] + nexport; 
          nexport++;
          P[i].Type&= 7;
        }
        else
          place= nbuffer + (ncount-nexport);
        
        P[i].NgbVelDisp = SidmDataResult[place].Ngb;
        ncount++;
      }
    }
      
    nstart=i; /* continue at this point in the timeline in next iteration */
    startcounter=timelinecounter;
      
    npleft-= ncount;
    ntotleft-= ntot;
  }

  free(nrecv);
  free(noffset);
}




/* the function below detects particles that have a number of neighbours 
 * outside the allowed tolerance range. For these, particles the smoothing
 * length is adjusted accordingly. Note that the smoothing length is
 */
void sidm_ensure_neighbours(int mode)
{
#define MAXITER 30

  int    i, ntot, last=0;
  float  *r2list;
  int    *ngblist, count, candidates;
  int    iter=0;
  double save;
  
  /*
  int    IndFirstUpdateBackup, NumForceUpdateBackup;
  
  IndFirstUpdateBackup= IndFirstUpdate;
  NumForceUpdateBackup= NumForceUpdate;
  for(i=IndFirstUpdate, count=0; count<NumForceUpdate; 
      i=P[i].ForceFlag, count++){
    P[i].ForceFlagBackup= P[i].ForceFlag;
  }
  */

  for(i=IndFirstUpdate, count=0, candidates=0; count<NumForceUpdate; 
      i=P[i].ForceFlag, count++) {
    if(P[i].Type>0) {
      if(P[i].NgbVelDisp < (All.DesNumNgb-All.MaxNumNgbDeviation) || 
	 (P[i].NgbVelDisp > (All.DesNumNgb+All.MaxNumNgbDeviation)))
        candidates++;
    }
  }

  MPI_Reduce(&candidates, &ntot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ntot, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  if(ntot > 0) {
    //#ifdef FINDNBRLOG
      //      if(ThisTask==0)
      //{
      //printf("\n%d particles have too few/too many neighbours in sidm calculation!\n", ntot);
      // printf("Now fixing that...\n"); 
      // }
    //#endif
      
    for(i=IndFirstUpdate, count=0; count<NumForceUpdate; 
	i=P[i].ForceFlag, count++)   
      P[i].Left= P[i].Right= 0;
      
    do {
      for(i=1+N_gas, NumForceUpdate= 0, NumSphUpdate=0; i<=NumPart; i++) { 
	if( P[i].NgbVelDisp < (All.DesNumNgb-All.MaxNumNgbDeviation) || 
	    P[i].NgbVelDisp > (All.DesNumNgb+All.MaxNumNgbDeviation)) {
          if(P[i].Left>0 && P[i].Right>0)
            if((P[i].Right-P[i].Left) < 1.0e-3 * P[i].Left)
              continue;
	  
          if(NumForceUpdate==0)
            IndFirstUpdate= i;
          else
            P[last].ForceFlag= i;

          NumForceUpdate++;
          last=i;
          
          if(P[i].NgbVelDisp < (All.DesNumNgb-All.MaxNumNgbDeviation))
            P[i].Left= dmax(P[i].HsmlVelDisp, P[i].Left);
          else
            if(P[i].Right!=0)
              {
            if(P[i].HsmlVelDisp<P[i].Right)
              P[i].Right= P[i].HsmlVelDisp;
              }
            else
              P[i].Right= P[i].HsmlVelDisp;
        }
      }
      
      MPI_Allreduce(&NumForceUpdate, &ntot, 1, MPI_INT, MPI_SUM, 
		    MPI_COMM_WORLD);
      
      if(ntot > 0) {
	//#ifdef FINDNBRLOG
	//  if(ThisTask==0)
	//printf("ngb iteration %d.  still %d particles\n", iter, ntot);
	//#endif
	for(i=IndFirstUpdate, count=0; count<NumForceUpdate; 
	    i=P[i].ForceFlag, count++) {
          if(iter >= 20) {
	    printf("i=%d ID=%d Hsml=%g Left=%g Right=%g Ngbs=%d Right-Left=%g\n   pos=(%g|%g|%g)\n",
		   i, P[i].ID, P[i].HsmlVelDisp, P[i].Left, P[i].Right, 
		   P[i].NgbVelDisp, P[i].Right-P[i].Left,
		   P[i].PosPred[0], P[i].PosPred[1], P[i].PosPred[2]);
	  }
          
          if(iter == MAXITER) {
	    printf("ThisTask=%d Mi=(%g|%g|%g) Ma=(%g|%g|%g)\n", ThisTask,
		   DomainMin[P[i].Type][0], DomainMin[P[i].Type][1], 
		   DomainMin[P[i].Type][2], DomainMax[P[i].Type][0], 
		   DomainMax[P[i].Type][1], DomainMax[P[i].Type][2]);
	    printf("i=%d ID=%d coord=(%g|%g|%g)\n", i, P[i].ID, 
		   P[i].PosPred[0], P[i].PosPred[1], P[i].PosPred[2]);
	    printf("ngb_treefind= %g\n", sqrt(ngb_treefind(P[i].PosPred , 
			  All.DesNumNgb, 0, P[i].Type, &ngblist, &r2list)));
	  }
          
          if(P[i].Left==0 || P[i].Right==0) {
	    if(P[i].Right==0 && P[i].NgbVelDisp<15 && 
	       NtypeLocal[P[i].Type]>All.DesNumNgb) {
              P[i].HsmlVelDisp= sqrt(ngb_treefind(P[i].PosPred, All.DesNumNgb, 
					    0, P[i].Type, &ngblist, &r2list));  
            }
	    else {
              P[i].HsmlVelDisp=  P[i].HsmlVelDisp*( 0.5 + 0.5*pow(P[i].NgbVelDisp/((double)All.DesNumNgb), -1.0/3));
            }
	  }
          else {
	    P[i].HsmlVelDisp=  0.5*(P[i].Left + P[i].Right);
	  }
        }
          
	sidm();

	iter++;

	if(iter > MAXITER) {
          fprintf(stdout, "failed to converge in function ensure_neighbours()\n");
          endrun(1155);
        }
      }
    } while(ntot > 0);
      
    if(mode == 0) {   /* restore timeline to active particles */
      /*
      IndFirstUpdate= IndFirstUpdateBackup;
      NumForceUpdate= NumForceUpdateBackup;
      for(i=IndFirstUpdateBackup, count=0; count<NumForceUpdateBackup; 
	  i=P[i].ForceFlagBackup, count++){
	P[i].ForceFlag= P[i].ForceFlagBackup;
      }
      */
      
      save= All.TimeStep;
      find_next_time();
      All.TimeStep= save;
    }
    else { /* make all particles active again */
      for(i=1; i<=NumPart; i++) 
        P[i].ForceFlag=i+1;
      
      P[NumPart].ForceFlag=1; 
      IndFirstUpdate=1; 
      NumForceUpdate=NumPart; 
      NumSphUpdate=N_gas;
    }
  }
#undef MAXITER
}

double getvmax()
{
  int j;
  double v2, vmax=0.0;

  for(j=1 ; j<= NumPart ; j++) {
    v2=   P[j].Vel[0]*P[j].Vel[0] 
        + P[j].Vel[1]*P[j].Vel[1]
        + P[j].Vel[2]*P[j].Vel[2];
    if(vmax < v2)
      vmax= v2;
  }
  vmax= sqrt(vmax);

#ifdef FINDNBRLOG
  if(ThisTask == 0)
    fprintf(stdout, "Vmax= %g Processor %d\n", vmax, ThisTask);
#endif

  return vmax;
}

void update_node_sidm(void)
{
  int k;
  for(k=0; k<n_scat_particles; k++)
    update_node_of_scat_particle(scat_particles[k]-1); 
}

#endif
