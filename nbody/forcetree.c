#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"

/*
 *  This file contains the computation of the gravitational force 
 *  by means of a tree. The type of tree implemented
 *  is a geometrical oct-tree, starting from a cube encompassing
 *  all particles. This cube is automatically found.
 *  For particle types with different softenings, separate trees
 *  are always constructed. If 'COMBINETREES' is defined during 
 *  compilation, different particle types with the same softening
 *  are put into the same tree. 
 *  The multipole moments of branches of the tree are dynamically
 *  updated during the tree walk, if needed.
 */


#define KERN_LEN   10000  /* length of arrays for tabulated functions  */ 

static struct NODE 
{ 
  float     len;                              /* sidelength of treenode */
  float     mass;                             /* total mass represented by node */       
  float     s[3], center[3];                  /* center-of-mass/geometrical center */
  float     Q11, Q22, Q33, Q12, Q13, Q23, P;  /* quadrupol tensor */  
  float     vs[3];                            /* velocity of center-of-mass */
  float     tilu;                             /* time of last update of node */
  float     oc;                               /* variable for opening criterion */
  float     hmax;                             /* maximum SPH smoothing length in node. Only used for gas */
#ifdef BMAX
  float     bmax2;                            /* holds the max. sqr. distance of the com to any point in the cell */
#endif
  int4byte  suns[8];                          /* ancestor nodes */
  int4byte  sibling;                          /* gives next node if this node needs not to be opened in tree walk */
  int4byte  partind;                          /* index of first particle in node */
  int4byte  count;                            /* total number of particles in node and below  */
  int4byte  cost;                             /* counts number of node interactions */
} *nodes, *nodes_base;

static int4byte  *next;       /* provides link-list for particle groups (to find all particles in a node) */
static int4byte  *nextnode;   /* gives next node in tree walk */
static int4byte  *father;     /* gives parent node in tree */

static struct particle_data *Part;  /* gives the particle data (starting at offset 0) */
static struct sph_particle_data *SphPart;  /* gives the sph particle data (starting at offset 0) */

static double S_of_a_inverse;
static double BoxHalf, Box;

static int trees[6];           /* holds first node of each tree */
static int numnodestree[6];    /* number of (internal) nodes in each tree */

static int numnodestotal;      /* total number of internal nodes */
static int MaxNodes;           /* maximum allowed number of internal nodes */

static int num_in_tree;

static int treecost[6];          /* sums the particle-particle interactions  (diagnostic) */
static int treecost_quadru[6];   /* sums the node-particle interactions      (diagnostic) */

static double  knlrad  [KERN_LEN+1],   /* tabulated functions for force/potential softening */
               knlforce[KERN_LEN+1],
               knlpot  [KERN_LEN+1],
               knlW2   [KERN_LEN+1],
               knlW3   [KERN_LEN+1],
               knlW4   [KERN_LEN+1];

static int last; /* auxialiary variable used to set-up non-recursive walk */



void dump_particles(void);



/* Constructs the gravitational oct-tree(s).
 * For each particle type, we normally construct a separate tree,
 * except when COMBINETREES is defined. Then types with identical
 * softening are combined into a single tree, but the gas particles
 * are nevertheless packed in their own tree (because of neighbour
 * search).
 */
int force_treebuild(void)
{

  int i, tr;
  int typelist[6];
  int nfree;
  int numnodes;
  int ntype[6];   /* holds the number of particles of each type */

  Box= All.BoxSize;
  BoxHalf= All.BoxSize/2;
  Part= &P[1];
  SphPart= &SphP[1];  /* map offset 0 to offset 1 */
  nodes= nodes_base - All.MaxPart;  
  num_in_tree= NumPart;

  for(tr=0; tr<6; tr++)
    ntype[tr]=0;

  for(i=0; i<NumPart; i++)
    ntype[ Part[i].Type ]++;

  for(tr=0; tr<6; tr++)
    NtypeLocal[tr]= ntype[tr];

  MPI_Allreduce(ntype, Ntype, 5, MPI_INT, MPI_SUM, MPI_COMM_WORLD);


  /* 'nfree' gives the next free internal tree node. 
   *  Note that the inidices are shifted by NumPart, i.e. 
   *  nodes[NumPart] is the first allocated internal node 
   */
  for(tr=0, nfree=All.MaxPart, numnodestotal=0; tr<6; tr++)
    {
      treecost[tr]= treecost_quadru[tr]=0;

      if(ntype[tr]>0)
	{
	  for(i=0; i<6; i++)
	    typelist[i]=0;

	  typelist[tr]= 1;
	  ntype[tr]= 0;

#ifdef COMBINETREES
	  if(tr>0) /* always make a separate tree for the gas particles */
	    for(i=tr+1; i<6; i++)
	      if(All.SofteningTable[tr] == All.SofteningTable[i] && ntype[i]>0)
		{
		  typelist[i]=1;
		  ntype[i]=0;
		}
#endif

	  trees[tr]= nfree;	  
	  
	  nfree = force_treebuild_single(trees[tr], typelist, &numnodes);

	  numnodestree[tr]= numnodes;
	}
      else
	{
	  trees[tr]= -1;
	  numnodestree[tr]= 0;
	}
    }

  return numnodestotal;
}




/* packs the particles of those types that are flagged in 'typelist' into an oct-tree
 * and computes the multipole moments for the internal nodes. 
 */
int force_treebuild_single(int startnode, int *typelist, int *creatednodes)
{
  int  i, j, p, subnode, fak, fp, parent,type;
  int  nfree, th, nn, ff, numnodes;
  struct NODE *nfreep;
  double xmin[3], xmax[3], len, epsilon;

  for(fp=0; fp<NumPart; fp++)          /* find first particle of this type */
    if(typelist[Part[fp].Type])
      break;
  
  epsilon= All.SofteningTable[Part[fp].Type];
  
  for(j=0; j<3; j++)                        /* find enclosing rectangle */
    xmin[j]= xmax[j]= Part[fp].PosPred[j];

  for(i=fp+1; i<NumPart; i++)
    if(typelist[Part[i].Type])
      for(j=0; j<3; j++)
	{
	  if(Part[i].PosPred[j]>xmax[j]) 
	    xmax[j]= Part[i].PosPred[j];
	  if(Part[i].PosPred[j]<xmin[j]) 
	    xmin[j]= Part[i].PosPred[j];
	}

  for(type=0; type<5; type++)
    if(typelist[type]) /* keep the extent */
      for(j=0; j<3; j++)  
	{
	  DomainMin[type][j]= xmin[j];
	  DomainMax[type][j]= xmax[j];
	}

  for(j=1, len=xmax[0]-xmin[0]; j<3; j++)  /* determine maxmimum extension */
    if((xmax[j]-xmin[j])>len)
      len=xmax[j]-xmin[j];

  len*=1.01;

  /* create a root node and insert first particle of correct type as its leaf */

  nfree= startnode;  nfreep= &nodes[nfree];

  for(j=0; j<3; j++)
    nfreep->center[j]= (xmax[j]+xmin[j])/2;
  nfreep->len= len;

  father[nfree]= -1; 

  for(i=0; i<8; i++)
    nfreep->suns[i]= -1;
  
  for(j=0,subnode=0,fak=1; j<3; j++,fak<<=1)
    if(Part[fp].PosPred[j] > nfreep->center[j])
      subnode+=fak;
  nfreep->suns[subnode]= fp;

  father[fp]= nfree;

  nfreep->count= 1;
  nfreep->partind= fp;
  nfreep->sibling= -1;

  next[fp]= -1;

  numnodes=1; numnodestotal++; nfree++; nfreep=&nodes[nfree];
  if(numnodestotal>=MaxNodes)
    {
      printf("task %d: maximum number %d of tree-nodes reached.\n", ThisTask, numnodestotal);
      printf("for particle %d\n", fp);
      dump_particles();
      endrun(1);
    }

  for(i=fp+1; i<NumPart; i++)    /* insert all other particles */
    {
      if(typelist[Part[i].Type] == 0)
	continue;

      th= startnode;
      
      while(1)
	{
	  if(th >= All.MaxPart)  /* we are dealing with an internal node */
	    {
	      nodes[th].count++;

	      for(j=0, subnode=0, fak=1; j<3; j++, fak<<=1)
		if(Part[i].PosPred[j] > nodes[th].center[j])
		  subnode+=fak;

	      nn= nodes[th].suns[subnode];

	      if(nn>=0)
		th=nn;
	      else
		{
		  /* here we have found an empty slot where we can 
		   * attach the new particle as sleaf 
		   */
		  
		  nodes[th].suns[subnode]= i;
		  father[i]= th;

		  /* find last particle of link-list in the given node. 
		   * Note: nodes[th].count is already increased 
		   */
                  for(j=0, p= nodes[th].partind; j<(nodes[th].count-2); j++)
		    p= next[p];

		  /* insert particle into link-list */
		  next[i]= next[p];
		  next[p]= i;
		    
		  break;
		}
	    }
	  else
	    {
	      /* node is occupied with only one particle (a leaf)   
	       * need to generate a new internal node at this point 
	       */
	      
	      p= th;
	      parent= father[p]; 

	      for(subnode=0; subnode<8; subnode++)
		if(nodes[parent].suns[subnode]==p)
		  break;

	      nodes[parent].suns[subnode]= nfree; 

	      nfreep->len= nodes[parent].len/2;
	      
	      for(j=0, fak=1; j<3; j++, fak<<=1)
		if(subnode&fak)
		  nfreep-> center[j]= nodes[parent].center[j] + nfreep->len/2;
		else
		  nfreep-> center[j]= nodes[parent].center[j] - nfreep->len/2;		  

	      
	      for(j=0;j<8;j++)
		nfreep->suns[j]= -1;

	      father[nfree]= parent;
	      nfreep->partind= p;
	      nfreep->count= 1;
	      nfreep->sibling= -1;

	      for(j=0, subnode=0, fak=1; j<3; j++, fak<<=1)
		if(Part[p].PosPred[j] > nodes[nfree].center[j])
		  subnode+=fak;

	      if(nfreep->len < 1.0e-3*epsilon)
		{
		  /* seems like we're dealing with particles   
		   * at identical locations. randomize subnode (well below the gravity scale). 
		   */
		  subnode= (int) (8.0*rand()/(RAND_MAX+1.0));	
		}

	      nfreep->suns[subnode]= p;
	      
	      father[p]= nfree;
			
	      th= nfree; /* resume trying to insert the new particle at the newly created internal node */
	  
	      numnodes++; numnodestotal++; nfree++; nfreep=&nodes[nfree];

	      if(numnodestotal>=MaxNodes)
		{
		  printf("task %d: maximum number %d of tree-nodes reached.\n",ThisTask, numnodestotal);
		  printf("for particle %d\n", i);
		  dump_particles();
		  endrun(1);
		}
	    }
	}
    }



  /* now finish up tree 
   */ 

#ifdef RECTREECONS
  force_update_node_recursive(startnode);
#endif

  for(i=0, th=startnode; i<numnodes; i++, th++)
    {
      nodes[th].tilu= All.Time;
      nodes[th].cost= 0; 
#ifndef RECTREECONS
      force_update_node(th, 0);  /* center-of-mass and quadrupol computation */
#endif

      for(j=7,nn=-1; j>=0; j--)  /* preparations for non-recursive walk */
	{
	  if(nodes[th].suns[j] >= 0)
	    {
	      if(nodes[th].suns[j] >= All.MaxPart) /* write it only into suns that are internal nodes */
		nodes[nodes[th].suns[j]].sibling= nn;
	      nn= nodes[th].suns[j];
	    }
	}
    }

  last= -1;
  force_setupnonrecursive(startnode);  	/* set up non-recursive walk */
  nextnode[last]= -1;

  for(i=0,th=startnode; i<numnodes; i++, th++)
    if(nodes[th].sibling < 0)
      {
	ff=th;
	nn= nodes[th].sibling;

	while(nn<0)
	  {
	    ff=father[ff];
	    if(ff<0)
	      break;
	    nn=nodes[ff].sibling;
	  }
	
	nodes[th].sibling= nn;
      }

  *creatednodes= numnodes;

  return nfree;
} 


/* This routine walks the full tree once and 
 * generates the indices nextnode[] for
 * each node. Together with the sibling pointers
 * this allows non-recursive tree-walk without
 * having to loop over the suns[]-array.
 */
void force_setupnonrecursive(int no)
{
  int i, nn;
  
  if(last>=0)
    nextnode[last]= no;

  last= no;

  if(no>=All.MaxPart)
    for(i=0; i<8; i++)
      if((nn=nodes[no].suns[i])>=0)
	force_setupnonrecursive(nn);

}



/* this routine computes the multipole moments for all particles
 * of the given internal node. Note that this computation is done
 * in double precision with positions taken relative to the node 
 * geometrical center/center-of-mass. Note that the geometrical
 * center will be overwritten by the center-of-mass in the first call.
 * Particle positions are predicted to the current time if 'flag' is set.
 */
void force_update_node(int no, int flag)
{
  int    i, k, p;
  struct particle_data *pa;
  struct NODE *nop;
  double dt_h0=0, dt, extmax;
  double rel[3], s[3], vs[3], mass, Q11, Q22, Q33, Q12, Q23, Q13, P;
  float  hmax;
#ifdef BMAX
  double r2, dx;
#endif

  nop= &nodes[no];
  
#ifdef DIAG
  if(flag)
    {
      Num_nodeupdates++;
      Num_nodeupdate_particles+= nop->count;
    }
#endif

  for(i=0; i<3; i++)
    s[i]= vs[i]= 0;
  
  mass= Q11= Q22= Q33= Q12= Q13= Q23= P = 0;

  hmax= extmax=0;

  for(k=0, p= nop->partind; k<nop->count; k++, p= next[p])
    {
      pa = &Part[p];
      
      if(flag)
	{
	  dt= All.Time - pa->CurrentTime; 
	  dt_h0= dt*S_of_a_inverse;  
	}

      for(i=0; i<3; i++)
	{
	  if(flag)
	    {
	      /* need to predict the particle (not guaranteed that this is already done) */
	      /* however, one cannot predict velocity by using Acc, since the latter may have 
		 already been calculated in an earlier particle group, then it would have 
		 the wrong value. Use most recently predicted velocity instead. */

	      pa->PosPred[i]= pa->Pos[i] + pa->Vel[i]*dt_h0;
	    }

	  rel[i]= pa->PosPred[i] - nop->center[i];

	  s[i]  += pa->Mass*rel[i];
	  vs[i] += pa->Mass*pa->VelPred[i];
	  if(extmax < fabs(rel[i]))
	    extmax = fabs(rel[i]);
	}

      mass += pa->Mass;

      Q11 += pa->Mass*rel[0]*rel[0];
      Q22 += pa->Mass*rel[1]*rel[1];
      Q33 += pa->Mass*rel[2]*rel[2];
      Q12 += pa->Mass*rel[0]*rel[1];
      Q13 += pa->Mass*rel[0]*rel[2];
      Q23 += pa->Mass*rel[1]*rel[2];

      P +=   pa->Mass*(rel[0]*rel[0] + rel[1]*rel[1] + rel[2]*rel[2]);


      if(pa->Type==0)
	if(SphPart[p].Hsml > hmax)
	  hmax= SphPart[p].Hsml;
    }
 
  if(flag)
    {
      nop->tilu= All.Time;

      extmax *= 2;
      if(extmax > nop->len)
	{
	  nop->len= extmax;
	  /* force_update_size_of_parent_node(no); */
	}
    }
#ifdef SIDM
      if(extmax > nop->len)
	{
	  nop->len= extmax;
	  /* force_update_size_of_parent_node(no); */
	}
#endif

  if(mass)
    for(i=0; i<3; i++)
      {
	s[i]  /= mass;
	vs[i] /= mass;
      }

  Q11 -= mass*s[0]*s[0];
  Q22 -= mass*s[1]*s[1];
  Q33 -= mass*s[2]*s[2];
  Q12 -= mass*s[0]*s[1];
  Q23 -= mass*s[1]*s[2];
  Q13 -= mass*s[0]*s[2];
  P   -= mass*(s[0]*s[0] + s[1]*s[1] + s[2]*s[2]); 
	  
  for(i=0; i<3; i++)
    {
      s[i]+= nop->center[i];
      nop->s[i]=   s[i];
      nop->vs[i]= vs[i];
    }
  
  nop->mass= mass;
  nop->Q11= Q11;
  nop->Q22= Q22;
  nop->Q33= Q33;
  nop->Q12= Q12;
  nop->Q13= Q13;
  nop->Q23= Q23;
  nop->P  = P;

  nop->oc = nop->len * nop->len;
  nop->oc = nop->mass * nop->oc * nop->oc;     /* used in cell-opening criterion */

  nop->hmax = hmax;
#ifdef BMAX
  for(i=0, r2=0; i<3; i++)
    {
      dx= fabs(s[i]- nop->center[i]) + 0.5*nop->len;
      r2+= dx*dx;
    }
  nop->bmax2= r2;
#endif
}




/* this routine computes the multipole moments for the given node
 * and all its subnodes using a recursive computation.
 * Note that the moments of the daughter nodes are already stored 
 * in single precision. For very large particle numbers, loss
 * of precision may results for certain particle distributions
 */
void force_update_node_recursive(int no)
{
  int    i, j, p;
  struct particle_data *pa;
  struct NODE *nop, *subnop;
  double rel[3], s[3], vs[3], mass, Q11, Q22, Q33, Q12, Q23, Q13, P;
  float  hmax;
#ifdef BMAX
  double r2, dx;
#endif

  nop= &nodes[no];
 
  for(i=0; i<3; i++)
    s[i]= vs[i]= 0;
  
  mass= Q11= Q22= Q33= Q12= Q13= Q23= P = 0;

  hmax= 0;

  for(j=0; j<8; j++)
    {
      if(nodes[no].suns[j] >= 0)
	{
	  if(nodes[no].suns[j] >= NumPart) /* a node */
	    {
	      force_update_node_recursive(nodes[no].suns[j]);
	      
	      subnop= &nodes[nodes[no].suns[j]];

	      mass+= subnop->mass;
	      for(i=0; i<3; i++)
		{
		  s[i]  += subnop->mass*subnop->s[i];
		  vs[i] += subnop->mass*subnop->vs[i];
		}
	      
	      if(subnop->hmax > hmax)
		hmax= subnop->hmax;
	    }
	  else
	    {
	      p  = nodes[no].suns[j];
	      pa = &Part[p];

	      mass+= pa->Mass;
	      for(i=0; i<3; i++)
		{
		  s[i]  += pa->Mass*pa->PosPred[i];
		  vs[i] += pa->Mass*pa->VelPred[i];
		}

	      if(pa->Type==0)
		if(SphPart[p].Hsml > hmax)
		  hmax= SphPart[p].Hsml;	
	    }
	}
    }
  
  if(mass)
    for(i=0; i<3; i++)
      {
	s[i] /= mass;
	vs[i]/= mass;
      }


  for(j=0; j<8; j++)
    {
      if(nodes[no].suns[j] >= 0)
	{
	  if(nodes[no].suns[j] >= NumPart) /* a subnode */
	    {
	      subnop= &nodes[nodes[no].suns[j]];
	      for(i=0; i<3; i++)
		rel[i]= subnop->s[i] - s[i];
	      
	      Q11 += subnop->mass*rel[0]*rel[0];
	      Q22 += subnop->mass*rel[1]*rel[1];
	      Q33 += subnop->mass*rel[2]*rel[2];
	      Q12 += subnop->mass*rel[0]*rel[1];
	      Q13 += subnop->mass*rel[0]*rel[2];
	      Q23 += subnop->mass*rel[1]*rel[2];
	      P   +=subnop->mass*(rel[0]*rel[0] + rel[1]*rel[1] + rel[2]*rel[2]);

	      Q11 += subnop->Q11;
	      Q22 += subnop->Q22;
	      Q33 += subnop->Q33;
	      Q12 += subnop->Q12;
	      Q13 += subnop->Q13;
	      Q23 += subnop->Q23;
	      P   += subnop->P;
	    }
	  else
	    {
	      pa = &Part[nodes[no].suns[j]];
	      for(i=0; i<3; i++)
		rel[i]= pa->PosPred[i] - s[i];

	      Q11 += pa->Mass*rel[0]*rel[0];
	      Q22 += pa->Mass*rel[1]*rel[1];
	      Q33 += pa->Mass*rel[2]*rel[2];
	      Q12 += pa->Mass*rel[0]*rel[1];
	      Q13 += pa->Mass*rel[0]*rel[2];
	      Q23 += pa->Mass*rel[1]*rel[2];
	      P   += pa->Mass*(rel[0]*rel[0] + rel[1]*rel[1] + rel[2]*rel[2]);
	    }
	}
    }

  for(i=0; i<3; i++)
    {
      nop->s[i]=   s[i];
      nop->vs[i]= vs[i];
    }
 
  nop->mass= mass;
  nop->Q11= Q11;
  nop->Q22= Q22;
  nop->Q33= Q33;
  nop->Q12= Q12;
  nop->Q13= Q13;
  nop->Q23= Q23;
  nop->P  = P;

  nop->oc = nop->len  * nop->len;
  nop->oc = nop->mass * nop->oc * nop->oc;     /* used in cell-opening criterion */

  nop->hmax = hmax;
#ifdef BMAX
  for(i=0, r2=0; i<3; i++)
    {
      dx= fabs(s[i]- nop->center[i]) + 0.5*nop->len;
      r2+= dx*dx;
    }
  nop->bmax2= r2;
#endif
}


/* this routine propagates a change of extension of a tree-node
 * upwards in the tree. 
 */
void force_update_size_of_parent_node(int no)
{
  int     fa;
  float   maxext;

  if(no>0)
    {
      fa= father[no];

      if(fa>0)
	{
	  maxext= 2*fabs(nodes[no].center[0]-nodes[fa].center[0]) + nodes[no].len;
	  
	  if(maxext > nodes[fa].len || nodes[no].hmax > nodes[fa].hmax)
	    {
	      if(maxext > nodes[fa].len) 
		nodes[fa].len=  1.0001*maxext;

	      if(nodes[no].hmax > nodes[fa].hmax)
		nodes[fa].hmax= 1.0001*nodes[no].hmax;
	      
	      force_update_size_of_parent_node(fa);
	    }
	}
    }
}



 
/* here we assign the accumulated cost of internal nodes to
 * the particles. This is done by walking backwords from each
 * leaf (particle) to the root, and by summing up contributions
 * from the internal nodes. In this way, the cost of each internal
 * node is evenly ditributed onto the particles that are represented
 * by it.
 */
void force_costevaluate(void) 
{
  int i, fa;
  
  for(i=0; i<num_in_tree; i++)  
    {  
      fa= father[i];  
	  
      while(fa>=0)  
	{  
	  Part[i].GravCost += ((float) nodes[fa].cost)/nodes[fa].count;   
	  fa= father[fa];  
	}  
    }  
}



/* This routine computes the gravitational force for a given particle
 * in the communication buffer.
 * Each separate tree is walked, and depending on the value of 
 * TypeOfOpeningCriterion either the geometrical BH cell-opening criterion, 
 * or our new `relative' criterion is selected.
 */
void force_treeevaluate(int target, double one_over_s_of_a) 
{
  int    i, tr;
  double epsilon;

  S_of_a_inverse= one_over_s_of_a;

  GravDataPotential[target]= 0;

  for(i=0; i<3; i++)
    GravDataResult[target].Acc[i]= 0;
  
  for(tr=0; tr<6; tr++)
    {
      if(numnodestree[tr])
	{
	  epsilon= dmax(All.SofteningTable[tr], All.SofteningTable[GravDataIn[target].Type]);
	  if(All.TypeOfOpeningCriterion==0 || GravDataIn[target].OldAcc==0)
	    force_treeevaluate_single_BH(tr, target, epsilon);
	  else
	    force_treeevaluate_single(tr, target, epsilon);
	}
    }
}




/*  evaluate the force contribution of a tree based on the
 *  ordinary BH-criterion 
 */
void force_treeevaluate_single_BH(int tree, int targetpart, double epsilon)  
{
  int   no;
  struct NODE *nop;
  int p,ii;
  double r2,dx,dy,dz,r,fac,u,h,h2_inv,h3_inv,h5_inv,h4_inv,h6_inv,ff;
  double wf,wp,w2,w3,w4,potq;
  float  *tppos;   /* target particle for treewalk */
  double *tppot, *tpaccel;
  double q11dx,q12dy,q13dz,q12dx,q22dy,q23dz,q13dx,q23dy,q33dz;
  double r_inv,r2_inv,r3_inv,r5_inv;
  double h_inv;
  double dt,dsx,dsy,dsz;
#ifdef PERIODIC
  double fcorr[3];
#endif
  
  h = 2.8*epsilon;
  h_inv=1/h;
  h2_inv=h_inv*h_inv;
  h3_inv=h2_inv*h_inv;
  h4_inv=h2_inv*h2_inv;
  h5_inv=h2_inv*h3_inv;
  h6_inv=h3_inv*h3_inv;


  tppos   = &GravDataIn[targetpart].Pos[0];
  tpaccel = &GravDataResult[targetpart].Acc[0];
  tppot   = &GravDataPotential[targetpart];

  no= trees[tree];  nop= &nodes[no];

  while(no>=0)
    {
      if(no < All.MaxPart)   /* single particle */
	{
	  p=no;  /* the index of the node is the index of the particle */

	  /* we do the prediction right here */

	  dt=(All.Time - Part[p].CurrentTime) * S_of_a_inverse;  /*  we need here 1.0/S(a) */	  
	  /* observe the sign ! */
	  /* this vector is -y in my thesis notation */


	  /* the following explicit type-casts are done because the computation in 
             brackets is carried out in double precision, while tppos is float...
             we want to guarantee dx=0 if p is the particle at tppos */

	  dx=(float)(Part[p].Pos[0]+Part[p].Vel[0]*dt) - tppos[0]; 
	  dy=(float)(Part[p].Pos[1]+Part[p].Vel[1]*dt) - tppos[1];     
	  dz=(float)(Part[p].Pos[2]+Part[p].Vel[2]*dt) - tppos[2];

#ifdef PERIODIC
	  while(dx>BoxHalf) dx-=Box;
	  while(dy>BoxHalf) dy-=Box;
	  while(dz>BoxHalf) dz-=Box;
	  while(dx<-BoxHalf) dx+=Box;
	  while(dy<-BoxHalf) dy+=Box;
	  while(dz<-BoxHalf) dz+=Box;
#endif

	  r2=dx*dx+dy*dy+dz*dz;

	  r=sqrt(r2);  
	   
	  u=r*h_inv;

	  Part[p].GravCost +=1; 
#ifdef DIAG
	  treecost[tree] +=1;	  	  
#endif
	  if(u>=1)
	    {
	      r_inv=1/r;

	      fac= Part[p].Mass*r_inv*r_inv*r_inv;
	  
	      tpaccel[0]+=dx*fac;
	      tpaccel[1]+=dy*fac;
	      tpaccel[2]+=dz*fac;
#ifdef POTENTIAL_SIMULTO_FORCE
	      *tppot-= Part[p].Mass*r_inv;
#endif
	    }
	  else
	    {
	      ii = (int)(u*KERN_LEN); ff=(u-knlrad[ii])*KERN_LEN;
	      wf=knlforce[ii]+(knlforce[ii+1]-knlforce[ii])*ff;
	      wp=knlpot[ii]+(knlpot[ii+1]-knlpot[ii])*ff;

	      if(u>1.0e-4) 
		{
		  fac=Part[p].Mass*h_inv*h_inv*h_inv*wf;
		  
		  tpaccel[0]+=dx*fac;
		  tpaccel[1]+=dy*fac;
		  tpaccel[2]+=dz*fac;
		}
#ifdef POTENTIAL_SIMULTO_FORCE
	      *tppot+= Part[p].Mass*h_inv*wp;
#endif
	    }

#ifdef PERIODIC
	  if(u>1.0e-4)
	    {
	      ewald_corr(dx,dy,dz,fcorr);
	      
	      tpaccel[0]+=  Part[p].Mass*fcorr[0];
	      tpaccel[1]+=  Part[p].Mass*fcorr[1];
	      tpaccel[2]+=  Part[p].Mass*fcorr[2];
	    }
#endif
	  no= nextnode[no];  nop= &nodes[no];
	}
      else
	{
	  dt=(All.Time - nop->tilu) * S_of_a_inverse;  /*  we need here 1.0/S(a) */	  

	  dsx=nop->vs[0]*dt;
	  dsy=nop->vs[1]*dt;
	  dsz=nop->vs[2]*dt;

	  if((dsx*dsx + dsy*dsy + dsz*dsz) > (All.MaxNodeMove*nop->len)*(All.MaxNodeMove*nop->len))
	    {
	      force_update_node(no, 1);

	      dx=nop->s[0] - tppos[0];     /* observe the sign ! */
	      dy=nop->s[1] - tppos[1];     /* this vector is -y in my thesis notation */
	      dz=nop->s[2] - tppos[2];
	    }
	  else
	    {
	      dx=nop->s[0] + dsx - tppos[0];     /* observe the sign ! */
	      dy=nop->s[1] + dsy - tppos[1];     /* this vector is -y in my thesis notation */
	      dz=nop->s[2] + dsz - tppos[2];
	    }

#ifdef PERIODIC
	  while(dx>BoxHalf) dx-=Box;
	  while(dy>BoxHalf) dy-=Box;
	  while(dz>BoxHalf) dz-=Box;
	  while(dx<-BoxHalf) dx+=Box;
	  while(dy<-BoxHalf) dy+=Box;
	  while(dz<-BoxHalf) dz+=Box;
#endif

	  r2=dx*dx+dy*dy+dz*dz;

	  if(nop->len*nop->len > r2*All.ErrTolTheta*All.ErrTolTheta)
	    {
	      no= nextnode[no];  /* open cell */
	      nop= &nodes[no];
	    }
	  else
	    {
#ifdef DIAG
	      treecost_quadru[tree]+=1;	  	  
#endif   
	      nop->cost += 1;

	      r=sqrt(r2);  
	  
	      u=r*h_inv;
	  
	      if(u>=1)  /* ordinary quadrupol moment */
		{
		  r_inv=1/r;
		  r2_inv=r_inv*r_inv;
		  r3_inv=r2_inv*r_inv;
		  r5_inv=r2_inv*r3_inv;

                  q11dx=nop->Q11*dx;
                  q12dy=nop->Q12*dy;
		  q13dz=nop->Q13*dz;
		  q12dx=nop->Q12*dx;
		  q22dy=nop->Q22*dy;
		  q23dz=nop->Q23*dz;
		  q13dx=nop->Q13*dx;
                  q23dy=nop->Q23*dy;
		  q33dz=nop->Q33*dz;

		  potq=0.5*(q11dx*dx+              /* 1/2 y^T Q y */
			    q22dy*dy+
			    q33dz*dz)+
		       q12dx*dy+
		       q13dx*dz+
		       q23dy*dz;
#ifdef POTENTIAL_SIMULTO_FORCE 
		  *tppot += -nop->mass*r_inv  /* monopole */
                           +r3_inv*( -3*potq*r2_inv + 0.5*nop->P);  /* quadrupole */
#endif
		  
		  fac=nop->mass*r3_inv  /* monopole force*/
		      +(15*potq*r2_inv -1.5*nop->P)*r5_inv;  /* radial quadrupole part */

		  tpaccel[0]+=dx*fac;
		  tpaccel[1]+=dy*fac;
		  tpaccel[2]+=dz*fac;

		  /* add tensor contribution */

		  ff=-3*r5_inv;
		  tpaccel[0] += ff*(q11dx + q12dy + q13dz);
		  tpaccel[1] += ff*(q12dx + q22dy + q23dz);
		  tpaccel[2] += ff*(q13dx + q23dy + q33dz);

		}
	      else    /* softened quadrupol moment */
		{
		  ii = (int)(u*KERN_LEN); ff=(u-knlrad[ii])*KERN_LEN;
		  wf=knlforce[ii] + (knlforce[ii+1]-knlforce[ii])*ff;
		  wp=knlpot[ii]   + (knlpot[ii+1]-knlpot[ii])*ff;
		  w2=knlW2[ii]    + (knlW2[ii+1]-knlW2[ii])*ff;
		  w3=knlW3[ii]    + (knlW3[ii+1]-knlW3[ii])*ff;
		  w4=knlW4[ii]    + (knlW4[ii+1]-knlW4[ii])*ff;

                  q11dx=nop->Q11*dx;
                  q12dy=nop->Q12*dy;
		  q13dz=nop->Q13*dz;
		  q12dx=nop->Q12*dx;
		  q22dy=nop->Q22*dy;
		  q23dz=nop->Q23*dz;
		  q13dx=nop->Q13*dx;
                  q23dy=nop->Q23*dy;
		  q33dz=nop->Q33*dz;

		  potq=0.5*(q11dx*dx+      /* 1/2 y^T Q y */
			    q22dy*dy+
			    q33dz*dz)+
		       q12dx*dy+
		       q13dx*dz+
		       q23dy*dz;
#ifdef POTENTIAL_SIMULTO_FORCE
		  *tppot += nop->mass*h_inv*wp +  /* monopole */   
		            potq*w2*h5_inv + 0.5*nop->P*wf*h2_inv*h_inv ; /* quadru contribution */
#endif
		  /* note: observe definition (sign!) of dx,dy,dz */

		  if(u>1.0e-4) 
		    {
		      r_inv=1/r;
		  
		      fac=nop->mass*h2_inv*h_inv*wf +  /* monopole force */ 
			+potq*h6_inv * w3*r_inv   + 0.5*nop->P * w4 *h4_inv*r_inv; /* radial contribution
									     of quadrupole */
		      tpaccel[0]+=dx*fac;
		      tpaccel[1]+=dy*fac;
		      tpaccel[2]+=dz*fac;
		  
		      /* add tensor contribution */
		      ff=w2*h5_inv;
		      tpaccel[0] += ff*(q11dx + q12dy + q13dz);
		      tpaccel[1] += ff*(q12dx + q22dy + q23dz);
		      tpaccel[2] += ff*(q13dx + q23dy + q33dz);
		    }
		}

#ifdef PERIODIC
	      ewald_corr(dx,dy,dz,fcorr);
	      
	      tpaccel[0]+= nop->mass*fcorr[0];
	      tpaccel[1]+= nop->mass*fcorr[1];
	      tpaccel[2]+= nop->mass*fcorr[2];
#endif

	      no=nop->sibling;   nop= &nodes[no];
	    }
	}
    }

}




/*  evaluate the force contribution of a tree based on the
 *  NEW RELATIVE OPENING CRITERION
 */
void force_treeevaluate_single(int tree, int targetpart, double epsilon) 
{
  int    no;
  struct NODE *nop;
  int    p,ii;
  double r2,dx,dy,dz,r,fac,u,h,h2_inv,h3_inv,h5_inv,h4_inv,h6_inv,ff;
  double wf,wp,w2,w3,w4,potq;
  float *tppos;   /* target particle for treewalk */
  double *tppot, *tpaccel;
  double q11dx,q12dy,q13dz,q12dx,q22dy,q23dz,q13dx,q23dy,q33dz;
  double r_inv,r2_inv,r3_inv,r5_inv;
  double h_inv;
  double dt,dsx,dsy,dsz;
  double oac;
#ifdef PERIODIC
  double fcorr[3];
#endif
  

  h = 2.8*epsilon;
  h_inv=1/h;
  h2_inv=h_inv*h_inv;
  h3_inv=h2_inv*h_inv;
  h4_inv=h2_inv*h2_inv;
  h5_inv=h2_inv*h3_inv;
  h6_inv=h3_inv*h3_inv;


  tppos   = &GravDataIn[targetpart].Pos[0];
  tpaccel = &GravDataResult[targetpart].Acc[0];
  tppot   = &GravDataPotential[targetpart];

  oac= GravDataIn[targetpart].OldAcc * All.ErrTolForceAcc;

  no=trees[tree];  nop= &nodes[no];

  while(no>=0)
    {
      if(no < All.MaxPart)   /* single particle */
	{
	  p= no;

	  /* we do the prediction right here */

	  dt=(All.Time - Part[p].CurrentTime) * S_of_a_inverse;  /*  we need here 1.0/S(a) */	  

	  /* observe the sign ! */
	  /* this vector is -y in my thesis notation */

	  /* the following explicit type-casts are done because the computation in 
             brackets is carried out in double precision, while tppos is float...
             we want to guarantee dx=0 if p is the particle at tppos */
	  dx=(float)(Part[p].Pos[0]+Part[p].Vel[0]*dt) - tppos[0];     
	  dy=(float)(Part[p].Pos[1]+Part[p].Vel[1]*dt) - tppos[1];     
	  dz=(float)(Part[p].Pos[2]+Part[p].Vel[2]*dt) - tppos[2];

#ifdef PERIODIC
	  while(dx>BoxHalf) dx-=Box;
	  while(dy>BoxHalf) dy-=Box;
	  while(dz>BoxHalf) dz-=Box;
	  while(dx<-BoxHalf) dx+=Box;
	  while(dy<-BoxHalf) dy+=Box;
	  while(dz<-BoxHalf) dz+=Box;
#endif
	  r2=dx*dx+dy*dy+dz*dz;

	  r=sqrt(r2);  
	   
	  u=r*h_inv;

	  Part[p].GravCost +=1; 

#ifdef DIAG
	  treecost[tree] +=1;	  	  
#endif

	  if(u>=1)
	    {
	      r_inv=1/r;

	      fac= Part[p].Mass*r_inv*r_inv*r_inv;
	  
	      tpaccel[0]+=dx*fac;
	      tpaccel[1]+=dy*fac;
	      tpaccel[2]+=dz*fac;
#ifdef POTENTIAL_SIMULTO_FORCE
	      *tppot-= Part[p].Mass*r_inv;
#endif
	    }
	  else
	    {
	      ii = (int)(u*KERN_LEN); ff=(u-knlrad[ii])*KERN_LEN;
	      wf=knlforce[ii]+(knlforce[ii+1]-knlforce[ii])*ff;
	      wp=knlpot[ii]+(knlpot[ii+1]-knlpot[ii])*ff;
	      
	      if(u>1.0e-4)
		{
		  fac= Part[p].Mass*h_inv*h_inv*h_inv*wf;
		  
		  tpaccel[0]+=dx*fac;
		  tpaccel[1]+=dy*fac;
		  tpaccel[2]+=dz*fac;
		}
#ifdef POTENTIAL_SIMULTO_FORCE
	      *tppot+= Part[p].Mass*h_inv*wp;
#endif
	    }

#ifdef PERIODIC
	  if(u>1.0e-4)
	    {
	      ewald_corr(dx,dy,dz,fcorr);
	      
	      tpaccel[0]+=  Part[p].Mass*fcorr[0];
	      tpaccel[1]+=  Part[p].Mass*fcorr[1];
	      tpaccel[2]+=  Part[p].Mass*fcorr[2];
	    }
#endif
	  no= nextnode[no];   nop= &nodes[no];
	}
      else
	{
	  dt=(All.Time - nop->tilu) * S_of_a_inverse;  /*  we need here 1.0/S(a) */	  

	  dsx=nop->vs[0]*dt;
	  dsy=nop->vs[1]*dt;
	  dsz=nop->vs[2]*dt;


	  if((dsx*dsx + dsy*dsy + dsz*dsz) > (All.MaxNodeMove*nop->len)*(All.MaxNodeMove*nop->len))
	    {
	      force_update_node(no, 1);

	      dx=nop->s[0] - tppos[0];     /* observe the sign ! */
	      dy=nop->s[1] - tppos[1];     /* this vector is -y in my thesis notation */
	      dz=nop->s[2] - tppos[2];
	    }
	  else
	    {
	      dx=nop->s[0] + dsx - tppos[0];     /* observe the sign ! */
	      dy=nop->s[1] + dsy - tppos[1];     /* this vector is -y in my thesis notation */
	      dz=nop->s[2] + dsz - tppos[2];
	    }

#ifdef PERIODIC
	  while(dx> BoxHalf) dx-=Box;
	  while(dy> BoxHalf) dy-=Box;
	  while(dz> BoxHalf) dz-=Box;
	  while(dx<-BoxHalf) dx+=Box;
	  while(dy<-BoxHalf) dy+=Box;
	  while(dz<-BoxHalf) dz+=Box;
#endif

	  r2=dx*dx+dy*dy+dz*dz;


#ifdef BMAX
	  if(nop->oc > oac*r2*r2*r2 || r2<nop->bmax2)
#else
	  if(nop->oc > oac*r2*r2*r2)
#endif
	    {
	      no= nextnode[no];  /* open cell */
	      nop= &nodes[no];
	    }
	  else
	    {
#ifdef DIAG
	      treecost_quadru[tree]+=1;	  	  
#endif   
	      nop->cost += 1;

	      r=sqrt(r2);  
	  
	      u=r*h_inv;
	  
	      if(u>=1)  /* ordinary quadrupol moment */
		{
		  r_inv=1/r;
		  r2_inv=r_inv*r_inv;
		  r3_inv=r2_inv*r_inv;
		  r5_inv=r2_inv*r3_inv;

                  q11dx=nop->Q11*dx;
                  q12dy=nop->Q12*dy;
		  q13dz=nop->Q13*dz;
		  q12dx=nop->Q12*dx;
		  q22dy=nop->Q22*dy;
		  q23dz=nop->Q23*dz;
		  q13dx=nop->Q13*dx;
                  q23dy=nop->Q23*dy;
		  q33dz=nop->Q33*dz;

		  potq=0.5*(q11dx*dx+              /* 1/2 y^T Q y */
			    q22dy*dy+
			    q33dz*dz)+
		       q12dx*dy+
		       q13dx*dz+
		       q23dy*dz;
#ifdef POTENTIAL_SIMULTO_FORCE 
		  *tppot += -nop->mass*r_inv  /* monopole */
                           +r3_inv*( -3*potq*r2_inv + 0.5*nop->P);  /* quadrupole */
#endif
		  
		  fac=nop->mass*r3_inv  /* monopole force*/
		      +(15*potq*r2_inv -1.5*nop->P)*r5_inv;  /* radial quadrupole part */

		  tpaccel[0]+=dx*fac;
		  tpaccel[1]+=dy*fac;
		  tpaccel[2]+=dz*fac;

		  /* add tensor contribution */

		  ff=-3*r5_inv;
		  tpaccel[0] += ff*(q11dx + q12dy + q13dz);
		  tpaccel[1] += ff*(q12dx + q22dy + q23dz);
		  tpaccel[2] += ff*(q13dx + q23dy + q33dz);

		}
	      else    /* softened quadrupol moment */
		{
		  ii = (int)(u*KERN_LEN); ff=(u-knlrad[ii])*KERN_LEN;
		  wf=knlforce[ii] + (knlforce[ii+1]-knlforce[ii])*ff;
		  wp=knlpot[ii]   + (knlpot[ii+1]-knlpot[ii])*ff;
		  w2=knlW2[ii]    + (knlW2[ii+1]-knlW2[ii])*ff;
		  w3=knlW3[ii]    + (knlW3[ii+1]-knlW3[ii])*ff;
		  w4=knlW4[ii]    + (knlW4[ii+1]-knlW4[ii])*ff;

                  q11dx=nop->Q11*dx;
                  q12dy=nop->Q12*dy;
		  q13dz=nop->Q13*dz;
		  q12dx=nop->Q12*dx;
		  q22dy=nop->Q22*dy;
		  q23dz=nop->Q23*dz;
		  q13dx=nop->Q13*dx;
                  q23dy=nop->Q23*dy;
		  q33dz=nop->Q33*dz;

		  potq=0.5*(q11dx*dx+      /* 1/2 y^T Q y */
			    q22dy*dy+
			    q33dz*dz)+
		       q12dx*dy+
		       q13dx*dz+
		       q23dy*dz;
#ifdef POTENTIAL_SIMULTO_FORCE
		  *tppot += nop->mass*h_inv*wp +  /* monopole */   
		            potq*w2*h5_inv + 0.5*nop->P*wf*h2_inv*h_inv ; /* quadru contribution */
#endif
		  /* note: observe definition (sign!) of dx,dy,dz */
		  if(u>1.0e-4) 
		    {
		      r_inv=1/r;

		      fac=nop->mass*h2_inv*h_inv*wf +  /* monopole force */ 
			+potq*h6_inv * w3*r_inv   + 0.5*nop->P * w4 *h4_inv*r_inv; /* radial contribution
										      of quadrupole */
		      tpaccel[0]+=dx*fac;
		      tpaccel[1]+=dy*fac;
		      tpaccel[2]+=dz*fac;
		  
		      /* add tensor contribution */
		      ff=w2*h5_inv;
		      tpaccel[0] += ff*(q11dx + q12dy + q13dz);
		      tpaccel[1] += ff*(q12dx + q22dy + q23dz);
		      tpaccel[2] += ff*(q13dx + q23dy + q33dz);
		    }
		}

#ifdef PERIODIC
	      ewald_corr(dx,dy,dz,fcorr);
	      
	      tpaccel[0]+=nop->mass*fcorr[0];
	      tpaccel[1]+=nop->mass*fcorr[1];
	      tpaccel[2]+=nop->mass*fcorr[2];
#endif
	      no= nop->sibling;    nop= &nodes[no];
	    }
	}
    }

}






/* This function computes the gravitational potential
 * for the specified particle in the communication
 * buffer. The contributions of all the separate trees
 * are added up.
 */
void force_treeevaluate_potential(int target) 
{
  int    tr;
  double epsilon;
  double dmax(double, double);


  GravDataPotential[target]=0;

  for(tr=0; tr<6; tr++)
    {
      if(numnodestree[tr])
	{
	  epsilon=dmax(All.SofteningTable[tr], All.SofteningTable[GravDataIn[target].Type]);

	  if(All.TypeOfOpeningCriterion==0 || GravDataIn[target].OldAcc==0)
	    force_treeevaluate_potential_single_BH(tr,target,epsilon);
	  else
	    force_treeevaluate_potential_single(tr,target,epsilon);
	}
    }
}



/* tree-walk for the gravitational potential using the normal 
 * BH opening criterion 
 */ 
void force_treeevaluate_potential_single_BH(int tree, int targetpart, double epsilon)  /* non-recursive walk */
{
  int    no;
  struct NODE *nop;
  int p,ii;
  double r2,dx,dy,dz,r,u,h,h2_inv,h3_inv,h5_inv,h4_inv,h6_inv,ff;
  double wf,wp,w2,potq;
  float  *tppos;   /* target particle for treewalk */
  double *tppot;
  double q11dx,q12dy,q13dz,q12dx,q22dy,q23dz,q13dx,q23dy,q33dz;
  double r_inv,r2_inv,r3_inv,r5_inv;
  double h_inv;

  h = 2.8*epsilon;
  h_inv=1/h;
  h2_inv=h_inv*h_inv;
  h3_inv=h2_inv*h_inv;
  h4_inv=h2_inv*h2_inv;
  h5_inv=h2_inv*h3_inv;
  h6_inv=h3_inv*h3_inv;


  tppos   = &GravDataIn[targetpart].Pos[0];
  tppot   = &GravDataPotential[targetpart];


  no= trees[tree];  nop= &nodes[no];

  while(no>=0)
    {
      if(no < All.MaxPart)   /* single particle */
	{
	  p=no;

	  dx= Part[p].PosPred[0]- tppos[0];     /* observe the sign ! */
	  dy= Part[p].PosPred[1]- tppos[1];     /* this vector is -y in my thesis notation */
	  dz= Part[p].PosPred[2]- tppos[2];
	  
#ifdef PERIODIC
	  while(dx> BoxHalf) dx-=Box;
	  while(dy> BoxHalf) dy-=Box;
	  while(dz> BoxHalf) dz-=Box;
	  while(dx<-BoxHalf) dx+=Box;
	  while(dy<-BoxHalf) dy+=Box;
	  while(dz<-BoxHalf) dz+=Box;
#endif
	  
	  r2=dx*dx+dy*dy+dz*dz;

	  r=sqrt(r2);  
	   
	  u=r*h_inv;

	  if(u>=1)
	    {
	      *tppot-= Part[p].Mass/r;
	    }
	  else
	    {
	      ii = (int)(u*KERN_LEN); ff=(u-knlrad[ii])*KERN_LEN;
	      wp=knlpot[ii]+(knlpot[ii+1]-knlpot[ii])*ff;
	      
	      *tppot+= Part[p].Mass*h_inv*wp;
	    }

#ifdef PERIODIC
	  *tppot+= Part[p].Mass * ewald_pot_corr(dx,dy,dz);
#endif
	  no= nextnode[no];  nop= &nodes[no];
	}
      else
	{
	  dx=nop->s[0]- tppos[0];     /* observe the sign ! */
	  dy=nop->s[1]- tppos[1];     /* this vector is -y in my thesis notation */
	  dz=nop->s[2]- tppos[2];
	  
#ifdef PERIODIC
	  while(dx> BoxHalf) dx-=Box;
	  while(dy> BoxHalf) dy-=Box;
	  while(dz> BoxHalf) dz-=Box;
	  while(dx<-BoxHalf) dx+=Box;
	  while(dy<-BoxHalf) dy+=Box;
	  while(dz<-BoxHalf) dz+=Box;
#endif
	  
	  r2=dx*dx+dy*dy+dz*dz;

	  if(nop->len*nop->len > r2*All.ErrTolTheta*All.ErrTolTheta)
	    {
	      no= nextnode[no];  /* open cell */
	      nop= &nodes[no];
	    }
	  else
	    {
	      r=sqrt(r2);  
	  
	      u=r*h_inv;
	  
	      if(u>=1)  /* ordinary quadrupol moment */
		{
		  r_inv=1/r;
		  r2_inv=r_inv*r_inv;
		  r3_inv=r2_inv*r_inv;
		  r5_inv=r2_inv*r3_inv;

                  q11dx=nop->Q11*dx;
                  q12dy=nop->Q12*dy;
		  q13dz=nop->Q13*dz;
		  q12dx=nop->Q12*dx;
		  q22dy=nop->Q22*dy;
		  q23dz=nop->Q23*dz;
		  q13dx=nop->Q13*dx;
                  q23dy=nop->Q23*dy;
		  q33dz=nop->Q33*dz;

		  potq=0.5*(q11dx*dx+              /* 1/2 y^T Q y */
			    q22dy*dy+
			    q33dz*dz)+
		       q12dx*dy+
		       q13dx*dz+
		       q23dy*dz;
 
		  *tppot += -nop->mass*r_inv  /* monopole */
                           +r3_inv*( -3*potq*r2_inv + 0.5*nop->P);  /* quadrupole */

		}
	      else    /* softened quadrupol moment */
		{
		  ii = (int)(u*KERN_LEN); ff=(u-knlrad[ii])*KERN_LEN;
		  wf=knlforce[ii] + (knlforce[ii+1]-knlforce[ii])*ff;
		  wp=knlpot[ii]   + (knlpot[ii+1]-knlpot[ii])*ff;
		  w2=knlW2[ii]    + (knlW2[ii+1]-knlW2[ii])*ff;

                  q11dx=nop->Q11*dx;
                  q12dy=nop->Q12*dy;
		  q13dz=nop->Q13*dz;
		  q12dx=nop->Q12*dx;
		  q22dy=nop->Q22*dy;
		  q23dz=nop->Q23*dz;
		  q13dx=nop->Q13*dx;
                  q23dy=nop->Q23*dy;
		  q33dz=nop->Q33*dz;

		  potq=0.5*(q11dx*dx+      /* 1/2 y^T Q y */
			    q22dy*dy+
			    q33dz*dz)+
		       q12dx*dy+
		       q13dx*dz+
		       q23dy*dz;

		  *tppot += nop->mass*h_inv*wp +  /* monopole */   
		            potq*w2*h5_inv + 0.5*nop->P*wf*h2_inv*h_inv ; /* quadru contribution */

		}
#ifdef PERIODIC
	      *tppot+= nop->mass * ewald_pot_corr(dx,dy,dz);
#endif
	      no=nop->sibling;      nop= &nodes[no];
	    }
	}
    }

}


/* tree-walk for the gravitational potential using the     
 * NEW opening criterion
 */
void force_treeevaluate_potential_single(int tree, int targetpart, double epsilon)  /* non-recursive walk */
{
  int no;
  struct NODE *nop;
  int p,ii;
  double r2,dx,dy,dz,r,u,h,h2_inv,h3_inv,h5_inv,h4_inv,h6_inv,ff;
  double wf,wp,w2,potq;
  float  *tppos;   /* target particle for treewalk */
  double *tppot;
  double q11dx,q12dy,q13dz,q12dx,q22dy,q23dz,q13dx,q23dy,q33dz;
  double r_inv,r2_inv,r3_inv,r5_inv;
  double h_inv;
  double oac;

  h = 2.8*epsilon;
  h_inv=1/h;
  h2_inv=h_inv*h_inv;
  h3_inv=h2_inv*h_inv;
  h4_inv=h2_inv*h2_inv;
  h5_inv=h2_inv*h3_inv;
  h6_inv=h3_inv*h3_inv;


  tppos   = &GravDataIn[targetpart].Pos[0];
  tppot   = &GravDataPotential[targetpart];


  oac= GravDataIn[targetpart].OldAcc * All.ErrTolForceAcc;  


  no=trees[tree];  nop= &nodes[no];

  while(no>=0)
    {
      if(no < All.MaxPart)   /* single particle */
	{
	  p=no;

	  dx=Part[p].PosPred[0]- tppos[0];     /* observe the sign ! */
	  dy=Part[p].PosPred[1]- tppos[1];     /* this vector is -y in my thesis notation */
	  dz=Part[p].PosPred[2]- tppos[2];
	  
#ifdef PERIODIC
	  while(dx> BoxHalf) dx-=Box;
	  while(dy> BoxHalf) dy-=Box;
	  while(dz> BoxHalf) dz-=Box;
	  while(dx<-BoxHalf) dx+=Box;
	  while(dy<-BoxHalf) dy+=Box;
	  while(dz<-BoxHalf) dz+=Box;
#endif
	  
	  r2=dx*dx+dy*dy+dz*dz;

	  r=sqrt(r2);  
	   
	  u=r*h_inv;

	  if(u>=1)
	    {
	      *tppot-= Part[p].Mass/r;
	    }
	  else
	    {
	      ii = (int)(u*KERN_LEN); ff=(u-knlrad[ii])*KERN_LEN;
	      wp=knlpot[ii]+(knlpot[ii+1]-knlpot[ii])*ff;
	      
	      *tppot+= Part[p].Mass*h_inv*wp;
	    }

#ifdef PERIODIC
	  *tppot+= Part[p].Mass * ewald_pot_corr(dx,dy,dz);
#endif

	  no= nextnode[no];  nop= &nodes[no];
	}
      else
	{
	  
	  dx=nop->s[0]- tppos[0];     /* observe the sign ! */
	  dy=nop->s[1]- tppos[1];     /* this vector is -y in my thesis notation */
	  dz=nop->s[2]- tppos[2];
	  
#ifdef PERIODIC
	  while(dx> BoxHalf) dx-=Box;
	  while(dy> BoxHalf) dy-=Box;
	  while(dz> BoxHalf) dz-=Box;
	  while(dx<-BoxHalf) dx+=Box;
	  while(dy<-BoxHalf) dy+=Box;
	  while(dz<-BoxHalf) dz+=Box;
#endif
	  
	  r2=dx*dx+dy*dy+dz*dz;
#ifdef BMAX
	  if(nop->oc > oac*r2*r2*r2 || r2<nop->bmax2)
#else
	  if(nop->oc > oac*r2*r2*r2)
#endif
	    {
	      no= nextnode[no];  /* open cell */
	      nop= &nodes[no];
	    }
	  else
	    {
	      r=sqrt(r2);  
	  
	      u=r*h_inv;
	  
	      if(u>=1)  /* ordinary quadrupol moment */
		{
		  r_inv=1/r;
		  r2_inv=r_inv*r_inv;
		  r3_inv=r2_inv*r_inv;
		  r5_inv=r2_inv*r3_inv;

                  q11dx=nop->Q11*dx;
                  q12dy=nop->Q12*dy;
		  q13dz=nop->Q13*dz;
		  q12dx=nop->Q12*dx;
		  q22dy=nop->Q22*dy;
		  q23dz=nop->Q23*dz;
		  q13dx=nop->Q13*dx;
                  q23dy=nop->Q23*dy;
		  q33dz=nop->Q33*dz;

		  potq=0.5*(q11dx*dx+              /* 1/2 y^T Q y */
			    q22dy*dy+
			    q33dz*dz)+
		       q12dx*dy+
		       q13dx*dz+
		       q23dy*dz;
 
		  *tppot += -nop->mass*r_inv  /* monopole */
                           +r3_inv*( -3*potq*r2_inv + 0.5*nop->P);  /* quadrupole */

		}
	      else    /* softened quadrupol moment */
		{
		  ii = (int)(u*KERN_LEN); ff=(u-knlrad[ii])*KERN_LEN;
		  wf=knlforce[ii] + (knlforce[ii+1]-knlforce[ii])*ff;
		  wp=knlpot[ii]   + (knlpot[ii+1]-knlpot[ii])*ff;
		  w2=knlW2[ii]    + (knlW2[ii+1]-knlW2[ii])*ff;

                  q11dx=nop->Q11*dx;
                  q12dy=nop->Q12*dy;
		  q13dz=nop->Q13*dz;
		  q12dx=nop->Q12*dx;
		  q22dy=nop->Q22*dy;
		  q23dz=nop->Q23*dz;
		  q13dx=nop->Q13*dx;
                  q23dy=nop->Q23*dy;
		  q33dz=nop->Q33*dz;

		  potq=0.5*(q11dx*dx+      /* 1/2 y^T Q y */
			    q22dy*dy+
			    q33dz*dz)+
		       q12dx*dy+
		       q13dx*dz+
		       q23dy*dz;

		  *tppot += nop->mass*h_inv*wp +  /* monopole */   
		            potq*w2*h5_inv + 0.5*nop->P*wf*h2_inv*h_inv ; /* quadru contribution */

		}
#ifdef PERIODIC
	      *tppot+= nop->mass * ewald_pot_corr(dx,dy,dz);
#endif
	      no=nop->sibling;   nop= &nodes[no];
	    }
	}
    }
}



/* this routine tabulates the kernel used for the softening
 * of the force and potential of monopole and quadrupole
 * moments
 */
void force_setkernel(void) 
{
  int i;
  double u;

  for(i=0;i<=KERN_LEN;i++)
    {
      u=((double)i)/KERN_LEN;

      knlrad[i] = u;

      if(u<=0.5)
	{
	  knlforce[i]=32 * (1.0/3 -6.0/5*pow(u,2) + pow(u,3));
	  knlpot[i]=16.0/3*pow(u,2)-48.0/5*pow(u,4)+32.0/5*pow(u,5)-14.0/5;

	  knlW2[i]= -384.0/5 +96.0*u;
	  knlW3[i]= 96.0;
	  knlW4[i]= 96.0/5*u*(5*u-4);
	}
      else
	{
	  knlforce[i]=64*(1.0/3 -3.0/4*u + 3.0/5*pow(u,2)-pow(u,3)/6) - 1.0/15/pow(u,3);
	  knlpot[i]=1.0/15/u +32.0/3*pow(u,2)-16.0*pow(u,3)+48.0/5*pow(u,4)-32.0/15*pow(u,5)-16.0/5;

	  knlW2[i]= 384.0/5 + 1/(5.0*pow(u,5)) -48.0/u -32*u;
	  knlW3[i]= -32-1/pow(u,6)+48/pow(u,2);
	  knlW4[i]= -48 +1/(5*pow(u,4)) +384.0/5*u -32*pow(u,2);
	}
    }
}


/* this function allocates memory used for storage of the tree
 * and auxiliary arrays for tree-walk and link-lists.
 */
void force_treeallocate(int maxnodes, int maxpart)  /* usually maxnodes=0.7*maxpart is sufficient */
{
  int bytes, allbytes=0;

  MaxNodes=maxnodes;

  if(!(nodes_base=malloc(bytes=(MaxNodes+1)*sizeof(struct NODE))))
    {
      printf("failed to allocate memory for %d tree-nodes (%d bytes).\n", MaxNodes, bytes);
      endrun(3);
    }
  allbytes+=bytes;

  if(!(next=malloc(bytes=maxpart*sizeof(int4byte))))
    {
      fprintf(stdout,"Failed to allocate %d spaces for 'next' array (%d bytes)\n", maxpart, bytes);
      exit(0);
    }
  allbytes+=bytes;

  if(!(nextnode=malloc(bytes=(maxpart+maxnodes)*sizeof(int4byte))))
    {
      fprintf(stdout,"Failed to allocate %d spaces for 'nextnode' array (%d bytes)\n", maxpart+maxnodes, bytes);
      exit(0);
    }
  allbytes+=bytes;


  if(!(father=malloc(bytes=(maxpart+maxnodes)*sizeof(int4byte))))
    {
      fprintf(stdout,"Failed to allocate %d spaces for 'father' array (%d bytes)\n", maxpart+maxnodes, bytes);
      exit(0);
    }
  allbytes+=bytes;

  if(ThisTask==0)
    printf("\nAllocated %g MByte for BH-tree.\n\n",allbytes/(1024.0*1024.0));

  force_setkernel();
}


/* free the allocated memory
 */
void force_treefree(void)
{
  free(father);
  free(nextnode);
  free(next);
  free(nodes_base);
}


/* set some counters to zero
 */
void force_resetcost(void)
{
  int tr;

  Num_nodeupdates= Num_nodeupdate_particles= 0;

  for(tr=0; tr<6; tr++) 
    treecost[tr]= treecost_quadru[tr]= 0;
}

/* return the total number of particle-particle 
 * interactions that have occured
 */
int force_getcost_single(void)
{
  int tr,cost;

  for(tr=cost=0;tr<6;tr++) 
    cost += treecost[tr];
  
  return cost;
}

/* return the total number of particle-node (internal nodes) 
 * interactions that have occured
 */
int force_getcost_quadru(void)
{
  int tr,cost;

  for(tr=cost=0;tr<6;tr++) 
    cost += treecost_quadru[tr];
  
  return cost;
}



/* do force computation with direct summation 
 * for the specified particle in the communication buffer
 * (useful for debugging).
 */
void force_treeevaluate_direct(int target, double one_over_s_of_a) 
{
  double epsilon;
  double h, h_inv, dx, dy, dz, r, r2, u, wf, ff, r_inv, fac, dt;
  int    i,ii;
#ifdef PERIODIC
  double fcorr[3];
#endif

  S_of_a_inverse= one_over_s_of_a; 
 
  for(i=0;i<3;i++)
    GravDataResult[target].Acc[i]=0;
  
  for(i=0; i<NumPart; i++)
    {
      epsilon= dmax(All.SofteningTable[Part[i].Type], All.SofteningTable[GravDataIn[target].Type]);
      
      h = 2.8*epsilon;
      h_inv=1/h;

      dt=(All.Time - Part[i].CurrentTime) * S_of_a_inverse;

      dx= Part[i].Pos[0] + Part[i].Vel[0]*dt - GravDataIn[target].Pos[0];
      dy= Part[i].Pos[1] + Part[i].Vel[1]*dt - GravDataIn[target].Pos[1];
      dz= Part[i].Pos[2] + Part[i].Vel[2]*dt - GravDataIn[target].Pos[2];

#ifdef PERIODIC
      while(dx>BoxHalf) dx-=Box;
      while(dy>BoxHalf) dy-=Box;
      while(dz>BoxHalf) dz-=Box;
      while(dx<-BoxHalf) dx+=Box;
      while(dy<-BoxHalf) dy+=Box;
      while(dz<-BoxHalf) dz+=Box;
#endif
 
      r2=dx*dx+dy*dy+dz*dz;

      r=sqrt(r2);  
	   
      u=r*h_inv;

      if(u>=1)
	{
	  r_inv=1/r;
	  
	  fac=Part[i].Mass*r_inv*r_inv*r_inv;
	  
          GravDataResult[target].Acc[0]+= dx*fac;
	  GravDataResult[target].Acc[1]+= dy*fac;
	  GravDataResult[target].Acc[2]+= dz*fac;
	}
      else
	{
	  ii = (int)(u*KERN_LEN); ff=(u-knlrad[ii])*KERN_LEN;
	  wf=knlforce[ii]+(knlforce[ii+1]-knlforce[ii])*ff;
	      
	  if(u>1.0e-4)
	    {
	      fac=Part[i].Mass*h_inv*h_inv*h_inv*wf;
	      
	      GravDataResult[target].Acc[0]+= dx*fac;
	      GravDataResult[target].Acc[1]+= dy*fac;
	      GravDataResult[target].Acc[2]+= dz*fac;
	    }
	}

#ifdef PERIODIC
      if(u>1.0e-4)
	{
	  ewald_corr(dx,dy,dz,fcorr);
	  
	  GravDataResult[target].Acc[0]+=  Part[i].Mass*fcorr[0];
	  GravDataResult[target].Acc[1]+=  Part[i].Mass*fcorr[1];
	  GravDataResult[target].Acc[2]+=  Part[i].Mass*fcorr[2];
	}
#endif

    }
}





/***************************************   
 *    Below we have the routines       *
 *    dealing with neighbour finding.  *
 *                                     *
 *    We use the gravity tree and a    *
 *    range-searching technique to     *
 *    find neighbours.                 *
 ***************************************/


static int    *ngblist, numngb;
static float  *r2list;
static float  searchmin[3], searchmax[3], searchcenter[3];
static float  h_i, h_i2;


/* this routine maps a coordinate difference
 * to the nearest periodic image
 */
float INLINE_FUNC ngb_periodic(float x)
{
  while(x > BoxHalf)
    x-=Box;
  while(x < -BoxHalf)
      x+=Box;
  return x;
}



/* This routine finds all neighbours `j' that can interact with the particle
 * `i' in the communication buffer. 
 *
 * Note that an interaction takes place if:   r_ij < h_i  OR  r_ij < h_j. 
 * 
 * In the range-search this is taken into account. For this purpose, each
 * node knows the maximum h occuring among the particles inside it. 
 *
 */
int ngb_treefind_pairs(float xyz[3], float hsml, int **ngblistback, float **r2listback)
{
  int    i,ind,k;
  float dx,dy,dz,r2;

  h_i = hsml;
  h_i2= h_i*h_i;
  
  for(k=0; k<3; k++) /* cube-box window */
    {
      searchcenter[k]=xyz[k];
      searchmin[k]=xyz[k]- hsml;
      searchmax[k]=xyz[k]+ hsml;
    }
      

  numngb=0;
  ngb_treesearch_pairs(trees[0]);  /* gas particles always have a separate tree */

  if(numngb==MAX_NGB)
    {
      printf("PROBLEM! In ngb_treefind_pairs() we reached MAX_NGB=%d ... STOPPING\n", MAX_NGB);
      endrun(77);
    }

  for(i=0; i<numngb; i++)
    {
      ind=ngblist[i];
      dx=Part[ind].PosPred[0] - xyz[0];
      dy=Part[ind].PosPred[1] - xyz[1];
      dz=Part[ind].PosPred[2] - xyz[2];

#ifdef PERIODIC
      dx= ngb_periodic(dx);
      dy= ngb_periodic(dy);
      dz= ngb_periodic(dz);
#endif
      r2=dx*dx+dy*dy+dz*dz;
      r2list[i]=r2;
    }

  *ngblistback= ngblist;
  *r2listback= r2list;


  return numngb;
}


/* recursive search routine for the function ngb_treefind_pairs()
 */
void ngb_treesearch_pairs(int no)
{
  int k,p;
  float hdiff;
  struct NODE *this;


  if(no < All.MaxPart) /* single particle */
    {
      for(k=0, p= no; k<3; k++)
	{
	  hdiff= SphPart[p].Hsml - h_i;
	  if(hdiff<0)
	    hdiff=0;
#ifdef PERIODIC
	  if(ngb_periodic(Part[p].PosPred[k] - searchcenter[k]) < (searchmin[k]-hdiff-searchcenter[k]))
	    return;
	  if(ngb_periodic(Part[p].PosPred[k] - searchcenter[k]) > (searchmax[k]+hdiff-searchcenter[k]))
	    return;
#else
	  if(Part[p].PosPred[k]<(searchmin[k]-hdiff))
	    return;
	  if(Part[p].PosPred[k]>(searchmax[k]+hdiff))
	    return;
#endif
	}

      if(numngb<MAX_NGB)
	ngblist[numngb++]= p;
    }
  else
    {
      this= &nodes[no];
      hdiff= this->hmax - h_i;
      if(hdiff<0)
	hdiff=0;

      for(k=0; k<3; k++)
	{
#ifdef PERIODIC
	  if((ngb_periodic(this->center[k]- searchcenter[k])+0.5*this->len) < (searchmin[k]-hdiff-searchcenter[k]))
	    return;
	  if((ngb_periodic(this->center[k]- searchcenter[k])-0.5*this->len) > (searchmax[k]+hdiff-searchcenter[k]))
	    return;
#else
	  if((this->center[k]+0.5*this->len)<(searchmin[k]-hdiff))
	    return;
	  if((this->center[k]-0.5*this->len)>(searchmax[k]+hdiff))
	    return;
#endif
	}

      for(k=0; k<3; k++)
	{
#ifdef PERIODIC
	  if((ngb_periodic(this->center[k]- searchcenter[k])+0.5*this->len) > (searchmax[k]-searchcenter[k]))
	    break;
	  if((ngb_periodic(this->center[k]- searchcenter[k])-0.5*this->len) < (searchmin[k]-searchcenter[k]))
	    break;
#else
	  if((this->center[k]+0.5*this->len)>searchmax[k])
	    break;
	  if((this->center[k]-0.5*this->len)<searchmin[k])
	    break;
#endif
	}

      if(k>=3) 	  /* node lies completely inside */
	{
	  for(k=0, p= this->partind; k<this->count; k++, p= next[p])
	    {
	      if(numngb<MAX_NGB)
		ngblist[numngb++]=p;
	    }
	}
      else
	{
	  for(k=0;k<8;k++) 
	    if((no=this->suns[k])>=0)
	      {
		ngb_treesearch_pairs(no);  /* open cell */
	      }
	}
    }
}




/*  ngb_treefind_variable() returns all neighbours (and only those) with distance <= hguess 
 *  and returns them in ngblistback and r2listback
 */
int ngb_treefind_variable(float xyz[3], float hguess, int parttype, int **ngblistback, float **r2listback)
{
  float sr, sr2;  /* search radius */
  int   i, ind, k;
  float dx, dy, dz, r2;
  

  sr=hguess;
  
  for(k=0;k<3;k++) /* cube-box window */
    {
      searchcenter[k]= xyz[k];
      searchmin[k]= xyz[k]-sr;
      searchmax[k]= xyz[k]+sr;
    }
      
  sr2= sr*sr;
  numngb= 0;
  ngb_treesearch(trees[parttype]);  


  if(numngb==MAX_NGB)
    {
      printf("PROBLEM! In ngb_treefind_variable() we reached MAX_NGB=%d ... STOPPING\n", MAX_NGB);
      endrun(78);
    }


  for(i=0; i<numngb; i++)
    {
      ind= ngblist[i];
      
      dx= Part[ind].PosPred[0] - xyz[0];
      dy= Part[ind].PosPred[1] - xyz[1];
      dz= Part[ind].PosPred[2] - xyz[2];
#ifdef PERIODIC
      dx= ngb_periodic(dx);
      dy= ngb_periodic(dy);
      dz= ngb_periodic(dz);
#endif
      r2= dx*dx+dy*dy+dz*dz;
      r2list[i]=r2;

      if(r2>=sr2)
	{
	  ngblist[i]= ngblist[numngb-1];
	  numngb--;
	  i--;
	}
    }

  *ngblistback= ngblist;
  *r2listback= r2list;

  return numngb;
}



/* recursive neighbour search for the routine ngb_treefind_variable().
 */
void ngb_treesearch(int no)
{
  int k, p;
  struct NODE *this;

  if(no<All.MaxPart)  /* single particle */
    {
      for(k=0, p=no; k<3; k++)
	{
#ifdef PERIODIC
	  if(ngb_periodic(Part[p].PosPred[k]-searchcenter[k]) < (searchmin[k]-searchcenter[k]))
	    return;
	  if(ngb_periodic(Part[p].PosPred[k]-searchcenter[k]) > (searchmax[k]-searchcenter[k]))
	    return;
#else
	  if(Part[p].PosPred[k]<searchmin[k])
	    return;
	  if(Part[p].PosPred[k]>searchmax[k])
	    return;
#endif
	}
      if(numngb<MAX_NGB)
	ngblist[numngb++]= p;
    }
  else
    {
      this= &nodes[no];

      for(k=0;k<3;k++)
	{
#ifdef PERIODIC
	  if((ngb_periodic(this->center[k]-searchcenter[k])+0.5*this->len) < (searchmin[k]-searchcenter[k]))
	    return;
	  if((ngb_periodic(this->center[k]-searchcenter[k])-0.5*this->len) > (searchmax[k]-searchcenter[k]))
	    return;
#else
	  if((this->center[k]+0.5*this->len)<searchmin[k])
	    return;
	  if((this->center[k]-0.5*this->len)>searchmax[k])
	    return;
#endif
	}

      for(k=0; k<3; k++)
	{
#ifdef PERIODIC
	  if((ngb_periodic(this->center[k]-searchcenter[k])+0.5*this->len) > (searchmax[k]-searchcenter[k]))
	    break;
	  if((ngb_periodic(this->center[k]-searchcenter[k])-0.5*this->len) < (searchmin[k]-searchcenter[k]))
	    break;
#else
	  if((this->center[k]+0.5*this->len)>searchmax[k])
	    break;
	  if((this->center[k]-0.5*this->len)<searchmin[k])
	    break;
#endif
	}

      if(k>=3) 	  /* cell lies completely inside */
	{
	  for(k=0, p=this->partind; k<this->count; k++, p= next[p])
	    if(numngb<MAX_NGB)
	      ngblist[numngb++]=p;
	}
      else
	{
	  for(k=0;k<8;k++) 
	    if((no=this->suns[k])>=0)
	      {
		ngb_treesearch(no);  /* open cell */
	      }
	}
    }
}


/* this routine finds exactly `desngb' neighbours on the local domain.
 * Ideally, the starting smoothing length is slightly larger than the 
 * smoothing radius actually needed, but the routine will find the 
 * correct soothing length in any case.
 *
 * For `hguess==', the routine estimates the local density by walking the tree until a
 * fair guess can be obtained.
 *
 * The square of the desired smoothing radius for `desngb' neighbours
 * is returned, and the actual neighbours are in the lists ngblistback and r2listback.
 */
float ngb_treefind(float xyz[3], int desngb, float hguess, int parttype, int **ngblistback, float **r2listback)
{
#define  SECFACTOR  1.2
#ifndef  PI
#define  PI               3.1415927
#endif
  float sr,sr2,h2max;  /* search radius */
  int   i,ind,j,subnode,fak,k,rep=0;
  float dx,dy,dz,r2;
  int   th, nnp;


  if(hguess>0)
    sr=hguess;
  else
    {
      /* determine estimate of local density */
      th= trees[parttype];
      while(nodes[th].count>200)
	{
	  for(j=0,subnode=0,fak=1;j<3;j++,fak<<=1)
	    if(xyz[j] > nodes[th].center[j])
	      subnode +=fak;
	  
	  if(nodes[th].suns[subnode]>=All.MaxPart)
	    {
	      nnp= nodes[th].suns[subnode];
	      if(nodes[nnp].count>200)
		th=nnp;
	      else
		break;
	    }
	  else
	    break;
	}

      sr= nodes[th].len * pow((3.0/(4*PI)*SECFACTOR)*desngb/((float)(nodes[th].count)),1.0/3);
    }

  do
    {
      for(k=0; k<3; k++)
	{
	  searchcenter[k]=xyz[k];
	  searchmin[k]=xyz[k]-sr;
	  searchmax[k]=xyz[k]+sr;
	}
      
      sr2=sr*sr;
      numngb=0;

      ngb_treesearch(trees[parttype]);  rep++;

      if(numngb==MAX_NGB)
	{
	  printf("PROBLEM! In ngb_treefind() we reached MAX_NGB=%d ... STOPPING\n", MAX_NGB);
	  endrun(79);
	}


      if(numngb<desngb)
	{
	  if(numngb>5)
	    sr*=pow((2.1*(float)desngb)/numngb,1.0/3);
	  else
	    sr*=2.0;

	  continue;
	}

      for(i=0;i<numngb;i++)
	{
	  ind=ngblist[i];
	  dx=Part[ind].PosPred[0]-xyz[0];
	  dy=Part[ind].PosPred[1]-xyz[1];
	  dz=Part[ind].PosPred[2]-xyz[2];
#ifdef PERIODIC
	  dx= ngb_periodic(dx);
	  dy= ngb_periodic(dy);
	  dz= ngb_periodic(dz);
#endif
	  r2=dx*dx+dy*dy+dz*dz;

	  r2list[i]= r2;
	}
      

      h2max= ngb_select_closest(desngb, numngb, r2list-1, ngblist-1);

      if(h2max<=sr2 && numngb>=desngb) 
	break;

      sr*=1.26;       /* 3th root of 2.0 */

      continue;
    }
  while(1);


  *ngblistback=ngblist;
  *r2listback=r2list;

  return h2max;
}







/* Allocate memory for the neighbour lists
 */
void ngb_treeallocate(int npart)  
{
  int totbytes=0, bytes;

  if(!(ngblist= malloc(bytes= npart*(long)sizeof(int))))
    {
      fprintf(stdout,"Failed to allocate %d bytes for ngblist array\n", bytes);
      endrun(78);
    }
  totbytes+= bytes;


  if(!(r2list= malloc(bytes= npart*(long)sizeof(float))))
    {
      fprintf(stdout,"Failed to allocate %d bytes for r2list array\n", bytes);
      endrun(78);
    }
  totbytes+= bytes;

  if(ThisTask==0)
    printf("allocated %f Mbyte for ngb search.\n", ((float)totbytes)/(1024.0*1024.0));
}


/* free memory allocated for neighbour lists
 */
void ngb_treefree(void)
{
  free(r2list);
  free(ngblist);
}

/* To construct the neighbour tree,
 * we actually need to construct the force-tree,
 * because we use it now for the neighbour search.
 * This routine is obsolute at this point.
 */
void ngb_treebuild(void) 
{
  if(ThisTask==0)
      printf("Begin Ngb-tree construction.\n");

  force_treebuild();

  if(ThisTask==0)
    printf("Ngb-Tree contruction finished \n");
}




/* this routine finds the maximum extent of the gas particle 
 * distribution on the local processors, and it updates the
 * sizes `len' of all nodes of the gas-tree such that they
 * really encompass all the particles of the node.
 * This is done to guarantee that the neighbour search will
 * always be correct, even if the tree is not reconstructed.
 * For the same reason, the maximum SPH smoothing for the 
 * particles represented by the node is redetermined.
 * Note: All SPH particles are guaranteed to be predicted to the
 * the current time when this function is called.
 */
void ngb_update_nodes(void) 
{
  int    i, j, fa, N;
  float  extmax;
  float  rel;

  /* 
   * Note: DomainMax/Min is set by the tree-construction, and is then shrunk
   * by determine_interior. Here it is again increased to the (new) enclosing box.
   */

  N= N_gas;  
#ifdef VELDISP
  N= NumPart;
#endif
#ifdef SIDM
  N= NumPart;
#endif

  for(i=0; i<N; i++)
    {
      for(j=0; j<3; j++)
	{
	  if(Part[i].PosPred[j] > DomainMax[Part[i].Type][j]) 
	    DomainMax[Part[i].Type][j]= Part[i].PosPred[j];
	  if(Part[i].PosPred[j] < DomainMin[Part[i].Type][j]) 
	    DomainMin[Part[i].Type][j]= Part[i].PosPred[j];
	}

      fa= father[i];  

      for(j=0, extmax=0; j<3; j++)
	{
	  rel=  Part[i].PosPred[j] - nodes[fa].center[j];
	  if(rel<0) rel= -rel;
	  
	  if(extmax < rel)
	    extmax = rel;
	}
      
      extmax *= 2;

      if(i<N_gas)
	{
	  if(extmax > nodes[fa].len || SphPart[i].Hsml > nodes[fa].hmax)
	    {
	      if(extmax > nodes[fa].len)
		nodes[fa].len= extmax;
	      if(SphPart[i].Hsml > nodes[fa].hmax)
		nodes[fa].hmax= SphPart[i].Hsml;
	      
	      force_update_size_of_parent_node(fa);
	    }
	}
      else
	{
	  if(extmax > nodes[fa].len)
	    {
	      nodes[fa].len= extmax;

	      force_update_size_of_parent_node(fa);
	    }
	}
    }
}




void dump_particles(void)
{
  FILE *fd;
  char buffer[200];
  int  i, j;
  double  dt, dt_h0, s_a, s_a_inverse;

  sprintf(buffer, "particles%d.dat", ThisTask);
  fd= fopen(buffer, "w");
  my_fwrite(&NumPart, 1, sizeof(int), fd);

  for(i=1; i<=NumPart; i++)
    my_fwrite(&P[i].PosPred[0], 3, sizeof(float), fd);

  for(i=1; i<=NumPart; i++)
    my_fwrite(&P[i].Pos[0], 3, sizeof(float), fd);

  for(i=1; i<=NumPart; i++)
    my_fwrite(&P[i].VelPred[0], 3, sizeof(float), fd);

  for(i=1; i<=NumPart; i++)
    my_fwrite(&P[i].Vel[0], 3, sizeof(float), fd);

  for(i=1; i<=NumPart; i++)
    my_fwrite(&P[i].CurrentTime, 1, sizeof(float), fd);

  for(i=1; i<=NumPart; i++)
    my_fwrite(&P[i].MaxPredTime, 1, sizeof(float), fd);

  for(i=1; i<=NumPart; i++)
    my_fwrite(&P[i].ID, 1, sizeof(int), fd);

  s_a=All.Hubble*sqrt(All.Omega0 + All.Time*(1-All.Omega0-All.OmegaLambda) 
		                 + All.Time*All.Time*All.Time*All.OmegaLambda);
  s_a_inverse=1/s_a;

  for(i=1; i<=NumPart; i++)
    {
      dt = All.Time - P[i].CurrentTime; 
      dt_h0 = dt*s_a_inverse;  
      
      for(j=0;j<3;j++)  
	{
	  P[i].PosPred[j] = P[i].Pos[j] + P[i].Vel[j]*dt_h0;
	  P[i].VelPred[j] = P[i].Vel[j] + P[i].Accel[j]*dt;
	}
    }
  
  for(i=1; i<=NumPart; i++)
    my_fwrite(&P[i].PosPred[0], 3, sizeof(float), fd);

  for(i=1; i<=NumPart; i++)
    my_fwrite(&P[i].VelPred[0], 3, sizeof(float), fd);
    
  
  my_fwrite(&s_a_inverse, 1, sizeof(double), fd);


  fclose(fd);


  printf("All.Time= %g\n", All.Time);
  printf("s_a_inverse= %g\n", s_a_inverse);
}






float ngb_select_closest(int k, int n, float *arr, int *ind)
{
#define SWAP(a,b)  temp =(a);(a)=(b);(b)=temp;
#define SWAPI(a,b) tempi=(a);(a)=(b);(b)=tempi;
  int   i, ir, j, l, mid, ai, tempi;
  float a,temp;
  
  l=1;
  ir=n;
  while(1)
    {
      if(ir <= l+1)
	{
	  if (ir == l+1 && arr[ir] < arr[l]) 
	    {
	      SWAP(arr[l],arr[ir]);
	      SWAPI(ind[l],ind[ir]);
	    }
	  return arr[k];
	} 
      else 
	{
	  mid=(l+ir) >> 1;
	  SWAP(arr[mid],arr[l+1]);
	  SWAPI(ind[mid],ind[l+1]);

	  if (arr[l] > arr[ir]) 
	    {
	      SWAP(arr[l],arr[ir]);
	      SWAPI(ind[l],ind[ir]);
	    }
	  if(arr[l+1] > arr[ir]) 
	    {
	      SWAP(arr[l+1],arr[ir]);
	      SWAPI(ind[l+1],ind[ir]);
	    }
	  if(arr[l] > arr[l+1])
	    {
	      SWAP(arr[l],arr[l+1]);
	      SWAPI(ind[l],ind[l+1]);
	    }
	  i=l+1;
	  j=ir;
	  a=arr[l+1];
	  ai=ind[l+1];

	  while(1)
	    {
	      do 
		i++; 
	      while(arr[i] < a);
	      
	      do 
		j--; 
	      while(arr[j] > a);
			
	      if(j < i) 
		break;
	      SWAP(arr[i],arr[j]);
	      SWAPI(ind[i],ind[j]);
	    }

	  arr[l+1]=arr[j];
	  arr[j]=a;
	  ind[l+1]=ind[j];
	  ind[j]=ai;

	  if(j >= k) 
	    ir=j-1;
	  if(j <= k) 
	    l=i;
	}
    }
#undef SWAP
#undef SWAPI
}





// SIDM
#ifdef SIDM
void update_node_of_scat_particle(int i) 
{
  // index i= particle number 0,1,2... P[i+1]

  i= father[i];
	  
  while(i>=0){
    force_update_node(i, 0 /*flag Prediction*/);
    // Assume prediction in gravtree.c every timestep.
    nodes[i].tilu= All.Time;
    if(nodes[i].count > 1000)
      break;
    i= father[i]; 
  }   
}
#endif











