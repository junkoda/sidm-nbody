#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"



static int    node_tail;
static double endofstrip;


/* increases Time to smallest max-prediction time and
 * determines which particles are grouped together
 * for force evaluation 
 */
void find_next_time(void)
{
  int i,node,count;
  double min, minGlobal, min_endofstrip;
  double t0, t1;
  
  t0=second();

  /* find the minimum in the tree */

  node=TimeTreeRoot;

  while(PTimeTree[node].left)
    node=PTimeTree[node].left;

  min=P[node].MaxPredTime;

  /* determine global minimum and broadcast it to all processors */
  MPI_Allreduce(&min, &minGlobal, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  All.TimeStep= minGlobal - All.Time;
  All.Time=     minGlobal;


  /* set-up a link-list of the particles in ascending order of
     MaxPredTime, in a strip subject to the condition that the 
     particle would be advanced at least half its timestep */

  NumForceUpdate=0;
  IndFirstUpdate=0;
  endofstrip=MAX_REAL_NUMBER;

  find_next_time_walk(TimeTreeRoot);


  /* now tailor the lists to a common strip */

  MPI_Allreduce(&endofstrip, &min_endofstrip, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  for(i=IndFirstUpdate, count=0, NumSphUpdate=0; count<NumForceUpdate; i=P[i].ForceFlag, count++)   
    {
      if(P[i].MaxPredTime > min_endofstrip)
	{
	  NumForceUpdate= count;
	  break;
	}
      if(P[i].Type==0)
	NumSphUpdate++;

      node_tail=i;
    }

  if(IndFirstUpdate)
    P[node_tail].ForceFlag= IndFirstUpdate; /* this ensures that ForceFlag!=0 for the last particle as well */

  t1=second();
  
  All.CPU_TimeLine+= timediff(t0,t1);
}


/* This routine walks the tree of the timeline, retrieving 
 * particles in the order of their maximum prediction time.
 * The walk is terminated when a particle is reached that can't
 * fo at least half its timestep if we continue.
 */
int find_next_time_walk(int node)
{
  if(PTimeTree[node].left)
    if(find_next_time_walk(PTimeTree[node].left))
      return 1;
  
  if((P[node].MaxPredTime - All.Time) <= 0.5*(P[node].MaxPredTime-P[node].CurrentTime) )
    {
      if(IndFirstUpdate==0)
	{
	  IndFirstUpdate=node;
	  node_tail=node;
	}
      else
	{
	  P[node_tail].ForceFlag= node;
	  node_tail=node;
	}

      NumForceUpdate++;
    }
  else
    {
      endofstrip = P[node].MaxPredTime; 

      return 1;   /* terminate tree walk */
    }

  if(PTimeTree[node].right)
    if(find_next_time_walk(PTimeTree[node].right))
      return 1;

  return 0;
}





/* This routine contructs a binary tree that contains the
 * particles ordered by their maximum prediction time.
 */
void construct_timetree(void)
{
  int i;
  int current_node,new_node;

  /* set-up root node */
  /* note: it is assumed that there is at least one particle */
  PTimeTree[1].left= PTimeTree[1].right=0;
  TimeTreeRoot=1;


  /* put in the other particles */
  for(i=2;i<=NumPart;i++)
    {
      current_node=TimeTreeRoot;
     
      do
	{
	  if(P[i].MaxPredTime < P[current_node].MaxPredTime)
	    {
	      if((new_node=PTimeTree[current_node].left))
		{
		  current_node=new_node;
		}
	      else
		{
		  PTimeTree[current_node].left= i; 
		  PTimeTree[i].left=PTimeTree[i].right=0;
		  break;
		}
	    }
	  else
	    {
	      if((new_node=PTimeTree[current_node].right))
		{
		  current_node=new_node;
		}
	      else
		{
		  PTimeTree[current_node].right= i; 
		  PTimeTree[i].left=PTimeTree[i].right=0;
		  break;
		}
	    }
	}
      while(new_node);
    }

#ifdef DEB
  if(ThisTask==0)
    { 	
      FdDEB=fopen("construct.deb","w");
      fprintf(FdDEB," %d \n",All.NumCurrentTiStep);
      fclose(FdDEB);
    }
#endif
}








/* This function inserts a new particle into the tree of the
 * timeline. This will be used in find_timesteps() to 
 * update the timesteps of particles.   
 */
void insert_node(int i)
{
  int current_node,new_node;


  current_node=TimeTreeRoot;
     
  do
    {
      if(P[i].MaxPredTime < P[current_node].MaxPredTime)
	{
	  if((new_node=PTimeTree[current_node].left))
	    {
	      current_node=new_node;
	    }
	  else
	    {
	      PTimeTree[current_node].left= i; 
	      PTimeTree[i].left=PTimeTree[i].right=0;
	      break;
	    }
	}
      else
	{
	  if((new_node=PTimeTree[current_node].right))
	    {
	      current_node=new_node;
	    }
	  else
	    {
	      PTimeTree[current_node].right= i; 
	      PTimeTree[i].left=PTimeTree[i].right=0;
	      break;
	    }
	}
      
    }
  while(new_node);
}


/*  This function delets a new particle from the tree of the
 *  timeline. 
 *  Note that at the time of the call of this function in find_timesteps()
 *  the variable MaxPredTime can still be used to find the way 
 *  to the particle along the tree.
 */
void delete_node(int i)
{
  int an,node,annode;

  an= find_ancestor(i);      /*  we need to find the ancestor.
			        an=0 if 'i' is the root node */

  if(PTimeTree[i].left>0 && PTimeTree[i].right>0)
    {
      /* ok, let's find the smallest node on the right side, and 
        also determine its ancestor */
      
      node=PTimeTree[i].right;
      annode=i;

      while(PTimeTree[node].left)
	{
	  annode=node;
	  node=PTimeTree[node].left;
	}
      
      /* cut the node out */
      if(PTimeTree[annode].left == node)
	PTimeTree[annode].left=PTimeTree[node].right;
      else
	PTimeTree[annode].right=PTimeTree[node].right;


      /* now put node 'node' instead of 'i' */
      
      if(an)
	{
	  if(PTimeTree[an].left == i)
	    PTimeTree[an].left=node;
	  else
	    PTimeTree[an].right=node;
	}
      else
	TimeTreeRoot= node;

      PTimeTree[node].left = PTimeTree[i].left;
      PTimeTree[node].right= PTimeTree[i].right;
    }
  else
    {
      if(PTimeTree[i].left) 
	{
	  if(an)
	    {
	      if(PTimeTree[an].left == i)
		PTimeTree[an].left=PTimeTree[i].left;
	      else
		PTimeTree[an].right=PTimeTree[i].left;
	    }
	  else
	    TimeTreeRoot=PTimeTree[i].left;
	}
      else
	{
	  if(PTimeTree[i].right)
	    {
	      if(an)
		{
		  if(PTimeTree[an].left == i)
		    PTimeTree[an].left=PTimeTree[i].right;
		  else
		    PTimeTree[an].right=PTimeTree[i].right;
		}
	      else
		TimeTreeRoot=PTimeTree[i].right;
	    }
	  else   /* a leaf */
	    {

	      /* note: we presume that this is not the root, i.e.
                 the tree will never be completely empty */ 
	      if(PTimeTree[an].left == i)
		PTimeTree[an].left=0;
	      else
		PTimeTree[an].right=0;
	    }
	}
    }
}



/*  This function finds the ancestor node for a
 *  particle in the timeline. 
 */
int find_ancestor(int i)
{
  int current_node, new_node, ancestor_node;

  if(i==TimeTreeRoot)   /* node is the root . no ancestor */
    return 0;         

  current_node=TimeTreeRoot;
  
  do
    {
      if(P[i].MaxPredTime < P[current_node].MaxPredTime)
	{
	  new_node=PTimeTree[current_node].left;
	}
      else
	{
	  new_node=PTimeTree[current_node].right;
	}
		 
      ancestor_node=current_node;
      current_node=new_node;
    }
  while(current_node != i);

  return ancestor_node;
}













