/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/fof/fof_findgroups.c
 * \date        MM/YYYY
 * \author
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <gsl/gsl_math.h>
#include <inttypes.h>

#include "../allvars.h"
#include "../proto.h"
#include "../domain.h"
#include "fof.h"
#include "../subfind/subfind.h"



#ifdef FOF


#ifdef GENERIC_ASYNC
#undef GENERIC_ASYNC            /* for the moment, this function cannot yet use the asynchronous schem */
#endif


static int fof_find_dmparticles_evaluate(int target, int mode, int threadid);
static int fof_treefind_fof_primary(MyDouble searchcenter[3], MyFloat hsml, int target, int numnodes, int *firstnode, int mode, int threadid);

static int *Tree_Head;

static MyIDType *MinID;
static int *Head, *Len, *Next, *Tail, *MinIDTask;


/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];

  MyIDType MinID;
  int MinIDTask;

  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;


/* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in * in, int i, int firstnode)
{
  in->Pos[0] = P[i].Pos[0];
  in->Pos[1] = P[i].Pos[1];
  in->Pos[2] = P[i].Pos[2];

  in->MinID = MinID[Head[i]];
  in->MinIDTask = MinIDTask[Head[i]];

  in->Firstnode = firstnode;
}

 /* local data structure that holds results acquired on remote processors */
typedef struct
{
  char link_count_flag;
} data_out;

static data_out *DataResult, *DataOut;


 /* routine to store or combine result data */
static void out2particle(data_out * out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES)      /* initial store */
    {
      terminate("here not used");
    }
  else                          /* combine */
    {
      if(out->link_count_flag)
        Flags[i].Marked = 1;
    }

}

#include "../generic_comm_helpers2.h"

static int link_across;
static int nprocessed;

static void kernel_local(void)
{
  int i;
#ifdef GENERIC_ASYNC
  int flag = 0;
#endif
  /* do local particles */
#pragma omp parallel private(i) reduction(+:nprocessed)
  {
    int j, threadid = get_thread_num();
#ifdef GENERIC_ASYNC
    int count = 0;
#endif

    for(j = 0; j < NTask; j++)
      Thread[threadid].Exportflag[j] = -1;

    while(1)
      {
        if(Thread[threadid].ExportSpace < MinSpace)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              if(generic_polling_primary(count, NumPart))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        i = NextParticle++;

        if(i >= NumPart)
          break;

        if(((1 << P[i].Type) & (FOF_PRIMARY_LINK_TYPES)))
          {
            if(Flags[i].Nonlocal && Flags[i].Changed)
              {
                fof_find_dmparticles_evaluate(i, MODE_LOCAL_PARTICLES, threadid);

                nprocessed++;
              }
          }
      }
  }
}

static void kernel_imported(void)
{
  /* now do the particles that were sent to us */
  int i, cnt = 0;
#pragma omp parallel private(i) reduction(+:link_across)
  {
    int threadid = get_thread_num();
#ifdef GENERIC_ASYNC
    int count = 0;
#endif

    while(1)
      {
#pragma omp atomic capture
        i = cnt++;

        if(i >= Nimport)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              generic_polling_secondary();
          }

        count++;
#endif

        link_across += fof_find_dmparticles_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }

}

double fof_find_groups(MyIDType * vMinID, int *vHead, int *vLen, int *vNext, int *vTail, int *vMinIDTask)
{
  MinID = vMinID;
  Head = vHead;
  Len = vLen;
  Next = vNext;
  Tail = vTail;
  MinIDTask = vMinIDTask;


  int i, npart, marked;
  long long totmarked, totnpart;
  long long link_across_tot, ntot;
  double t0, t1, tstart, tend;

  tstart = second();

  mpi_printf("FOF: Start linking particles (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

  /* allocate a flag field that is used to mark nodes that are fully inside the linking length */
  flag_node_inside_linkinglength = (unsigned char *) mymalloc("flag_node_inside_linkinglength", Tree_MaxNodes * sizeof(unsigned char));
  memset(flag_node_inside_linkinglength, 0, Tree_MaxNodes * sizeof(unsigned char));
  flag_node_inside_linkinglength -= Tree_MaxPart;

  Flags = (struct bit_flags *) mymalloc("Flags", NumPart * sizeof(struct bit_flags));

  generic_set_MaxNexport();


  Tree_Head = mymalloc("Tree_Head", Tree_NumNodes * sizeof(int));
  Tree_Head -= Tree_MaxPart;


  /* allocate buffers to arrange communication */
  generic_alloc_partlist_nodelist_ngblist_threadbufs();

  t0 = second();

  /* first, link only among local particles */
  for(i = 0, marked = 0, npart = 0; i < NumPart; i++)
    {
      if(((1 << P[i].Type) & (FOF_PRIMARY_LINK_TYPES)))
        {
          fof_find_dmparticles_evaluate(i, MODE_LOCAL_NO_EXPORT, 0);

          npart++;

          if(Flags[i].Nonlocal)
            marked++;
        }
    }

  sumup_large_ints(1, &marked, &totmarked);
  sumup_large_ints(1, &npart, &totnpart);
  t1 = second();
  mpi_printf("FOF: links on local processor done (took %g sec).\nFOF: Marked=%lld out of the %lld primaries which are linked\n", timediff(t0, t1), totmarked, totnpart);

  generic_free_partlist_nodelist_ngblist_threadbufs();


  t0 = second();
  fof_check_for_full_nodes_recursive(Tree_MaxPart);
  t1 = second();
  mpi_printf("FOF: fully linked nodes determined (took %g sec).\n", timediff(t0, t1));


  mpi_printf("FOF: begin linking across processors (presently allocated=%g MB) \n", AllocatedBytes / (1024.0 * 1024.0));

  for(i = 0; i < NumPart; i++)
    Flags[i].Marked = 1;

  do
    {
      t0 = second();

      for(i = 0; i < NumPart; i++)
        {
          Flags[i].Changed = Flags[i].Marked;
          Flags[i].Marked = 0;
          Flags[i].MinIDChanged = 0;
        }

      NextParticle = 0;         /* begin with this index */

      link_across = 0;
      nprocessed = 0;

      generic_comm_pattern(NumPart, kernel_local, kernel_imported);

      sumup_large_ints(1, &link_across, &link_across_tot);
      sumup_large_ints(1, &nprocessed, &ntot);

      t1 = second();

      mpi_printf("FOF: have done %15lld cross links (processed %14lld, took %g sec)\n", link_across_tot, ntot, timediff(t0, t1));

      /* let's check out which particles have changed their MinID */
      for(i = 0; i < NumPart; i++)
        if(Flags[i].Nonlocal)
          {
            if(Flags[Head[i]].MinIDChanged)
              Flags[i].Marked = 1;
          }

    }
  while(link_across_tot > 0);

  Tree_Head += Tree_MaxPart;
  myfree(Tree_Head);
  myfree(Flags);
  /* free flag */
  myfree(flag_node_inside_linkinglength + Tree_MaxPart);

  mpi_printf("FOF: Local groups found.\n");

  tend = second();
  return timediff(tstart, tend);
}


static int fof_find_dmparticles_evaluate(int target, int mode, int threadid)
{
  int j, n, links, p, s, ss, numnodes, *firstnode;
  int numngb;
  MyDouble *pos;
  data_in local, *target_data;

  links = 0;

  if(mode == MODE_LOCAL_NO_EXPORT || mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      target_data = &local;

      numnodes = 1;
      firstnode = NULL;
    }
  else
    {
      target_data = &DataGet[target];

      generic_get_numnodes(target, &numnodes, &firstnode);
    }


  pos = target_data->Pos;

  numngb = fof_treefind_fof_primary(pos, LinkL, target, numnodes, firstnode, mode, threadid);


  if(mode == MODE_LOCAL_PARTICLES || mode == MODE_LOCAL_NO_EXPORT)
    for(n = 0; n < numngb; n++)
      {
        j = Thread[threadid].Ngblist[n];

        if(Head[target] != Head[j])     /* only if not yet linked */
          {
            if(Len[Head[target]] > Len[Head[j]])        /* p group is longer */
              {
                p = target;
                s = j;
              }
            else
              {
                p = j;
                s = target;
              }
            Next[Tail[Head[p]]] = Head[s];

            Tail[Head[p]] = Tail[Head[s]];

            Len[Head[p]] += Len[Head[s]];

            if(MinID[Head[s]] < MinID[Head[p]])
              {
                MinID[Head[p]] = MinID[Head[s]];
                MinIDTask[Head[p]] = MinIDTask[Head[s]];
              }

            ss = Head[s];
            do
              Head[ss] = Head[p];
            while((ss = Next[ss]) >= 0);
          }
      }

  if(mode == MODE_IMPORTED_PARTICLES)
    {
      if(numngb > 0)
        DataResult[target].link_count_flag = 1;
      else
        DataResult[target].link_count_flag = 0;
    }

  links += numngb;

  return links;
}


/*! This function finds the neighbors among the primary link types which are within a certain distance.
 *
 * @param searchcenter
 * @param hsml
 * @param target
 * @param startnode
 * @param mode
 *          -1: only local particles should be found and no export occurs
 *          0:  export occurs, but local particles are ignored
 *          1:  particles are found for an imported point
 * @param nexport
 * @param nsend_local
 * @return
 */
static int fof_treefind_fof_primary(MyDouble searchcenter[3], MyFloat hsml, int target, int numnodes, int *firstnode, int mode, int threadid)
{
  int k, numngb, no, p, nexport_flag = 0;
  MyDouble dx, dy, dz, dist, r2;

#define FACT2 0.866025403785    /* sqrt(3)/2 */
#define FACT3 (2.0*FACT2)       /* sqrt(3)   */
#ifdef PERIODIC
  MyDouble xtmp, ytmp, ztmp;
#endif

  numngb = 0;

  for(k = 0; k < numnodes; k++)
    {
      if(mode == MODE_LOCAL_PARTICLES || mode == MODE_LOCAL_NO_EXPORT)
        {
          no = Tree_MaxPart;    /* root node */
        }
      else
        {
          no = firstnode[k];
          no = Nodes[no].u.d.nextnode;  /* open it */
        }

      while(no >= 0)
        {
          if(no < Tree_MaxPart) /* single particle */
            {
              p = no;
              no = Nextnode[no];

              if(!((1 << P[p].Type) & (FOF_PRIMARY_LINK_TYPES)))
                continue;

              if(mode == MODE_LOCAL_PARTICLES)
                continue;

              dist = hsml;
              dx = FOF_NEAREST_LONG_X(Tree_Pos_list[3 * p + 0] - searchcenter[0]);
              if(dx > dist)
                continue;
              dy = FOF_NEAREST_LONG_Y(Tree_Pos_list[3 * p + 1] - searchcenter[1]);
              if(dy > dist)
                continue;
              dz = FOF_NEAREST_LONG_Z(Tree_Pos_list[3 * p + 2] - searchcenter[2]);
              if(dz > dist)
                continue;
              if(dx * dx + dy * dy + dz * dz > dist * dist)
                continue;

              if(mode == MODE_IMPORTED_PARTICLES)
                {
                  if(MinID[Head[p]] > DataGet[target].MinID)
                    {
                      MinID[Head[p]] = DataGet[target].MinID;
                      MinIDTask[Head[p]] = DataGet[target].MinIDTask;
                      Flags[Head[p]].MinIDChanged = 1;
                      numngb++;
                    }
                }
              else
                {
                  /* this will only be done for MODE_LOCAL_NO_EXPORT */
                  Thread[threadid].Ngblist[numngb++] = p;
                }
            }
          else if(no < Tree_MaxPart + Tree_MaxNodes)    /* internal node */
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                {
                  if(no < Tree_FirstNonTopLevelNode)    /* we reached a top-level node again, which means that we are done with the branch */
                    break;

                  if(Tree_Head[no] >= 0)
                    if(MinID[Tree_Head[no]] <= DataGet[target].MinID)
                      {
                        no = Nodes[no].u.d.sibling;     /* the node can be discarded */
                        continue;
                      }
                }

              struct NODE *current = &Nodes[no];
              int nocur = no;
              no = current->u.d.sibling;        /* in case the node can be discarded */

              if(mode == MODE_LOCAL_PARTICLES)
                {
                  if(nocur >= Tree_FirstNonTopLevelNode)
                    {
                      /* we have a node with only local particles, hence we can skip it for mode == 0 */
                      continue;
                    }
                }

              dist = hsml + 0.5 * current->len;
              dx = FOF_NEAREST_LONG_X(current->center[0] - searchcenter[0]);
              if(dx > dist)
                continue;
              dy = FOF_NEAREST_LONG_Y(current->center[1] - searchcenter[1]);
              if(dy > dist)
                continue;
              dz = FOF_NEAREST_LONG_Z(current->center[2] - searchcenter[2]);
              if(dz > dist)
                continue;

              /* now test against the minimal sphere enclosing everything */
              dist += FACT1 * current->len;
              r2 = dx * dx + dy * dy + dz * dz;
              if(r2 > dist * dist)
                continue;

              if(mode != MODE_LOCAL_PARTICLES)
                {
                  /* test whether the node is contained within the sphere */
                  dist = hsml - FACT2 * current->len;
                  if(dist > 0)
                    if(r2 < dist * dist && hsml > FACT3 * current->len)
                      {
                        if(flag_node_inside_linkinglength[nocur] & (1 << BITFLAG_INSIDE_LINKINGLENGTH)) /* already flagged */
                          {
                            /* sufficient to return only one particle inside this cell */
                            p = fof_return_a_particle_in_cell_recursive(nocur);

                            if(p >= 0)
                              {
                                if(mode == MODE_IMPORTED_PARTICLES)
                                  {
                                    if(MinID[Head[p]] > DataGet[target].MinID)
                                      {
                                        MinID[Head[p]] = DataGet[target].MinID;
                                        MinIDTask[Head[p]] = DataGet[target].MinIDTask;
                                        Flags[Head[p]].MinIDChanged = 1;
                                        numngb++;
                                      }
                                  }
                                else
                                  Thread[threadid].Ngblist[numngb++] = p;
                              }

                            continue;
                          }
                        else
                          {
                            /* flag it now */
                            flag_node_inside_linkinglength[nocur] |= (1 << BITFLAG_INSIDE_LINKINGLENGTH);
                          }
                      }
                }

              no = current->u.d.nextnode;       /* ok, we need to open the node */
            }
          else if(no >= Tree_ImportedNodeOffset)        /* point from imported nodelist */
            {
              terminate("do not expect imported points here");
            }
          else
            {
              if(mode == MODE_LOCAL_PARTICLES)
                {
                  if(target >= 0)
                    tree_treefind_export_node_threads(no, target, threadid);
                }
              else if(mode == MODE_LOCAL_NO_EXPORT)
                {
                  nexport_flag = 1;
                }
              else if(mode == MODE_IMPORTED_PARTICLES)
                terminate("stop no=%d Tree_MaxPart=%d Tree_MaxNodes=%d", no, Tree_MaxPart, Tree_MaxNodes);

              no = Nextnode[no - Tree_MaxNodes];
              continue;
            }
        }
    }

  if(mode == MODE_LOCAL_NO_EXPORT)
    {
      if(nexport_flag == 0)
        Flags[target].Nonlocal = 0;
      else
        Flags[target].Nonlocal = 1;
    }

  return numngb;
}



void fof_check_for_full_nodes_recursive(int no)
{
  if(no >= Tree_MaxPart && no < Tree_MaxPart + Tree_MaxNodes)   /* internal node */
    {
      int head = -1;            /* no particle yet */

      int p = Nodes[no].u.d.nextnode;

      while(p != Nodes[no].u.d.sibling)
        {
          if(p < Tree_MaxPart)  /* a particle */
            {
              if((1 << P[p].Type) & (FOF_PRIMARY_LINK_TYPES))
                {
                  if(head == -1)
                    head = Head[p];
                  else if(head >= 0)
                    {
                      if(head != Head[p])
                        head = -2;
                    }
                }

              p = Nextnode[p];
            }
          else if(p < Tree_MaxPart + Tree_MaxNodes)     /* an internal node  */
            {
              fof_check_for_full_nodes_recursive(p);

              if(head == -1)
                head = Tree_Head[p];
              else if(head >= 0)
                {
                  if(head != Tree_Head[p])
                    head = -2;
                }

              p = Nodes[p].u.d.sibling;
            }
          else                  /* a pseudo particle */
            p = Nextnode[p - Tree_MaxNodes];
        }

      Tree_Head[no] = head;
    }
}



int fof_return_a_particle_in_cell_recursive(int no)
{
  if(no >= Tree_MaxPart && no < Tree_MaxPart + Tree_MaxNodes)   /* internal node */
    {
      int p = Nodes[no].u.d.nextnode;

      while(p != Nodes[no].u.d.sibling)
        {
          if(p < Tree_MaxPart)  /* a particle */
            {
              if((1 << P[p].Type) & (FOF_PRIMARY_LINK_TYPES))
                {
                  return p;
                }

              p = Nextnode[p];
            }
          else if(p < Tree_MaxPart + Tree_MaxNodes)     /* an internal node  */
            {
              int ret = fof_return_a_particle_in_cell_recursive(p);

              if(ret >= 0)
                return ret;

              p = Nodes[p].u.d.sibling;
            }
          else                  /* a pseudo particle */
            p = Nextnode[p - Tree_MaxNodes];
        }
    }

  return -1;
}





#endif
