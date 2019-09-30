/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/modgrav/modgrav_forcetree.c
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
#include <time.h>


#include "../allvars.h"
#include "../proto.h"
#include "../domain.h"

/*! Constructs the full gravitational oct-tree down to All.MaxLevelFullTree;\n
 *  necessary for the AMR solver.
 *
 *  \return if successful returns 0\n
 *          -1 if the number of allocated tree nodes is too small
 */

int modgrav_force_construct_full_tree()
{
  int parent = -1, redblack;
  MyFloat lenhalf;


  if(All.MaxLevelFullTree > All.MinLevelTopLeaf)        /* add empty local nodes down to level specified in MaxLevelFullTree */
    {
      int nextfirstlocal = Tree_NextFreeNode;
      for(int i = 0; i < NTopleaves; i++)
        {
          parent = DomainNodeIndex[i];

          if(DomainTask[i] == ThisTask)
            {
              for(int j = 0; j < 8; j++)    /* create 8 subnodes */
                {
                  Nodes[parent].u.suns[j] = Tree_NextFreeNode;
                  struct NODE *nfreep = &Nodes[Tree_NextFreeNode];

                  nfreep->len = 0.5 * Nodes[parent].len;
                  lenhalf = 0.25 * Nodes[parent].len;

                  redblack = 0;

                  if(j & 1)
                    {
                      nfreep->center[0] = Nodes[parent].center[0] + lenhalf;
                      ++redblack;
                    }
                  else
                    nfreep->center[0] = Nodes[parent].center[0] - lenhalf;

                  if(j & 2)
                    {
                      nfreep->center[1] = Nodes[parent].center[1] + lenhalf;
                      ++redblack;
                    }
                  else
                    nfreep->center[1] = Nodes[parent].center[1] - lenhalf;

                  if(j & 4)
                    {
                      nfreep->center[2] = Nodes[parent].center[2] + lenhalf;
                      ++redblack;
                    }
                  else
                    nfreep->center[2] = Nodes[parent].center[2] - lenhalf;

                  nfreep->mg_bitflag.red = 1 - (redblack % 2);

                  nfreep->u.suns[0] = -1;
                  nfreep->u.suns[1] = -1;
                  nfreep->u.suns[2] = -1;
                  nfreep->u.suns[3] = -1;
                  nfreep->u.suns[4] = -1;
                  nfreep->u.suns[5] = -1;
                  nfreep->u.suns[6] = -1;
                  nfreep->u.suns[7] = -1;

                  Tree_NumNodes++;
                  Tree_NextFreeNode++;

                  if((Tree_NumNodes) >= Tree_MaxNodes)
                    {
                      if(All.TreeAllocFactor > MAX_TREE_ALLOC_FACTOR)
                        {
                          printf("task %d: cannot fit local empty nodes\n", ThisTask);
                          terminate("Something wrong\n");
                        }
                      else
                        {
                          return -1;
                        }
                    }
                }
            }
        }

      int firstlocal, lastlocal = Tree_NextFreeNode - 1;
      for(int curlevel = All.MinLevelTopLeaf + 1; curlevel <= All.MaxLevelFullTree - 1; ++curlevel)
        {
          firstlocal = nextfirstlocal;
          nextfirstlocal = Tree_NextFreeNode;
          for(parent = firstlocal; parent <= lastlocal; ++parent)
            {
              for(int j = 0; j < 8; j++)    /* create 8 subnodes */
                {
                  Nodes[parent].u.suns[j] = Tree_NextFreeNode;
                  struct NODE *nfreep = &Nodes[Tree_NextFreeNode];

                  nfreep->len = 0.5 * Nodes[parent].len;
                  lenhalf = 0.25 * Nodes[parent].len;

                  redblack = 0;

                  if(j & 1)
                    {
                      nfreep->center[0] = Nodes[parent].center[0] + lenhalf;
                      ++redblack;
                    }
                  else
                    nfreep->center[0] = Nodes[parent].center[0] - lenhalf;

                  if(j & 2)
                    {
                      nfreep->center[1] = Nodes[parent].center[1] + lenhalf;
                      ++redblack;
                    }
                  else
                    nfreep->center[1] = Nodes[parent].center[1] - lenhalf;

                  if(j & 4)
                    {
                      nfreep->center[2] = Nodes[parent].center[2] + lenhalf;
                      ++redblack;
                    }
                  else
                    nfreep->center[2] = Nodes[parent].center[2] - lenhalf;

                  nfreep->mg_bitflag.red = 1 - (redblack % 2);

                  nfreep->u.suns[0] = -1;
                  nfreep->u.suns[1] = -1;
                  nfreep->u.suns[2] = -1;
                  nfreep->u.suns[3] = -1;
                  nfreep->u.suns[4] = -1;
                  nfreep->u.suns[5] = -1;
                  nfreep->u.suns[6] = -1;
                  nfreep->u.suns[7] = -1;

                  Tree_NumNodes++;
                  Tree_NextFreeNode++;

                  if(Tree_NumNodes >= Tree_MaxNodes)
                    {
                      if(All.TreeAllocFactor > MAX_TREE_ALLOC_FACTOR)
                        {
                          printf("task %d: cannot fit local empty nodes\n", ThisTask);
                          terminate("Something wrong\n");
                        }
                      else
                        {
                          return -1;
                        }
                    }
                }
            }
          lastlocal = Tree_NextFreeNode - 1;
        }
    }
  return 0;
}

#ifdef MODGRAV
/*! inserts a single particle into the gravitational tree\n
 *  in contrast to the standard version of this function all 8 subnodes are constructed once a node is opened in this modgrav version
 *
 *  \return 0 if successful \n
 *          -1 if too few nodes have been allocated in the Nodes array
 */

int force_treebuild_insert_single_point(int i /*!< index of particle */ ,
                                        unsigned long long *intpos /*!< integer representation of particle position */ ,
                                        int th /*!< target node */ ,
                                        unsigned char levels /*!< level of target node */ )
{
  int j, parent = -1;
  unsigned char subnode = 0;
  unsigned long long xxb = intpos[0];
  unsigned long long yyb = intpos[1];
  unsigned long long zzb = intpos[2];
  unsigned long long mask = ((unsigned long long) 1) << ((52 - 1) - levels);
  unsigned char shiftx = (52 - 1) - levels;
  unsigned char shifty = (52 - 2) - levels;
  unsigned char shiftz = (52 - 3) - levels;
  signed long long centermask = (0xFFF0000000000000llu);
  unsigned long long *intppos;
  centermask >>= levels;


  while(1)
    {
      if(th >= Tree_MaxPart && th < Tree_ImportedNodeOffset)    /* we are dealing with an internal node */
        {
          subnode = (((unsigned char) ((xxb & mask) >> (shiftx--))) | ((unsigned char) ((yyb & mask) >> (shifty--))) | ((unsigned char) ((zzb & mask) >> (shiftz--))));

          centermask >>= 1;
          mask >>= 1;
          levels++;

          if(levels > MAX_TREE_LEVEL)
            {
              /* seems like we're dealing with particles at identical (or extremely close)
               * locations. Shift subnode index to allow tree construction. Note: Multipole moments
               * of tree are still correct, but one should MAX_TREE_LEVEL large enough to have
               *      DomainLen/2^MAX_TREE_LEEL  < gravitational softening length
               */
              for(j = 0; j < 8; j++)
                {
                  if(Nodes[th].u.suns[subnode] < 0)
                    break;

                  subnode++;
                  if(subnode >= 8)
                    subnode = 7;
                }
            }

          int nn = Nodes[th].u.suns[subnode];

          if(nn >= 0)           /* ok, something is in the daughter slot already, need to continue */
            {
              parent = th;
              th = nn;
            }
          else
            {
              /* here we have found an empty slot where we can attach
               * the new particle as a leaf.
               */
              Nodes[th].u.suns[subnode] = i;
              break;            /* done for this particle */
            }
        }
      else
        {
          /* We try to insert into a leaf with a single particle.  Need
           * to generate a new internal node at this point.
           */

          int old_subnode = subnode; /* store the subnode where the particle was attached */

          for(int j = 0; j <8; j++) /*create 8 subnodes */
            {
              if(Nodes[parent].u.suns[j] < All.MaxPart)
                {
                  int attached_particle;

                  if(Nodes[parent].u.suns[j] >= 0)
                    attached_particle = Nodes[parent].u.suns[j];
                  else
                    attached_particle = -1;

                  Nodes[parent].u.suns[j] = Tree_NextFreeNode;
                  struct NODE *nfreep = &Nodes[Tree_NextFreeNode];

                  /* one possibility is:
                            double len = 2 * ((force_int_to_double(mask) - 1.0) * DomainLen);
                            double cx = (force_int_to_double((xxb & centermask) | mask) - 1.0) * DomainLen + DomainCorner[0];
                            double cy = (force_int_to_double((yyb & centermask) | mask) - 1.0) * DomainLen + DomainCorner[1];
                            double cz = (force_int_to_double((zzb & centermask) | mask) - 1.0) * DomainLen + DomainCorner[2];
                   */

                  /* the other is:
                            double len = ((double) (mask << 1)) * DomainBigFac;
                            double cx = ((double) ((xxb & centermask) | mask)) * DomainBigFac + DomainCorner[0];
                            double cy = ((double) ((yyb & centermask) | mask)) * DomainBigFac + DomainCorner[1];
                            double cz = ((double) ((zzb & centermask) | mask)) * DomainBigFac + DomainCorner[2];
                   */

                  /* Can we use either of the methods above here? Let's stick with the conventional version first. */
                  double len = 0.5 * Nodes[parent].len;
                  double lenhalf  = 0.5 * len;
                  int redblack = 0;

                  double cx, cy, cz;
                  if(j & 1)
                    {
                      cx = Nodes[parent].center[0] + lenhalf;
                      ++redblack;
                    }
                  else
                    cx = Nodes[parent].center[0] - lenhalf;

                  if(j & 2)
                    {
                      cy = Nodes[parent].center[1] + lenhalf;
                      ++redblack;
                    }
                  else
                    cy = Nodes[parent].center[1] - lenhalf;

                  if(j & 4)
                    {
                      cz = Nodes[parent].center[2] + lenhalf;
                      ++redblack;
                    }
                  else
                    cz = Nodes[parent].center[2] - lenhalf;

                  nfreep->mg_bitflag.red = 1 - (redblack % 2);

                  nfreep->len = len;
                  nfreep->center[0] = cx;
                  nfreep->center[1] = cy;
                  nfreep->center[2] = cz;

                  for(j = 0; j < 8; j++)
                    nfreep->u.suns[j] = -1;

                  if(attached_particle >= 0)
                    {
                      if(th >= Tree_ImportedNodeOffset)
                        intppos = Tree_Points[th - Tree_ImportedNodeOffset].IntPos;
                      else
                        intppos = &Tree_IntPos_list[3 * th];

                      subnode = (((unsigned char) ((intppos[0] & mask) >> shiftx)) | ((unsigned char) ((intppos[1] & mask) >> shifty)) | ((unsigned char) ((intppos[2] & mask) >> shiftz)));

                      nfreep->u.suns[subnode] = th;
                    }

                  if(j == old_subnode)
                    th = Tree_NextFreeNode;       /* resume trying to insert the new particle the newly created internal node */

                  Tree_NumNodes++;
                  Tree_NextFreeNode++;

                  if(Tree_NumNodes >= Tree_MaxNodes)
                    {
                      if(All.TreeAllocFactor > MAX_TREE_ALLOC_FACTOR)
                        {
                          char buf[500];
                          sprintf(buf,
                                  "task %d: looks like a serious problem for particle %d, stopping with particle dump.  Tree_NumNodes=%d Tree_MaxNodes=%d  Tree_NumPartImported=%d NumPart=%d\n",
                                  ThisTask, i, Tree_NumNodes, Tree_MaxNodes, Tree_NumPartImported, NumPart);
                          dump_particles();
                          terminate(buf);
                        }
                      return -1;
                    }
                }

              else
                terminate("should not happen! Why is there already a node?");
            }
        }
    }

  return 0;
}
#endif

#ifdef MODGRAV_EFF_MASS
/*! computes the modified gravity acceleration due to
 * the effective mass of the node no on the given target particle.\n
 */
void modgrav_tree_fifth_force(double pos_x, double pos_y, double pos_z, int no, double asmthfac, double hmin, float * shortrange_table, double * acc_x, double * acc_y, double * acc_z
#ifdef EVALPOTENTIAL
                 , MyLongDouble * pot
#endif
  )
{
  struct NODE *nop = &Nodes[no];

  double r2, r, dx, dy, dz, mass, fac, h, h_inv, h3_inv, u;

#if (defined(PERIODIC) && !defined(GRAVITY_NOT_PERIODIC)) || defined(GRAVITY_TALLBOX)
  double xtmp, ytmp, ztmp;
#endif
#ifdef EVALPOTENTIAL
  double facpot, wp;
#endif


  dx = GRAVITY_NEAREST_X(nop->mg.eff.s_eff[0] - pos_x);
  dy = GRAVITY_NEAREST_Y(nop->mg.eff.s_eff[1] - pos_y);
  dz = GRAVITY_NEAREST_Z(nop->mg.eff.s_eff[2] - pos_z);

  r2 = dx * dx + dy * dy + dz * dz;
  r = sqrt(r2);

  mass = nop->mg.eff.eff_mass;

  h = 1.0 * nop->len;      /* set softening to multiple of node size */
  if(h < hmin)
    h = hmin;

  if(r >= h)
    {
      fac = mass / (r2 * r);
#ifdef EVALPOTENTIAL
      facpot = -mass / r;
#endif
    }
  else
    {
      h_inv = 1.0 / h;
      h3_inv = h_inv * h_inv * h_inv;
      u = r * h_inv;

      if(u < 0.5)
        fac = mass * h3_inv * (10.666666666667 + u * u * (32.0 * u - 38.4));
      else
        fac = mass * h3_inv * (21.333333333333 - 48.0 * u + 38.4 * u * u - 10.666666666667 * u * u * u - 0.066666666667 / (u * u * u));
#ifdef EVALPOTENTIAL
      if(u < 0.5)
        wp = -2.8 + u * u * (5.333333333333 + u * u * (6.4 * u - 9.6));
      else
        wp = -3.2 + 0.066666666667 / u + u * u * (10.666666666667 + u * (-16.0 + u * (9.6 - 2.133333333333 * u)));

      facpot = mass * h_inv * wp;
#endif
    }

  int tabindex = (int) (asmthfac * r);

  if(tabindex < NTAB)
    {
      fac *= shortrange_table[tabindex];

      *acc_x += FLT(dx * fac);
      *acc_y += FLT(dy * fac);
      *acc_z += FLT(dz * fac);
#ifdef EVALPOTENTIAL
      *pot += FLT(facpot * shortrange_table_potential[tabindex]);
#endif
    }
}
#endif

void modgrav_update_node(int no, int father, int * suns)
{
  if(father < 0)
    Nodes[no].mg_bitflag.level = 1;
  else
    Nodes[no].mg_bitflag.level = Nodes[father].mg_bitflag.level + 1;

  Nodes[no].mg_bitflag.num_daughters =  0;
  for(int j = 0; j < 8; j++)
    if(suns[j] >= All.MaxPart && suns[j] < All.MaxPart + Tree_MaxNodes)
      {
        Nodes[no].mg_bitflag.num_daughters++;
      }

  /* now determine the AMR node flag\n
   * the highest top node gets flag 0\n
   * all nodes that have 8 daughters and are below MaxAMRLevel get 0 as well\n
   * the lowest nodes which have 7 siblings or are at MaxAMRLevel get flag 1\n
   * all cells below get flag 2
   *
   */
  if(father < 0)
    {
      Nodes[no].mg_bitflag.amr_node = 0;       /* highest top node */
      if(Nodes[no].mg_bitflag.num_daughters != 8)
        terminate("root node has less than 8 subnodes");
    }
  else if(Nodes[father].mg_bitflag.amr_node == 0 && (Nodes[no].mg_bitflag.num_daughters < 8 || Nodes[no].mg_bitflag.level == All.MaxAMRLevel))
    {
      if(Nodes[father].mg_bitflag.num_daughters != 8)
        terminate("something wrong when choosing AMR cells!");
      Nodes[no].mg_bitflag.amr_node = 1;
    }
  else if(Nodes[father].mg_bitflag.amr_node > 0)
    Nodes[no].mg_bitflag.amr_node = 2;
  else
    Nodes[no].mg_bitflag.amr_node = 0;

}
