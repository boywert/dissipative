/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/subfind/subfind_serial.c
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

#include "../allvars.h"
#include "../proto.h"

#ifdef SUBFIND
#include "subfind.h"
#include "../fof/fof.h"

/* this file processes the local groups in serial mode */



static int *Head, *Next, *Tail, *Len;



int subfind_process_group_serial(int gr, int Offs, int nsubgroups_cat)
{
  int i, j, k, p, count_cand, count, len, len_non_gas, N, nsubs, part_index, subnr, totlen;
  static struct unbind_data *ud;

  while(PS[Offs].GrNr != Group[gr].GrNr)
    {
      Offs++;
      if(Offs >= NumPart)
        {
          char buf[1000];
          sprintf(buf, "don't find a particle for groupnr=%d\n", Group[gr].GrNr);

          for(int i = 0; i < NumPart; i++)
            printf("task=%d i=%d PS[i].GrNr=%d\n", ThisTask, i, PS[i].GrNr);

          terminate(buf);
        }
    }

  N = Group[gr].Len;
  GrNr = Group[gr].GrNr;

#ifdef TRACER_PARTICLE
  // since we are sorted by GrNr and then DM_Density, tracers (with DM_Density=-1) should be at the end
  // modify the count so that they are not processed in the serial subfind
  N -= Group[gr].LenType[TRACER_PARTICLE];

  // sanity checks
  for(i = 0; i < N; i++)
    if(PS[Offs + i].GrNr != Group[gr].GrNr || P[Offs + i].Type == TRACER_PARTICLE)
      terminate("752340");

  for(i = N; i < N + Group[gr].LenType[TRACER_PARTICLE]; i++)
    if(PS[Offs + i].GrNr != Group[gr].GrNr || P[Offs + i].Type != TRACER_PARTICLE)
      terminate("235098");
#endif

  subfind_loctree_treeallocate((int) (All.TreeAllocFactor * N) + NTopnodes, NumPart);

  for(int i = 0; i < N; i++)
    if(PS[Offs + i].GrNr != Group[gr].GrNr)
      terminate("task=%d, gr=%d: don't have the number of particles for GrNr=%d i=%d group-len:N=%d found=%d before=%d\n", ThisTask, gr, Group[gr].GrNr, i, N, PS[Offs + i].GrNr, PS[Offs - 1].GrNr);

  candidates = (struct cand_dat *) mymalloc_movable(&candidates, "candidates", N * sizeof(struct cand_dat));

  Head = (int *) mymalloc_movable(&Head, "Head", N * sizeof(int));
  Next = (int *) mymalloc_movable(&Next, "Next", N * sizeof(int));
  Tail = (int *) mymalloc_movable(&Tail, "Tail", N * sizeof(int));
  Len = (int *) mymalloc_movable(&Len, "Len", N * sizeof(int));
  ud = (struct unbind_data *) mymalloc_movable(&ud, "ud", N * sizeof(struct unbind_data));

  Head -= Offs;
  Next -= Offs;
  Tail -= Offs;
  Len -= Offs;

  for(int i = 0; i < N; i++)
    ud[i].index = Offs + i;

  subfind_loctree_findExtent(N, ud);

  subfind_loctree_treebuild(N, &ud);    /* build tree for all particles of this group */

#ifdef SUBFIND_EXTENDED_PROPERTIES
  // compute the binding energy of FOF group
  double Epot = 0;
  for(int i = 0; i < N; i++)
    {
      int p = ud[i].index;
      double pot = subfind_loctree_treeevaluate_potential(p);

      // note: add self-energy
      pot += P[p].Mass / (All.ForceSoftening[P[p].SofteningType] / 2.8);        // (P[p].Soft / 2.8);

      // multiply with G, scale by scale factor
      pot *= All.G / All.cf_atime;

      Epot += (P[p].Mass / 2) * pot;
    }
  Group[gr].Epot = Epot;
#endif

  for(int i = Offs; i < Offs + N; i++)
    Head[i] = Next[i] = Tail[i] = -1;

  /* note: particles are already ordered in the order of decreasing density */


#ifndef ADD_GROUP_PROPERTIES
  int ss, ngbs, ndiff, head = 0, head_attach;
  int listofdifferent[2], prev;
  int ngb_index, rank;
  int desngb = All.DesLinkNgb;


  for(i = 0, count_cand = 0; i < N; i++)
    {
      part_index = Offs + i;

      MyDouble *pos;
#ifdef CELL_CENTER_GRAVITY
      if(P[part_index].Type == 0)
        pos = PS[part_index].Center;
      else
#endif
        pos = P[part_index].Pos;

      subfind_locngb_treefind(pos, desngb, PS[part_index].Hsml);

      /* note: returned neighbours are already sorted by distance */


      for(k = 0, ndiff = 0, ngbs = 0; k < desngb && ngbs < 2 && ndiff < 2; k++)
        {
          ngb_index = R2list[k].index;

          if(ngb_index != part_index)   /* to exclude the particle itself */
            {
              /* we only look at neighbours that are denser */
              if(PS[ngb_index].Density > PS[part_index].Density)
                {
                  ngbs++;

                  if(Head[ngb_index] >= 0)      /* neighbor is attached to a group */
                    {
                      if(ndiff == 1)
                        if(listofdifferent[0] == Head[ngb_index])
                          continue;

                      /* a new group has been found */
                      listofdifferent[ndiff++] = Head[ngb_index];
                    }
                  else
                    terminate
                      ("this may not occur: ThisTask=%d gr=%d k=%d i=%d part_index=%d ngb_index = %d  head[ngb_index]=%d P[part_index].DM_Density=%g %g GrNrs= %d %d \n",
                       ThisTask, gr, k, i, part_index, ngb_index, Head[ngb_index], PS[part_index].Density, PS[ngb_index].Density, PS[part_index].GrNr, PS[ngb_index].GrNr);
                }
            }
        }

      switch (ndiff)            /* treat the different possible cases */
        {
        case 0:                /* this appears to be a lonely maximum -> new group */
          head = part_index;
          Head[part_index] = Tail[part_index] = part_index;
          Len[part_index] = 1;
          Next[part_index] = -1;
          break;

        case 1:                /* the particle is attached to exactly one group */
          head = listofdifferent[0];
          Head[part_index] = head;
          Next[Tail[head]] = part_index;
          Tail[head] = part_index;
          Len[head]++;
          Next[part_index] = -1;
          break;

        case 2:                /* the particle merges two groups together */
          head = listofdifferent[0];
          head_attach = listofdifferent[1];
          if(Len[head_attach] > Len[head] || (Len[head_attach] == Len[head] && head_attach < head))     /* other group is longer, swap them. for equal length, take the larger head value */
            {
              head = listofdifferent[1];
              head_attach = listofdifferent[0];
            }

          /* only in case the attached group is long enough we bother to register is 
             as a subhalo candidate */

          if(Len[head_attach] >= All.DesLinkNgb)
            {
              candidates[count_cand].len = Len[head_attach];
              candidates[count_cand].head = Head[head_attach];
              count_cand++;
            }

          /* now join the two groups */
          Next[Tail[head]] = head_attach;
          Tail[head] = Tail[head_attach];
          Len[head] += Len[head_attach];

          ss = head_attach;
          do
            {
              Head[ss] = head;
            }
          while((ss = Next[ss]) >= 0);

          /* finally, attach the particle */
          Head[part_index] = head;
          Next[Tail[head]] = part_index;
          Tail[head] = part_index;
          Len[head]++;
          Next[part_index] = -1;
          break;

        default:
          terminate("can't be!");
          break;
        }
    }


  /* add the full thing as a subhalo candidate */
  for(i = 0, prev = -1; i < N; i++)
    {
      if(Head[Offs + i] == Offs + i)
        if(Next[Tail[Offs + i]] == -1)
          {
            if(prev < 0)
              head = Offs + i;
            if(prev >= 0)
              Next[prev] = Offs + i;

            prev = Tail[Offs + i];
          }
    }

  candidates[count_cand].len = N;
  candidates[count_cand].head = head;
  count_cand++;

  /* go through them once and assign the rank */
  for(i = 0, p = head, rank = 0; i < N; i++)
    {
      Len[p] = rank++;
      p = Next[p];
    }

  /* for each candidate, we now pull out the rank of its head */
  for(k = 0; k < count_cand; k++)
    candidates[k].rank = Len[candidates[k].head];

  for(i = Offs; i < Offs + N; i++)
    Tail[i] = -1;

  for(k = 0, nsubs = 0; k < count_cand; k++)
    {
      for(i = 0, p = candidates[k].head, len = 0; i < candidates[k].len; i++, p = Next[p])
        if(Tail[p] < 0)
          ud[len++].index = p;

      if(len >= All.DesLinkNgb)
        len = subfind_unbind(ud, len, &len_non_gas);

      if(len >= All.DesLinkNgb)
        {
          /* ok, we found a substructure */

          for(i = 0; i < len; i++)
            Tail[ud[i].index] = nsubs;  /* we use this to flag the substructures */

          candidates[k].nsub = nsubs;
          candidates[k].bound_length = len;
          nsubs++;
        }
      else
        {
          candidates[k].nsub = -1;
          candidates[k].bound_length = 0;
        }
    }

#ifdef VERBOSE
  printf("\nGroupLen=%d  (gr=%d)\n", N, gr);
  printf("Number of substructures: %d (before unbinding: %d)\n", nsubs, count_cand);
#endif

  mysort(candidates, count_cand, sizeof(struct cand_dat), subfind_compare_serial_candidates_boundlength);

  /* now we determine the parent subhalo for each candidate */
  for(k = 0; k < count_cand; k++)
    {
      candidates[k].subnr = k;
      candidates[k].parent = 0;
    }

  mysort(candidates, count_cand, sizeof(struct cand_dat), subfind_compare_serial_candidates_rank);


  for(k = 0; k < count_cand; k++)
    {
      for(j = k + 1; j < count_cand; j++)
        {
          if(candidates[j].rank > candidates[k].rank + candidates[k].len)
            break;

          if(candidates[k].rank + candidates[k].len >= candidates[j].rank + candidates[j].len)
            {
              if(candidates[k].bound_length >= All.DesLinkNgb)
                candidates[j].parent = candidates[k].subnr;
            }
          else
            {
              char buf[1000];
              sprintf(buf, "k=%d|%d has rank=%d and len=%d.  j=%d has rank=%d and len=%d bound=%d\n",
                      k, count_cand, (int) candidates[k].rank, candidates[k].len, (int) candidates[k].bound_length, candidates[j].rank, (int) candidates[j].len, candidates[j].bound_length);
              terminate(buf);
            }
        }
    }

  mysort(candidates, count_cand, sizeof(struct cand_dat), subfind_compare_serial_candidates_subnr);

#endif


#ifdef ADD_GROUP_PROPERTIES
  int previous_subnr = -1;
  int previous_partindex = -1;
  count_cand = 0;
  nsubs = 0;

  for(int i = 0; i < N; i++)
    {
      part_index = Offs + i;

      if(P[part_index].OriginalSubNr != previous_subnr)
        {
          previous_subnr = P[part_index].OriginalSubNr;
          candidates[count_cand].head = part_index;
          count_cand++;
          nsubs++;

          candidates[count_cand - 1].bound_length = 0;
          candidates[count_cand - 1].len = 0;
          candidates[count_cand - 1].nsub = nsubs;
        }
      else
        {
          Next[previous_partindex] = part_index;
        }

      Tail[part_index] = nsubs;

      candidates[count_cand - 1].bound_length++;
      candidates[count_cand - 1].len++;

      previous_partindex = part_index;
    }

  if(P[candidates[count_cand - 1].head].OriginalSubNr == 2000000000)
    {
      count_cand--;
      nsubs--;
    }

  /* now calculate binding energies for each subhalo */
  for(k = 0, subnr = 0, totlen = 0; k < nsubs; k++)
    {
      len = candidates[k].bound_length;

      for(i = 0, p = candidates[k].head, count = 0; i < candidates[k].len; i++)
        {
          if(Tail[p] == candidates[k].nsub)
            ud[count++].index = p;

          p = Next[p];
        }

      len = subfind_unbind(ud, len, &len_non_gas);
    }

  if(nsubs != GroupCat[GrNr].Nsubs)
    terminate("nsubs=%d != GroupCat[GrNr=%d].Nsubs=%d", nsubs, GrNr, GroupCat[GrNr].Nsubs);

#endif


  /* now determine the properties */
  Group[gr].Nsubs = nsubs;
  Group[gr].Pos[0] = Group[gr].CM[0];
  Group[gr].Pos[1] = Group[gr].CM[1];
  Group[gr].Pos[2] = Group[gr].CM[2];

  for(k = 0, subnr = 0, totlen = 0; k < nsubs; k++)
    {
      len = candidates[k].bound_length;

#ifdef VERBOSE
      printf("subnr=%d  SubLen=%d\n", subnr, len);
#endif

      totlen += len;

      for(i = 0, p = candidates[k].head, count = 0; i < candidates[k].len; i++)
        {
          if(Tail[p] == candidates[k].nsub)
            ud[count++].index = p;

          p = Next[p];
        }

      if(count != len)
        terminate("count=%d != len=%d  k=%d subnr=%d  nsubs=%d", count, len, k, subnr, nsubs);

      if(Nsubgroups > MaxNsubgroups)
        terminate("Nsubgroups = %d >= MaxNsubgroups = %d", Nsubgroups, MaxNsubgroups);

      subfind_determine_sub_halo_properties(ud, len, &SubGroup[Nsubgroups], GrNr, subnr, 0, nsubgroups_cat);

      SubGroup[Nsubgroups].SubParent = candidates[k].parent;
      SubGroup[Nsubgroups].SubNr = subnr;
      SubGroup[Nsubgroups].GrNr = Group[gr].GrNr;

      if(subnr == 0)
        {
          for(j = 0; j < 3; j++)
            Group[gr].Pos[j] = SubGroup[Nsubgroups].Pos[j];
        }

      Nsubgroups++;

      /* Let's now assign the subgroup number */

      for(i = 0; i < len; i++)
        PS[ud[i].index].SubNr = subnr;

      subnr++;
    }

#ifdef VERBOSE
  printf("Fuzz=%d\n", N - totlen);
#endif

  myfree(ud);
  myfree(Len + Offs);
  myfree(Tail + Offs);
  myfree(Next + Offs);
  myfree(Head + Offs);

  myfree(candidates);

  subfind_loctree_treefree();

  return Offs;
}


int subfind_unbind(struct unbind_data *ud, int len, int *len_non_gas)
{
  double *bnd_energy, energy_limit, weakly_bound_limit = 0;
  int i, j, p, minindex, unbound, phaseflag, iter = 0;
  double ddxx, s[3], dx[3], v[3], dv[3], pos[3];
  double vel_to_phys, H_of_a, atime, pot, minpot = 0;
  double boxsize, xtmp;
  double TotMass;

  boxsize = All.BoxSize;

  if(All.ComovingIntegrationOn) // TODO simplyfy by using All.cf_atime
    {
      vel_to_phys = 1.0 / All.Time;
      H_of_a = hubble_function(All.Time);
      atime = All.Time;
    }
  else
    {
      vel_to_phys = atime = 1;
      H_of_a = 0;
    }

  bnd_energy = (double *) mymalloc("bnd_energy", len * sizeof(double));

  phaseflag = 0;                /* this means we will recompute the potential for all particles */

  do
    {
      subfind_loctree_treebuild(len, &ud);

      /* let's compute the potential  */

      if(phaseflag == 0)        /* redo it for all the particles */
        {
          for(i = 0, minindex = -1, minpot = 1.0e30; i < len; i++)
            {
              p = ud[i].index;

              pot = subfind_loctree_treeevaluate_potential(p);

              PS[p].Potential = All.G / All.cf_atime * pot;

              if(PS[p].Potential < minpot || minindex == -1)
                {
                  minpot = PS[p].Potential;
                  minindex = p;
                }
            }

#ifdef CELL_CENTER_GRAVITY
          if(P[minindex].Type == 0)
            {
              for(j = 0; j < 3; j++)
                pos[j] = PS[minindex].Center[j];        /* position of minimum potential */
            }
          else
#endif
            {
              for(j = 0; j < 3; j++)
                pos[j] = P[minindex].Pos[j];    /* position of minimum potential */
            }
        }
      else
        {
          /* we only repeat for those close to the unbinding threshold */
          for(i = 0; i < len; i++)
            {
              p = ud[i].index;

              if(PS[p].BindingEnergy >= weakly_bound_limit)
                {
                  pot = subfind_loctree_treeevaluate_potential(p);

                  PS[p].Potential *= All.G / All.cf_atime;
                }
            }
        }

      /* let's get bulk velocity and the center-of-mass */

      v[0] = v[1] = v[2] = 0;
      s[0] = s[1] = s[2] = 0;

      for(i = 0, TotMass = 0; i < len; i++)
        {
          p = ud[i].index;

          for(j = 0; j < 3; j++)
            {
#ifdef CELL_CENTER_GRAVITY
              if(P[p].Type == 0)
                ddxx = GRAVITY_NEAREST_X(PS[p].Center[j] - pos[j]);
              else
#endif
                ddxx = GRAVITY_NEAREST_X(P[p].Pos[j] - pos[j]);
              s[j] += P[p].Mass * ddxx;
              v[j] += P[p].Mass * P[p].Vel[j];
            }
          TotMass += P[p].Mass;
        }

      for(j = 0; j < 3; j++)
        {
          v[j] /= TotMass;
          s[j] /= TotMass;      /* center-of-mass */

          s[j] += pos[j];

#ifdef PERIODIC
          while(s[j] < 0)
            s[j] += boxsize;
          while(s[j] >= boxsize)
            s[j] -= boxsize;
#endif
        }

      for(i = 0; i < len; i++)
        {
          p = ud[i].index;

          for(j = 0; j < 3; j++)
            {
              dv[j] = vel_to_phys * (P[p].Vel[j] - v[j]);
#ifdef CELL_CENTER_GRAVITY
              if(P[p].Type == 0)
                dx[j] = atime * GRAVITY_NEAREST_X(PS[p].Center[j] - s[j]);
              else
#endif
                dx[j] = atime * GRAVITY_NEAREST_X(P[p].Pos[j] - s[j]);

              dv[j] += H_of_a * dx[j];
            }

          PS[p].BindingEnergy = PS[p].Potential + 0.5 * (dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]);
          PS[p].BindingEnergy += All.G / All.cf_atime * P[p].Mass / (All.ForceSoftening[P[p].SofteningType] / 2.8);  /* note: add self-energy */

          if(P[p].Type == 0)
            PS[p].BindingEnergy += PS[p].Utherm;

          bnd_energy[i] = PS[p].BindingEnergy;
        }

      mysort(bnd_energy, len, sizeof(double), subfind_compare_binding_energy);  /* largest comes first! */

      energy_limit = bnd_energy[(int) (0.25 * len)];

      for(i = 0, unbound = 0; i < len - 1; i++)
        {
          if(bnd_energy[i] > 0)
            unbound++;
          else
            unbound--;

          if(unbound <= 0)
            break;
        }
      weakly_bound_limit = bnd_energy[i];

      /* now omit unbound particles,  but at most 1/4 of the original size */

      for(i = 0, unbound = 0, *len_non_gas = 0; i < len; i++)
        {
          p = ud[i].index;
          if(PS[p].BindingEnergy > 0 && PS[p].BindingEnergy > energy_limit)
            {
              unbound++;
              ud[i] = ud[len - 1];
              i--;
              len--;
            }
          else if(P[p].Type != 0)
            (*len_non_gas)++;
        }

      if(len < All.DesLinkNgb)
        break;

      if(phaseflag == 0)
        {
          if(unbound > 0)
            phaseflag = 1;
        }
      else
        {
          if(unbound == 0)
            {
              phaseflag = 0;    /* this will make us repeat everything once more for all particles */
              unbound = 1;
            }
        }

      if(iter++ > MAXITER)
        terminate("iter > MAXITER = %d", MAXITER);
    }
  while(unbound > 0);

  myfree(bnd_energy);

  return (len);
}


#ifdef SUBFIND_EXTENDED_PROPERTIES
int subfind_fof_calc_am_serial(int gr, int Offs, int snapnr, int ngroups_cat)
{
  long long index;
  int len, i, k;
  double Pos_pbc[3], Vel_tot[3], gr_Jtot[3], gr_Jdm[3], gr_Jgas[3], gr_Jstars[3], jpart[3];
  double gr_CMFrac, gr_CMFracType[NTYPES], gr_Ekin, gr_Ethr;
  int gr_len_dm;
  double gr_mass, gr_mass_gas, gr_mass_stars;
  int ptype;

  while(PS[Offs].GrNr != Group[gr].GrNr)
    {
      Offs++;
      if(Offs >= NumPart)
        {
          char buf[1000];
          sprintf(buf, "don't find a particle for groupnr=%d\n", Group[gr].GrNr);

          for(i = 0; i < NumPart; i++)
            printf("task=%d i=%d PS[i].GrNr=%d\n", ThisTask, i, PS[i].GrNr);

          terminate(buf);
        }
    }

  len = Group[gr].Len;

#ifdef TRACER_PARTICLE
  // since we are sorted by GrNr and then DM_Density, tracers (with DM_Density=-1) should be at the end
  // modify the count so that they are not processed in the serial subfind
  len -= Group[gr].LenType[TRACER_PARTICLE];

  // sanity checks
  for(i = 0; i < N; i++)
    if(PS[Offs + i].GrNr != Group[gr].GrNr || P[Offs + i].Type == TRACER_PARTICLE)
      terminate("752340");

  for(i = N; i < N + Group[gr].LenType[TRACER_PARTICLE]; i++)
    if(PS[Offs + i].GrNr != Group[gr].GrNr || P[Offs + i].Type != TRACER_PARTICLE)
      terminate("235098");
#endif

  struct unbind_data *ud = (struct unbind_data *) mymalloc("ud", len * sizeof(struct unbind_data));

  // get all fof particles
  for(i = 0; i < len; i++)
    ud[i].index = Offs + i;

  // initialize
  gr_CMFrac = 0;
  gr_Ekin = 0;
  gr_Ethr = 0;

  for(k = 0; k < 3; k++)
    {
      gr_Jtot[k] = 0;
      gr_Jdm[k] = 0;
      gr_Jgas[k] = 0;
      gr_Jstars[k] = 0;
    }
  for(k = 0; k < NTYPES; k++)
    {
      gr_CMFracType[k] = 0;
    }


  // calc angular momentum for dm, gas, stars
  for(k = 0; k < len; k++)
    {
      index = ud[k].index;
      ptype = P[index].Type;

#ifdef GFM_WINDS
      /* count wind as gas for mass, but not for LenType since we use this to construct offset tables */
      if(P[index].Type == 4 && STP(index).BirthTime < 0)
        ptype = 0;
#endif

      for(i = 0; i < 3; i++)
        Pos_pbc[i] = P[index].Pos[i] - Group[gr].Pos[i];

      for(i = 0; i < 3; i++)
        Pos_pbc[i] = fof_periodic(Pos_pbc[i]);

      for(i = 0; i < 3; i++)
        Pos_pbc[i] = Pos_pbc[i] * All.cf_atime; // units: phys kpc/h 

      for(i = 0; i < 3; i++)
        Vel_tot[i] = P[index].Vel[i] / All.cf_atime - Group[gr].Vel[i] / All.cf_atime + All.cf_Hrate * Pos_pbc[i];

      gr_Ekin += (P[index].Mass / 2) * (Vel_tot[0] * Vel_tot[0] + Vel_tot[1] * Vel_tot[1] + Vel_tot[2] * Vel_tot[2]);
      if(P[index].Type == 0)
        gr_Ethr += P[index].Mass * SphP[PS[index].OldIndex].Utherm;

      gr_Jtot[0] += P[index].Mass * (Pos_pbc[1] * Vel_tot[2] - Pos_pbc[2] * Vel_tot[1]);
      gr_Jtot[1] += P[index].Mass * (Pos_pbc[2] * Vel_tot[0] - Pos_pbc[0] * Vel_tot[2]);
      gr_Jtot[2] += P[index].Mass * (Pos_pbc[0] * Vel_tot[1] - Pos_pbc[1] * Vel_tot[0]);

      if(ptype == 1)            // dm illustris
        {
          gr_Jdm[0] += P[index].Mass * (Pos_pbc[1] * Vel_tot[2] - Pos_pbc[2] * Vel_tot[1]);
          gr_Jdm[1] += P[index].Mass * (Pos_pbc[2] * Vel_tot[0] - Pos_pbc[0] * Vel_tot[2]);
          gr_Jdm[2] += P[index].Mass * (Pos_pbc[0] * Vel_tot[1] - Pos_pbc[1] * Vel_tot[0]);
        }
      if(ptype == 0)            // gas (incl. winds)
        {
          gr_Jgas[0] += P[index].Mass * (Pos_pbc[1] * Vel_tot[2] - Pos_pbc[2] * Vel_tot[1]);
          gr_Jgas[1] += P[index].Mass * (Pos_pbc[2] * Vel_tot[0] - Pos_pbc[0] * Vel_tot[2]);
          gr_Jgas[2] += P[index].Mass * (Pos_pbc[0] * Vel_tot[1] - Pos_pbc[1] * Vel_tot[0]);
        }
      if(ptype == 4)            // stars
        {
          gr_Jstars[0] += P[index].Mass * (Pos_pbc[1] * Vel_tot[2] - Pos_pbc[2] * Vel_tot[1]);
          gr_Jstars[1] += P[index].Mass * (Pos_pbc[2] * Vel_tot[0] - Pos_pbc[0] * Vel_tot[2]);
          gr_Jstars[2] += P[index].Mass * (Pos_pbc[0] * Vel_tot[1] - Pos_pbc[1] * Vel_tot[0]);
        }
    }

  Group[gr].Ekin = gr_Ekin;
  Group[gr].Ethr = gr_Ethr;
  for(i = 0; i < 3; i++)
    {
      Group[gr].J[i] = gr_Jtot[i];
      Group[gr].JDM[i] = gr_Jdm[i];
      Group[gr].JGas[i] = gr_Jgas[i];
      Group[gr].JStars[i] = gr_Jstars[i];
    }

  // calc counter-rotating fractions
  gr_len_dm = 0;
  gr_mass = gr_mass_gas = gr_mass_stars = 0;

  for(k = 0; k < len; k++)
    {
      index = ud[k].index;
      ptype = P[index].Type;
#ifdef GFM_WINDS
      /* count wind as gas for mass, but not for LenType since we use this to construct offset tables */
      if(P[index].Type == 4 && STP(index).BirthTime < 0)
        ptype = 0;
#endif

      for(i = 0; i < 3; i++)
        Pos_pbc[i] = P[index].Pos[i] - Group[gr].Pos[i];

      for(i = 0; i < 3; i++)
        Pos_pbc[i] = fof_periodic(Pos_pbc[i]);

      for(i = 0; i < 3; i++)
        Pos_pbc[i] = Pos_pbc[i] * All.cf_atime; // units: phys kpc/h 

      for(i = 0; i < 3; i++)
        Vel_tot[i] = P[index].Vel[i] / All.cf_atime - Group[gr].Vel[i] / All.cf_atime + All.cf_Hrate * Pos_pbc[i];

      jpart[0] = P[index].Mass * (Pos_pbc[1] * Vel_tot[2] - Pos_pbc[2] * Vel_tot[1]);
      jpart[1] = P[index].Mass * (Pos_pbc[2] * Vel_tot[0] - Pos_pbc[0] * Vel_tot[2]);
      jpart[2] = P[index].Mass * (Pos_pbc[0] * Vel_tot[1] - Pos_pbc[1] * Vel_tot[0]);

      gr_mass += P[index].Mass;
      if((gr_Jtot[0] * jpart[0] + gr_Jtot[1] * jpart[1] + gr_Jtot[2] * jpart[2]) < 0.)
        gr_CMFrac += P[index].Mass;     // / Group[gr].Mass;

      if(ptype == 1)            // dm illustris
        {
          gr_len_dm++;
          if((gr_Jdm[0] * jpart[0] + gr_Jdm[1] * jpart[1] + gr_Jdm[2] * jpart[2]) < 0.)
            gr_CMFracType[1]++;
        }
      if(ptype == 0)            // gas (incl. winds)
        {
          gr_mass_gas += P[index].Mass;
          if((gr_Jgas[0] * jpart[0] + gr_Jgas[1] * jpart[1] + gr_Jgas[2] * jpart[2]) < 0.)
            gr_CMFracType[0] += P[index].Mass;  // / Group[gr].MassType[0];
        }
      if(ptype == 4)            // stars
        {
          gr_mass_stars += P[index].Mass;
          if((gr_Jstars[0] * jpart[0] + gr_Jstars[1] * jpart[1] + gr_Jstars[2] * jpart[2]) < 0.)
            gr_CMFracType[4] += P[index].Mass;  // / Group[gr].MassType[4];
        }
    }

  gr_CMFrac /= gr_mass;         //Group[gr].Mass;
  gr_CMFracType[1] /= gr_len_dm;
  gr_CMFracType[0] /= gr_mass_gas;      //Group[gr].MassType[0];
  gr_CMFracType[4] /= gr_mass_stars;    //Group[gr].MassType[4];

  Group[gr].CMFrac = gr_CMFrac;
  for(i = 0; i < NTYPES; i++)
    Group[gr].CMFracType[i] = gr_CMFracType[i];

  myfree(ud);
  return Offs;
}
#endif


#endif
