/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/subfind/subfind_properties.c
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

/* if parallel_flag is set, the code calculates the properties for a subhalo distributed onto several processors
 */
void subfind_determine_sub_halo_properties(struct unbind_data *d, int num, struct subgroup_properties *subgroup, int grnr, int subnr, int parallel_flag, int nsubgroups_cat)
{
  int i, j, p, len_type[NTYPES], len_type_loc[NTYPES], totlen;
  double s[3], v[3], pos[3], vel[3], spin[3], cm[3], veldisp, max, vel_to_phys, H_of_a, minpot;
#ifdef MHD
  double bfld_halo, bfld_disk, bfld_vol_halo, bfld_vol_disk;
#endif
#ifdef SUBFIND_EXTENDED_PROPERTIES
  double Ekin = 0, Epot = 0, Ethr = 0, Jdm[3], Jgas[3], Jstars[3], CMFrac, CMFracType[NTYPES];
  double Jdm_inHalfRad[3], Jgas_inHalfRad[3], Jstars_inHalfRad[3], CMFrac_inHalfRad, CMFracType_inHalfRad[NTYPES];
  double Jdm_inRad[3], Jgas_inRad[3], Jstars_inRad[3], CMFrac_inRad, CMFracType_inRad[NTYPES];
  double jpart[3], Jtot[3], Jtot_inRad[3], Jtot_inHalfRad[3];
  double sinrad[3], sinhalfrad[3], vinrad[3], vinhalfrad[3];
#endif
  double lx, ly, lz, dv[3], dx[3], disp, rr_tmp, disp_tmp, halfmassrad = 0, halfmassradtype[NTYPES];
  double boxsize, ddxx, vmax, vmaxrad, maxrad;
  double mass, massinrad, massinhalfrad, massinmaxrad;
  double mass_tab[NTYPES], massinrad_tab[NTYPES], massinhalfrad_tab[NTYPES], massinmaxrad_tab[NTYPES];
  double xtmp;
#ifdef GFM_STELLAR_PHOTOMETRICS
  double ytmp, ztmp;
#endif

  sort_r2list *rr_list = 0;
  int minindex;
  MyIDType mostboundid;

#ifdef USE_SFR
  double sfr = 0, sfrinrad = 0, sfrinhalfrad = 0, sfrinmaxrad = 0, gasMassSfr = 0;
#endif
#ifdef GFM_STELLAR_EVOLUTION
  double gasMassMetallicity = 0;
  double gasMassMetals[GFM_N_CHEM_ELEMENTS];
  for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
    gasMassMetals[j] = 0;

  double stellarMassMetallicity = 0;
  double stellarMassMetals[GFM_N_CHEM_ELEMENTS];
  for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
    stellarMassMetals[j] = 0;

  double gasMassMetallicityHalfRad = 0;
  double gasMassMetalsHalfRad[GFM_N_CHEM_ELEMENTS];
  for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
    gasMassMetalsHalfRad[j] = 0;

  double stellarMassMetallicityHalfRad = 0;
  double stellarMassMetalsHalfRad[GFM_N_CHEM_ELEMENTS];
  for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
    stellarMassMetalsHalfRad[j] = 0;

  double gasMassMetallicityMaxRad = 0;
  double gasMassMetalsMaxRad[GFM_N_CHEM_ELEMENTS];
  for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
    gasMassMetalsMaxRad[j] = 0;

  double stellarMassMetallicityMaxRad = 0;
  double stellarMassMetalsMaxRad[GFM_N_CHEM_ELEMENTS];
  for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
    stellarMassMetalsMaxRad[j] = 0;

#ifdef GFM_DUST
  int chan = 0;
  double gasMassDustMetallicity = 0;
  double gasMassDustMetallicityHalfRad = 0;
  double gasMassDustMetallicityMaxRad = 0;
#endif

#ifdef USE_SFR
  double gasMassMetallicitySfr = 0;
  double gasMassMetalsSfr[GFM_N_CHEM_ELEMENTS];
  for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
    gasMassMetalsSfr[j] = 0;

  double gasMassMetallicitySfrWeighted = 0;
  double gasMassMetalsSfrWeighted[GFM_N_CHEM_ELEMENTS];
  for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
    gasMassMetalsSfrWeighted[j] = 0;

#ifdef GFM_DUST
  double gasMassDustMetallicitySfr = 0;
  double gasMassDustMetallicitySfrWeighted = 0;
#endif

#endif
#endif
#ifdef GFM_STELLAR_PHOTOMETRICS
  double max_stellar_rad = 0;
  double brightness_limit_rad = 0;
  double stellar_mass_in_phot_rad = 0;
#endif
#ifdef BLACK_HOLES
  double bh_Mdot = 0, bh_Mass = 0;
#endif
#ifdef GFM_WINDS
  double windMass = 0;
#endif

  boxsize = All.BoxSize;

  vel_to_phys = 1.0 / All.cf_atime;

  if(All.ComovingIntegrationOn) // TODO?
    H_of_a = hubble_function(All.Time);
  else
    H_of_a = 0;

  mass = massinrad = massinhalfrad = massinmaxrad = 0;
  for(j = 0; j < NTYPES; j++)
    {
      len_type[j] = 0;
      mass_tab[j] = halfmassradtype[j] = massinrad_tab[j] = massinhalfrad_tab[j] = massinmaxrad_tab[j] = 0;
    }

  for(i = 0, minindex = -1, minpot = 1.0e30; i < num; i++)
    {
      p = d[i].index;
      if(PS[p].Potential < minpot || minindex == -1)
        {
          minpot = PS[p].Potential;
          minindex = p;
        }

#ifdef TRACER_MC
      len_type[TRACER_MC] += get_number_of_tracers(p);
#endif

#ifdef GFM_WINDS_SAVE_PARTTYPE
      /* count wind as new particle type for LenType since we save these separately */
      if((P[p].Type == 4) && (STP(p).BirthTime < 0))
        len_type[GFM_WINDS_SAVE_PARTTYPE]++;
      else
#endif
      len_type[P[p].Type]++;

#ifdef USE_SFR
      if(P[p].Type == 0)
        sfr += SphP[PS[p].OldIndex].Sfr;        /* note: the SphP[] array has not been reordered */
#endif
#ifdef BLACK_HOLES
      if(P[p].Type == 5)
        {
          bh_Mass += BPP(p).BH_Mass;
          bh_Mdot += BPP(p).BH_Mdot;
        }
#endif
#ifdef GFM_WINDS
      if(P[p].Type == 4 && STP(p).BirthTime < 0)
        windMass += P[p].Mass;
#endif
    }

  for(j = 0; j < NTYPES; j++)
    len_type_loc[j] = len_type[j];


  if(parallel_flag)
    {
      int len_typetot[NTYPES];
      MPI_Allreduce(len_type, len_typetot, NTYPES, MPI_INT, MPI_SUM, SubComm);
      for(j = 0; j < NTYPES; j++)
        len_type[j] = len_typetot[j];

      double *minpotlist = mymalloc("minpotlist", SubNTask * sizeof(double));
      MPI_Allgather(&minpot, 1, MPI_DOUBLE, minpotlist, 1, MPI_DOUBLE, SubComm);
      int mincpu;

      for(i = 0, mincpu = -1, minpot = 1.0e30; i < SubNTask; i++)
        if(minpotlist[i] < minpot)
          {
            mincpu = i;
            minpot = minpotlist[mincpu];
          }

      myfree(minpotlist);

      if(mincpu < 0)
        terminate("mincpu < 0");

      if(SubThisTask == mincpu)
        for(j = 0; j < 3; j++)
          {
#ifdef CELL_CENTER_GRAVITY
            if(P[minindex].Type == 0)
              pos[j] = SphP[PS[minindex].OldIndex].Center[j];
            else
#endif
              pos[j] = P[minindex].Pos[j];
          }

      MPI_Bcast(pos, 3, MPI_DOUBLE, mincpu, SubComm);

#ifdef USE_SFR
      double sfrtot;
      MPI_Allreduce(&sfr, &sfrtot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      sfr = sfrtot;
#endif
#ifdef BLACK_HOLES
      double bh_Masstot, bh_Mdottot;
      MPI_Allreduce(&bh_Mass, &bh_Masstot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(&bh_Mdot, &bh_Mdottot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      bh_Mass = bh_Masstot;
      bh_Mdot = bh_Mdottot;
#endif
#ifdef GFM_WINDS
      double windMasstot;
      MPI_Allreduce(&windMass, &windMasstot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      windMass = windMasstot;
#endif
    }
  else
    {
      if(minindex == -1)
        terminate("minindex == -1");

      for(j = 0; j < 3; j++)
        {
#ifdef CELL_CENTER_GRAVITY
          if(P[minindex].Type == 0)
            pos[j] = SphP[PS[minindex].OldIndex].Center[j];
          else
#endif
            pos[j] = P[minindex].Pos[j];
        }
    }


#ifdef ADD_GROUP_PROPERTIES
  /* override minimum potential position with previous position in subfind catalogue */
  /* printf("TASK=%d:  grnr=%d subnr=%d  (%g|%g|%g)  (%g|%g|%g)\n", ThisTask, grnr, subnr,
         pos[0], pos[1], pos[2],
         SubGroupPos[3 * (GroupCat[grnr].FirstSub + subnr) + 0], SubGroupPos[3 * (GroupCat[grnr].FirstSub + subnr) + 1], SubGroupPos[3 * (GroupCat[grnr].FirstSub + subnr) + 2]);
  */
  for(j = 0; j < 3; j++)
    pos[j] = SubGroupPos[3 * (GroupCat[grnr].FirstSub + subnr) + j];
#endif



  /* pos[] now holds the position of minimum potential */
  /* we'll take it that as the center */


  /* determine the particle ID with the smallest binding energy */
  for(i = 0, minindex = -1, minpot = 1.0e30; i < num; i++)
    {
      p = d[i].index;
      if(PS[p].BindingEnergy < minpot || minindex == -1)
        {
          minpot = PS[p].BindingEnergy;
          minindex = p;
        }
    }

  if(parallel_flag)
    {
      double *minpotlist = mymalloc("minpotlist", SubNTask * sizeof(double));
      MPI_Allgather(&minpot, 1, MPI_DOUBLE, minpotlist, 1, MPI_DOUBLE, SubComm);
      int mincpu;

      for(i = 0, mincpu = -1, minpot = 1.0e30; i < SubNTask; i++)
        if(minpotlist[i] < minpot)
          {
            mincpu = i;
            minpot = minpotlist[mincpu];
          }

      myfree(minpotlist);

      if(mincpu < 0)
        terminate("mincpu < 0");

      if(SubThisTask == mincpu)
        {
          mostboundid = P[minindex].ID;
#ifdef FOF_FUZZ_SORT_BY_NEAREST_GROUP
          if(subnr == 0)
            PS[minindex].GroupNr = grnr + 1;
#endif
        }

      MPI_Bcast(&mostboundid, sizeof(mostboundid), MPI_BYTE, mincpu, SubComm);
    }
  else
    {
      if(minindex == -1)
        terminate("minindex == -1");

      mostboundid = P[minindex].ID;
#ifdef FOF_FUZZ_SORT_BY_NEAREST_GROUP
      if(subnr == 0)
        PS[minindex].GroupNr = grnr + 1;
#endif
    }

  /* let's get bulk velocity and the center-of-mass */
  /* here we still take all particles */

  for(j = 0; j < 3; j++)
    s[j] = v[j] = 0;

  for(i = 0; i < num; i++)
    {
      p = d[i].index;
      for(j = 0; j < 3; j++)
        {
          ddxx = GRAVITY_NEAREST_X(P[p].Pos[j] - pos[j]);
          s[j] += P[p].Mass * ddxx;
          v[j] += P[p].Mass * P[p].Vel[j];
        }
      mass += P[p].Mass;

      int ptype = P[p].Type;
#ifdef GFM_WINDS
      /* count wind as gas for mass, but not for LenType since we use this to construct offset tables */
      if(P[p].Type == 4 && STP(p).BirthTime < 0)
        ptype = 0;
#endif
      mass_tab[ptype] += P[p].Mass;
    }

  if(parallel_flag)
    {
      double stot[3], vtot[3], masstot, mass_tabtot[NTYPES];

      MPI_Allreduce(s, stot, 3, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(v, vtot, 3, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(mass_tab, mass_tabtot, NTYPES, MPI_DOUBLE, MPI_SUM, SubComm);

      mass = masstot;
      for(j = 0; j < 3; j++)
        {
          s[j] = stot[j];
          v[j] = vtot[j];
        }

      for(j = 0; j < NTYPES; j++)
        mass_tab[j] = mass_tabtot[j];
    }

  for(j = 0; j < 3; j++)
    {
      s[j] /= mass;             /* center of mass */
      v[j] /= mass;
      vel[j] = vel_to_phys * v[j];
    }

  for(j = 0; j < 3; j++)
    {
      s[j] += pos[j];

#ifdef PERIODIC
      while(s[j] < 0)
        s[j] += boxsize;
      while(s[j] >= boxsize)
        s[j] -= boxsize;
#endif
      cm[j] = s[j];             // this is in comoving coordinates
    }

  disp = lx = ly = lz = 0;
#ifdef SUBFIND_EXTENDED_PROPERTIES
  Jtot[0] = Jtot[1] = Jtot[2] = 0;
  Jdm[0] = Jdm[1] = Jdm[2] = 0;
  Jgas[0] = Jgas[1] = Jgas[2] = 0;
  Jstars[0] = Jstars[1] = Jstars[2] = 0;
#endif

  rr_list = mymalloc("rr_list", sizeof(sort_r2list) * (num + 1));

  for(i = 0; i < num; i++)
    {
      p = d[i].index;

      for(j = 0, rr_tmp = 0, disp_tmp = 0; j < 3; j++)
        {
          ddxx = GRAVITY_NEAREST_X(P[p].Pos[j] - s[j]);
          dx[j] = All.cf_atime * ddxx;
          dv[j] = vel_to_phys * (P[p].Vel[j] - v[j]);
          dv[j] += H_of_a * dx[j];

          disp_tmp += P[p].Mass * dv[j] * dv[j];
          /* for rotation curve computation, take minimum of potential as center */
          ddxx = GRAVITY_NEAREST_X(P[p].Pos[j] - pos[j]);
          ddxx = All.cf_atime * ddxx;
          rr_tmp += ddxx * ddxx;
        }

      lx += P[p].Mass * (dx[1] * dv[2] - dx[2] * dv[1]);
      ly += P[p].Mass * (dx[2] * dv[0] - dx[0] * dv[2]);
      lz += P[p].Mass * (dx[0] * dv[1] - dx[1] * dv[0]);


#ifdef SUBFIND_EXTENDED_PROPERTIES
      for(j = 0; j < 3; j++)    // hubble drifts in velocity now with respect to pot min which we consider as the centre of rotation
        {
          ddxx = GRAVITY_NEAREST_X(P[p].Pos[j] - pos[j]);
          dx[j] = All.cf_atime * ddxx;
          dv[j] = vel_to_phys * (P[p].Vel[j] - v[j]);
          dv[j] += H_of_a * dx[j];
        }

      int ptype = P[p].Type;
#ifdef GFM_WINDS
      /* count wind as gas for mass, but not for LenType since we use this to construct offset tables */
      if(P[p].Type == 4 && STP(p).BirthTime < 0)
        ptype = 0;
#endif

      Ekin += (P[p].Mass / 2) * (dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]);
      Epot += (P[p].Mass / 2) * PS[p].Potential;
      if(P[p].Type == 0)
        Ethr += P[p].Mass * SphP[PS[p].OldIndex].Utherm;

      Jtot[0] += P[p].Mass * (dx[1] * dv[2] - dx[2] * dv[1]);
      Jtot[1] += P[p].Mass * (dx[2] * dv[0] - dx[0] * dv[2]);
      Jtot[2] += P[p].Mass * (dx[0] * dv[1] - dx[1] * dv[0]);

      if(ptype == 1)            // dm illustris
        {
          Jdm[0] += P[p].Mass * (dx[1] * dv[2] - dx[2] * dv[1]);
          Jdm[1] += P[p].Mass * (dx[2] * dv[0] - dx[0] * dv[2]);
          Jdm[2] += P[p].Mass * (dx[0] * dv[1] - dx[1] * dv[0]);
        }
      if(ptype == 0)            // gas (incl. winds!)
        {
          Jgas[0] += P[p].Mass * (dx[1] * dv[2] - dx[2] * dv[1]);
          Jgas[1] += P[p].Mass * (dx[2] * dv[0] - dx[0] * dv[2]);
          Jgas[2] += P[p].Mass * (dx[0] * dv[1] - dx[1] * dv[0]);
        }
      if(ptype == 4)            // stars (previously: StarP[P[p].AuxDataID].BirthTime)
        {
          Jstars[0] += P[p].Mass * (dx[1] * dv[2] - dx[2] * dv[1]);
          Jstars[1] += P[p].Mass * (dx[2] * dv[0] - dx[0] * dv[2]);
          Jstars[2] += P[p].Mass * (dx[0] * dv[1] - dx[1] * dv[0]);
        }
#endif

      rr_tmp = sqrt(rr_tmp);

      rr_list[i].mass = P[p].Mass;
      rr_list[i].r = rr_tmp;
      disp += disp_tmp;
    }

  if(parallel_flag)
    {
      double spintot[3], disptot;
      spin[0] = lx;
      spin[1] = ly;
      spin[2] = lz;
      MPI_Allreduce(spin, spintot, 3, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(&disp, &disptot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      disp = disptot;
      lx = spintot[0];
      ly = spintot[1];
      lz = spintot[2];
#ifdef SUBFIND_EXTENDED_PROPERTIES
      MPI_Allreduce(MPI_IN_PLACE, &Ekin, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, &Epot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, &Ethr, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, Jtot, 3, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, Jdm, 3, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, Jgas, 3, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, Jstars, 3, MPI_DOUBLE, MPI_SUM, SubComm);
#endif
    }

  spin[0] = lx / mass;
  spin[1] = ly / mass;
  spin[2] = lz / mass;

  veldisp = sqrt(disp / (3 * mass));    /* convert to 1d velocity dispersion */

#ifdef SUBFIND_EXTENDED_PROPERTIES
  // counter rotating mass fractions
  CMFrac = 0;
  for(i = 0; i < NTYPES; i++)
    CMFracType[i] = 0;

  for(i = 0; i < num; i++)
    {
      /* identify particle type */
      p = d[i].index;

      /* calculate particle radius */
      for(j = 0; j < 3; j++)
        {
          ddxx = GRAVITY_NEAREST_X(P[p].Pos[j] - pos[j]);       // counter-rotating mass calc with respect to pot min
          dx[j] = All.cf_atime * ddxx;
          dv[j] = vel_to_phys * (P[p].Vel[j] - v[j]);
          dv[j] += H_of_a * dx[j];
        }

      int ptype = P[p].Type;
#ifdef GFM_WINDS
      /* count wind as gas for mass, but not for LenType since we use this to construct offset tables */
      if(P[p].Type == 4 && STP(p).BirthTime < 0)
        ptype = 0;
#endif

      jpart[0] = P[p].Mass * (dx[1] * dv[2] - dx[2] * dv[1]);
      jpart[1] = P[p].Mass * (dx[2] * dv[0] - dx[0] * dv[2]);
      jpart[2] = P[p].Mass * (dx[0] * dv[1] - dx[1] * dv[0]);

      if((Jtot[0] * jpart[0] + Jtot[1] * jpart[1] + Jtot[2] * jpart[2]) < 0.)
        CMFrac += P[p].Mass / mass;

      if(ptype == 1)            // dm illustris
        if((Jdm[0] * jpart[0] + Jdm[1] * jpart[1] + Jdm[2] * jpart[2]) < 0.)
          CMFracType[1] += P[p].Mass / mass_tab[1];
      if(ptype == 0)            // gas (incl. winds!)
        if((Jgas[0] * jpart[0] + Jgas[1] * jpart[1] + Jgas[2] * jpart[2]) < 0.)
          CMFracType[0] += P[p].Mass / mass_tab[0];
      if(ptype == 4)            // stars
        if((Jstars[0] * jpart[0] + Jstars[1] * jpart[1] + Jstars[2] * jpart[2]) < 0.)
          CMFracType[4] += P[p].Mass / mass_tab[4];
    }

  if(parallel_flag)
    {
      MPI_Allreduce(MPI_IN_PLACE, &CMFrac, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, CMFracType, NTYPES, MPI_DOUBLE, MPI_SUM, SubComm);
    }

#endif



  if(parallel_flag)
    parallel_sort_comm(rr_list, num, sizeof(sort_r2list), subfind_compare_dist_rotcurve, SubComm);
  else
    mysort(rr_list, num, sizeof(sort_r2list), subfind_compare_dist_rotcurve);

  /* calculate cumulative mass */
  for(i = 1; i < num; i++)
    rr_list[i].mass += rr_list[i - 1].mass;


  if(parallel_flag)
    {
      double mass_part = 0;
      if(num)
        mass_part = rr_list[num - 1].mass;
      double *masslist = mymalloc("masslist", SubNTask * sizeof(double));
      MPI_Allgather(&mass_part, 1, MPI_DOUBLE, masslist, 1, MPI_DOUBLE, SubComm);

      double massbefore = 0;
      for(i = 0; i < SubThisTask; i++)
        massbefore += masslist[i];

      for(i = 0; i < num; i++)
        rr_list[i].mass += massbefore;

      myfree(masslist);

      /* now calculate rotation curve maximum and half mass radius */

      double halfmassrad_loc = 0;
      sort_r2list *rr_lowlist = mymalloc("rr_lowlist", SubNTask * sizeof(sort_r2list));
      sort_r2list low_element;
      if(num > 0)
        low_element = rr_list[0];
      else
        {
          low_element.mass = 0;
          low_element.r = 0;
        }
      MPI_Allgather(&low_element, sizeof(sort_r2list), MPI_BYTE, rr_lowlist, sizeof(sort_r2list), MPI_BYTE, SubComm);

      rr_list[num].mass = 0;
      rr_list[num].r = 0;

      for(j = SubThisTask + 1; j < SubNTask; j++)
        if(rr_lowlist[j].mass > 0)
          {
            rr_list[num] = rr_lowlist[j];
            break;
          }

      myfree(rr_lowlist);

      int *numlist = mymalloc("numlist", SubNTask * sizeof(int));
      MPI_Allgather(&num, 1, MPI_INT, numlist, 1, MPI_INT, SubComm);

      int nbefore = 0;
      for(i = 0; i < SubThisTask; i++)
        nbefore += numlist[i];

      for(i = num - 1, max = 0, maxrad = 0; i >= 0; i--)
        {
          if((i + nbefore) > 5 && rr_list[i].mass > max * rr_list[i].r)
            {
              max = rr_list[i].mass / rr_list[i].r;
              maxrad = rr_list[i].r;
            }

          if(rr_list[i].mass < 0.5 * mass && rr_list[i + 1].mass >= 0.5 * mass)
            halfmassrad_loc = 0.5 * (rr_list[i].r + rr_list[i + 1].r);
        }

      myfree(numlist);

      MPI_Allreduce(&halfmassrad_loc, &halfmassrad, 1, MPI_DOUBLE, MPI_MAX, SubComm);
      double *maxlist = mymalloc("maxlist", SubNTask * sizeof(double));
      double *maxradlist = mymalloc("maxradlist", SubNTask * sizeof(double));
      MPI_Allgather(&max, 1, MPI_DOUBLE, maxlist, 1, MPI_DOUBLE, SubComm);
      MPI_Allgather(&maxrad, 1, MPI_DOUBLE, maxradlist, 1, MPI_DOUBLE, SubComm);
      for(i = 0, max = maxrad = 0; i < SubNTask; i++)
        {
          if(maxlist[i] > max)
            {
              max = maxlist[i];
              maxrad = maxradlist[i];
            }
        }
      myfree(maxradlist);
      myfree(maxlist);
    }
  else
    {
      for(i = num - 1, max = 0, maxrad = 0; i >= 0; i--)
        {
          if(i > 5 && rr_list[i].mass > max * rr_list[i].r)
            {
              max = rr_list[i].mass / rr_list[i].r;
              maxrad = rr_list[i].r;
            }

          if(i < num - 1)
            if(rr_list[i].mass < 0.5 * mass && rr_list[i + 1].mass >= 0.5 * mass)
              halfmassrad = 0.5 * (rr_list[i].r + rr_list[i + 1].r);
        }
    }

  halfmassrad /= All.cf_atime;
  vmax = sqrt(All.G * max);
  vmaxrad = maxrad / All.cf_atime;

  myfree(rr_list);


  /* half mass radii for different types */
  /* need to recalculate len_type_loc first, because of special particle treatment in GFM */
  for(j = 0; j < NTYPES; j++)
    len_type_loc[j] = 0;

  for(i = 0; i < num; i++)
    {
      p = d[i].index;
      int ptype = P[p].Type;
#ifdef GFM_WINDS
      /* count wind as gas for mass, but not for LenType since we use this to construct offset tables */
      if(P[p].Type == 4 && STP(p).BirthTime < 0)
        ptype = 0;
#endif
      len_type_loc[ptype]++;
    }

  int itmp, type;
  for(type = 0; type < NTYPES; type++)
    {
      rr_list = mymalloc("rr_list", sizeof(sort_r2list) * (len_type_loc[type] + 1));
      itmp = 0;
      for(i = 0; i < num; i++)
        {
          p = d[i].index;

          int ptype = P[p].Type;
#ifdef GFM_WINDS
          /* count wind as gas for mass, but not for LenType since we use this to construct offset tables */
          if(P[p].Type == 4 && STP(p).BirthTime < 0)
            ptype = 0;
#endif

          if(ptype == type)
            {
              for(j = 0, rr_tmp = 0; j < 3; j++)
                {
                  ddxx = GRAVITY_NEAREST_X(P[p].Pos[j] - pos[j]);
                  rr_tmp += ddxx * ddxx;
                }

              rr_tmp = sqrt(rr_tmp);

              rr_list[itmp].mass = P[p].Mass;
              rr_list[itmp].r = rr_tmp;
              itmp++;
            }
        }

      if(itmp != len_type_loc[type])
        terminate("should not occur: %d %d", itmp, len_type_loc[type]);

      if(parallel_flag)
        parallel_sort_comm(rr_list, len_type_loc[type], sizeof(sort_r2list), subfind_compare_dist_rotcurve, SubComm);
      else
        mysort(rr_list, len_type_loc[type], sizeof(sort_r2list), subfind_compare_dist_rotcurve);

      /* calculate cumulative mass */
      for(i = 1; i < len_type_loc[type]; i++)
        rr_list[i].mass = rr_list[i - 1].mass + rr_list[i].mass;

      if(parallel_flag)
        {
          double mass_part = 0;
          if(len_type_loc[type])
            mass_part = rr_list[len_type_loc[type] - 1].mass;
          double *masslist = mymalloc("masslist", SubNTask * sizeof(double));
          MPI_Allgather(&mass_part, 1, MPI_DOUBLE, masslist, 1, MPI_DOUBLE, SubComm);

          double massbefore = 0;
          for(i = 0; i < SubThisTask; i++)
            massbefore += masslist[i];

          for(i = 0; i < len_type_loc[type]; i++)
            rr_list[i].mass += massbefore;

          myfree(masslist);
        }

      /* now calculate half mass radii */
      if(parallel_flag)
        {
          double halfmassrad_loc = 0;
          sort_r2list *rr_lowlist = mymalloc("rr_lowlist", SubNTask * sizeof(sort_r2list));
          sort_r2list low_element;
          if(len_type_loc[type] > 0)
            low_element = rr_list[0];
          else
            {
              low_element.mass = 0;
              low_element.r = 0;
            }

          MPI_Allgather(&low_element, sizeof(sort_r2list), MPI_BYTE, rr_lowlist, sizeof(sort_r2list), MPI_BYTE, SubComm);

          rr_list[len_type_loc[type]].mass = 0;
          rr_list[len_type_loc[type]].r = 0;
          for(j = SubThisTask + 1; j < SubNTask; j++)
            if(rr_lowlist[j].mass > 0)
              {
                rr_list[len_type_loc[type]] = rr_lowlist[j];
                break;
              }

          myfree(rr_lowlist);

          for(i = len_type_loc[type] - 1; i >= 0; i--)
            {
              if(rr_list[i].mass < 0.5 * mass_tab[type] && rr_list[i + 1].mass >= 0.5 * mass_tab[type])
                halfmassrad_loc = 0.5 * (rr_list[i].r + rr_list[i + 1].r);
            }

          MPI_Allreduce(&halfmassrad_loc, &halfmassradtype[type], 1, MPI_DOUBLE, MPI_MAX, SubComm);
        }
      else
        {
          for(i = len_type_loc[type] - 1; i >= 0; i--)
            {
              if(i < len_type_loc[type] - 1)
                if(rr_list[i].mass < 0.5 * mass_tab[type] && rr_list[i + 1].mass >= 0.5 * mass_tab[type])
                  halfmassradtype[type] = 0.5 * (rr_list[i].r + rr_list[i + 1].r);
            }
        }

#ifdef GFM_STELLAR_PHOTOMETRICS
      /* calculate maximum stellar radius */
      if(type == 4)
        {
          if(len_type_loc[type])
            max_stellar_rad = rr_list[len_type_loc[type] - 1].r;
          if(parallel_flag)
            {
              double max_stellar_rad_all;
              MPI_Allreduce(&max_stellar_rad, &max_stellar_rad_all, 1, MPI_DOUBLE, MPI_MAX, SubComm);
              max_stellar_rad = max_stellar_rad_all;
            }
          max_stellar_rad *= 1.01;      /* protect against round-off errors resulting in rxy > max_stellar_rad later */
        }
#endif

      myfree(rr_list);
    }


  /* properties of 'central galaxies', defined in several ways as particles within some radius:
     either (stellar half mass radius) or SUBFIND_GAL_RADIUS_FAC*(stellar half mass radius) or (radius of Vmax) */
#ifdef SUBFIND_EXTENDED_PROPERTIES
  // centre of mass /velocity of particles in half/ stellar mass rad
  sinrad[0] = sinrad[1] = sinrad[2] = 0;
  sinhalfrad[0] = sinhalfrad[1] = sinhalfrad[2] = 0;
  vinrad[0] = vinrad[1] = vinrad[2] = 0;
  vinhalfrad[0] = vinhalfrad[1] = vinhalfrad[2] = 0;
#endif

  for(i = 0; i < num; i++)
    {
      /* identify particle type */
      p = d[i].index;
      int ptype = P[p].Type;
#ifdef GFM_WINDS
      /* count wind as gas for mass, but not for LenType since we use this to construct offset tables */
      if(P[p].Type == 4 && STP(p).BirthTime < 0)
        ptype = 0;
#endif

      /* calculate particle radius */
      for(j = 0, rr_tmp = 0; j < 3; j++)
        {
          ddxx = GRAVITY_NEAREST_X(P[p].Pos[j] - pos[j]);
          rr_tmp += ddxx * ddxx;
        }
      rr_tmp = sqrt(rr_tmp);

      /* properties inside SUBFIND_GAL_RADIUS_FAC*(stellar half mass radius) */
      if(rr_tmp < SUBFIND_GAL_RADIUS_FAC * halfmassradtype[4])
        {
          massinrad += P[p].Mass;
          massinrad_tab[ptype] += P[p].Mass;

#ifdef SUBFIND_EXTENDED_PROPERTIES
          for(j = 0; j < 3; j++)
            {
              ddxx = GRAVITY_NEAREST_X(P[p].Pos[j] - pos[j]);   // comoving (as it should be.)
              sinrad[j] += P[p].Mass * ddxx;
              vinrad[j] += P[p].Mass * P[p].Vel[j];
            }
#endif

          if(ptype == 0)
            {
              if(P[p].Type == 0)
                {
#ifdef USE_SFR
                  sfrinrad += SphP[PS[p].OldIndex].Sfr; /* note: the SphP[] array has not been reordered */
#endif
#ifdef GFM_STELLAR_EVOLUTION
                  gasMassMetallicity += SphP[PS[p].OldIndex].Metallicity * P[p].Mass;
                  for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
                    gasMassMetals[j] += SphP[PS[p].OldIndex].MetalsFraction[j] * P[p].Mass;

#ifdef GFM_DUST
                  for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
                    for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
                      gasMassDustMetallicity += SphP[PS[p].OldIndex].MetalsDustFraction[chan][j] * P[p].Mass;
#endif
#endif
                }
              if(P[p].Type == 4)
                {
#ifdef GFM_STELLAR_EVOLUTION
                  gasMassMetallicity += STP(p).Metallicity * P[p].Mass;
                  for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
                    gasMassMetals[j] += STP(p).MassMetals[j];

#ifdef GFM_DUST
                  for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
                    for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
                      gasMassDustMetallicity += STP(p).MassMetals[j] * STP(p).InitialDustFractions[chan][j];
#endif
#endif
                }
            }
          if(ptype == 4)
            {
#ifdef GFM_STELLAR_EVOLUTION
              stellarMassMetallicity += STP(p).Metallicity * P[p].Mass;
              for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
                stellarMassMetals[j] += STP(p).MassMetals[j];
#endif
            }
        }

      /* properties inside (stellar half mass radius) */
      if(rr_tmp < 1.0 * halfmassradtype[4])
        {
          massinhalfrad += P[p].Mass;
          massinhalfrad_tab[ptype] += P[p].Mass;

#ifdef SUBFIND_EXTENDED_PROPERTIES
          for(j = 0; j < 3; j++)
            {
              ddxx = GRAVITY_NEAREST_X(P[p].Pos[j] - pos[j]);   // comoving (as it should be.)
              sinhalfrad[j] += P[p].Mass * ddxx;
              vinhalfrad[j] += P[p].Mass * P[p].Vel[j];
            }
#endif

          if(ptype == 0)
            {
              if(P[p].Type == 0)
                {
#ifdef USE_SFR
                  sfrinhalfrad += SphP[PS[p].OldIndex].Sfr;     /* note: the SphP[] array has not been reordered */
#endif
#ifdef GFM_STELLAR_EVOLUTION
                  gasMassMetallicityHalfRad += SphP[PS[p].OldIndex].Metallicity * P[p].Mass;
                  for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
                    gasMassMetalsHalfRad[j] += SphP[PS[p].OldIndex].MetalsFraction[j] * P[p].Mass;

#ifdef GFM_DUST
                  for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
                    for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
                      gasMassDustMetallicityHalfRad += SphP[PS[p].OldIndex].MetalsDustFraction[chan][j] * P[p].Mass;
#endif
#endif
                }
              if(P[p].Type == 4)
                {
#ifdef GFM_STELLAR_EVOLUTION
                  gasMassMetallicityHalfRad += STP(p).Metallicity * P[p].Mass;
                  for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
                    gasMassMetalsHalfRad[j] += STP(p).MassMetals[j];

#ifdef GFM_DUST
                  for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
                    for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
                      gasMassDustMetallicityHalfRad += STP(p).MassMetals[j] * STP(p).InitialDustFractions[chan][j];
#endif
#endif
                }
            }
          if(ptype == 4)
            {
#ifdef GFM_STELLAR_EVOLUTION
              stellarMassMetallicityHalfRad += STP(p).Metallicity * P[p].Mass;
              for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
                stellarMassMetalsHalfRad[j] += STP(p).MassMetals[j];
#endif
            }
        }

      /* properties inside (radius of Vmax) */
      if(rr_tmp < 1.0 * vmaxrad)
        {
          massinmaxrad += P[p].Mass;
          massinmaxrad_tab[ptype] += P[p].Mass;

          if(ptype == 0)
            {
              if(P[p].Type == 0)
                {
#ifdef USE_SFR
                  sfrinmaxrad += SphP[PS[p].OldIndex].Sfr;      /* note: the SphP[] array has not been reordered */
#endif
#ifdef GFM_STELLAR_EVOLUTION
                  gasMassMetallicityMaxRad += SphP[PS[p].OldIndex].Metallicity * P[p].Mass;
                  for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
                    gasMassMetalsMaxRad[j] += SphP[PS[p].OldIndex].MetalsFraction[j] * P[p].Mass;

#ifdef GFM_DUST
                  for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
                    for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
                      gasMassDustMetallicityMaxRad += SphP[PS[p].OldIndex].MetalsDustFraction[chan][j] * P[p].Mass;
#endif
#endif
                }
              if(P[p].Type == 4)
                {
#ifdef GFM_STELLAR_EVOLUTION
                  gasMassMetallicityMaxRad += STP(p).Metallicity * P[p].Mass;
                  for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
                    gasMassMetalsMaxRad[j] += STP(p).MassMetals[j];

#ifdef GFM_DUST
                  for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
                    for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
                      gasMassDustMetallicityMaxRad += STP(p).MassMetals[j] * STP(p).InitialDustFractions[chan][j];
#endif
#endif
                }
            }
          if(ptype == 4)
            {
#ifdef GFM_STELLAR_EVOLUTION
              stellarMassMetallicityMaxRad += STP(p).Metallicity * P[p].Mass;
              for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
                stellarMassMetalsMaxRad[j] += STP(p).MassMetals[j];
#endif
            }
        }
    }

  /* properties of star forming gas */
#ifdef USE_SFR
  for(i = 0; i < num; i++)
    {
      p = d[i].index;

      if(P[p].Type == 0)
        {
          if(SphP[PS[p].OldIndex].Sfr > 0)
            {
              gasMassSfr += P[p].Mass;
#ifdef GFM_STELLAR_EVOLUTION
              gasMassMetallicitySfr += SphP[PS[p].OldIndex].Metallicity * P[p].Mass;
              for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
                gasMassMetalsSfr[j] += SphP[PS[p].OldIndex].MetalsFraction[j] * P[p].Mass;

              gasMassMetallicitySfrWeighted += SphP[PS[p].OldIndex].Metallicity * SphP[PS[p].OldIndex].Sfr;
              for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
                gasMassMetalsSfrWeighted[j] += SphP[PS[p].OldIndex].MetalsFraction[j] * SphP[PS[p].OldIndex].Sfr;

#ifdef GFM_DUST
              for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
                for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
                  gasMassDustMetallicitySfr += SphP[PS[p].OldIndex].MetalsDustFraction[chan][j] * P[p].Mass;

              for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
                for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
                  gasMassDustMetallicitySfrWeighted += SphP[PS[p].OldIndex].MetalsDustFraction[chan][j] * SphP[PS[p].OldIndex].Sfr;
#endif
#endif
            }
        }
    }
#endif


#if defined(MHD) && defined(ADD_MAGNETIC_GROUP_PROPERTIES)  /* compute magnetic half-energy radius */
  struct rlist_mhd
  {
    double r;
    double b_egy;
  };

  struct rlist_mhd  *my_mhd_list;

  /* let's allocate a list of our tuples */
  my_mhd_list = mymalloc("my_mhd_list", sizeof(struct rlist_mhd) * (num + 1));

  /* let's fill the list */

  for(int i = 0; i < num; i++)  /* loop over particles */
     {
       int p = d[i].index;     /* get index of particle */

       /* calculate distance to minimum of potential  */

       double xx = GRAVITY_NEAREST_X(P[p].Pos[0] - pos[0]);
       double yy = GRAVITY_NEAREST_Y(P[p].Pos[1] - pos[1]);
       double zz = GRAVITY_NEAREST_Z(P[p].Pos[2] - pos[2]);

       double rr = sqrt(xx * xx + yy * yy + zz * zz);

       double bfld2 = (SphP[PS[p].OldIndex].B[0] * SphP[PS[p].OldIndex].B[0]) + (SphP[PS[p].OldIndex].B[1] * SphP[PS[p].OldIndex].B[1]) + (SphP[PS[p].OldIndex].B[2] * SphP[PS[p].OldIndex].B[2]);
       double vol = SphP[PS[p].OldIndex].Volume;

       my_mhd_list[i].r = rr;
       my_mhd_list[i].b_egy = bfld2 * vol;
     }

   /* now sort by distance */
  if(parallel_flag)
     parallel_sort_comm(my_mhd_list, num, sizeof(struct rlist_mhd), subfind_compare_rlist_mhd, SubComm);
   else
     mysort(my_mhd_list, num, sizeof(struct rlist_mhd), subfind_compare_rlist_mhd);

  /* calculate cumulative magnetic egy */
   for(int i = 1; i < num; i++)
     my_mhd_list[i].b_egy += my_mhd_list[i - 1].b_egy;

  double egy_tot = 0;

  if(parallel_flag)
    {
      double egy_part = 0;
      if(num)
	egy_part = my_mhd_list[num - 1].b_egy;
      double *egylist = mymalloc("egylist", SubNTask * sizeof(double));
      MPI_Allgather(&egy_part, 1, MPI_DOUBLE, egylist, 1, MPI_DOUBLE, SubComm);

      double egybefore = 0;
      for(int i = 0; i < SubThisTask; i++)
	egybefore += egylist[i];

      for(int i = 0; i < SubNTask; i++)
      	egy_tot += egylist[i];

      for(int i = 0; i < num; i++)
	my_mhd_list[i].b_egy += egybefore;

      myfree(egylist);
    }
  else
    {
      egy_tot = my_mhd_list[num - 1].b_egy;  /* last enry and we know list is not empty */
    }

  /* now let's find half energy radius */
  double halfegyrad = 0;

  for(int i = 0; i < num - 1; i++)
    {
      if(my_mhd_list[i].b_egy < 0.5 * egy_tot && my_mhd_list[i+1].b_egy >= 0.5 * egy_tot)
	halfegyrad = 0.5 * (my_mhd_list[i].r + my_mhd_list[i + 1].r);
    }

  if(parallel_flag)
    {
      double halfegyrad_all;
      MPI_Allreduce(&halfegyrad, &halfegyrad_all, 1, MPI_DOUBLE, MPI_MAX, SubComm);
      halfegyrad = halfegyrad_all;
    }

  /* the result is now contained in halfegyrad */

  myfree(my_mhd_list);

#endif





#ifdef MHD
  bfld_halo = bfld_disk = bfld_vol_halo = bfld_vol_disk = 0;

  for(i = 0; i < num; i++)
    {
      p = d[i].index;

      if(P[p].Type == 0)
        {
          double bfld2 = (SphP[PS[p].OldIndex].B[0] * SphP[PS[p].OldIndex].B[0]) + (SphP[PS[p].OldIndex].B[1] * SphP[PS[p].OldIndex].B[1]) + (SphP[PS[p].OldIndex].B[2] * SphP[PS[p].OldIndex].B[2]);
          double vol = SphP[PS[p].OldIndex].Volume;

          bfld_halo += bfld2 * vol;
          bfld_vol_halo += vol;

          /* calculate particle radius */
          for(j = 0, rr_tmp = 0; j < 3; j++)
            {
              ddxx = GRAVITY_NEAREST_X(P[p].Pos[j] - pos[j]);
              rr_tmp += ddxx * ddxx;
            }
          rr_tmp = sqrt(rr_tmp);

          if(rr_tmp < SUBFIND_GAL_RADIUS_FAC * halfmassradtype[4])
            {
              bfld_disk += bfld2 * vol;
              bfld_vol_disk += vol;
            }
        }
    }
#endif

  if(parallel_flag)
    {
      double massinradtot, massinrad_tabtot[NTYPES];
      MPI_Allreduce(&massinrad, &massinradtot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(massinrad_tab, massinrad_tabtot, NTYPES, MPI_DOUBLE, MPI_SUM, SubComm);
      massinrad = massinradtot;
      for(j = 0; j < NTYPES; j++)
        massinrad_tab[j] = massinrad_tabtot[j];

      double massinhalfradtot, massinhalfrad_tabtot[NTYPES];
      MPI_Allreduce(&massinhalfrad, &massinhalfradtot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(massinhalfrad_tab, massinhalfrad_tabtot, NTYPES, MPI_DOUBLE, MPI_SUM, SubComm);
      massinhalfrad = massinhalfradtot;
      for(j = 0; j < NTYPES; j++)
        massinhalfrad_tab[j] = massinhalfrad_tabtot[j];

      double massinmaxradtot, massinmaxrad_tabtot[NTYPES];
      MPI_Allreduce(&massinmaxrad, &massinmaxradtot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(massinmaxrad_tab, massinmaxrad_tabtot, NTYPES, MPI_DOUBLE, MPI_SUM, SubComm);
      massinmaxrad = massinmaxradtot;
      for(j = 0; j < NTYPES; j++)
        massinmaxrad_tab[j] = massinmaxrad_tabtot[j];

#ifdef SUBFIND_EXTENDED_PROPERTIES
      MPI_Allreduce(MPI_IN_PLACE, sinrad, 3, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, vinrad, 3, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, sinhalfrad, 3, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, vinhalfrad, 3, MPI_DOUBLE, MPI_SUM, SubComm);
#endif

#ifdef MHD
      double bfld_halo_tot, bfld_disk_tot, bfld_vol_halo_tot, bfld_vol_disk_tot;
      MPI_Allreduce(&bfld_halo, &bfld_halo_tot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(&bfld_vol_halo, &bfld_vol_halo_tot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(&bfld_disk, &bfld_disk_tot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(&bfld_vol_disk, &bfld_vol_disk_tot, 1, MPI_DOUBLE, MPI_SUM, SubComm);

      bfld_halo = bfld_halo_tot;
      bfld_vol_halo = bfld_vol_halo_tot;
      bfld_disk = bfld_disk_tot;
      bfld_vol_disk = bfld_vol_disk_tot;
#endif

#ifdef USE_SFR
      double sfrinradtot;
      MPI_Allreduce(&sfrinrad, &sfrinradtot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      sfrinrad = sfrinradtot;

      double sfrinhalfradtot;
      MPI_Allreduce(&sfrinhalfrad, &sfrinhalfradtot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      sfrinhalfrad = sfrinhalfradtot;

      double sfrinmaxradtot;
      MPI_Allreduce(&sfrinmaxrad, &sfrinmaxradtot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      sfrinmaxrad = sfrinmaxradtot;

      double gasMassSfrtot;
      MPI_Allreduce(&gasMassSfr, &gasMassSfrtot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      gasMassSfr = gasMassSfrtot;
#endif
#ifdef GFM_STELLAR_EVOLUTION
      double gasMassMetallicitytot;
      double gasMassMetalstot[GFM_N_CHEM_ELEMENTS];
      MPI_Allreduce(&gasMassMetallicity, &gasMassMetallicitytot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      gasMassMetallicity = gasMassMetallicitytot;
      MPI_Allreduce(gasMassMetals, gasMassMetalstot, GFM_N_CHEM_ELEMENTS, MPI_DOUBLE, MPI_SUM, SubComm);
      for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
        gasMassMetals[j] = gasMassMetalstot[j];

      double stellarMassMetallicitytot;
      double stellarMassMetalstot[GFM_N_CHEM_ELEMENTS];
      MPI_Allreduce(&stellarMassMetallicity, &stellarMassMetallicitytot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      stellarMassMetallicity = stellarMassMetallicitytot;
      MPI_Allreduce(stellarMassMetals, stellarMassMetalstot, GFM_N_CHEM_ELEMENTS, MPI_DOUBLE, MPI_SUM, SubComm);
      for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
        stellarMassMetals[j] = stellarMassMetalstot[j];

      double gasMassMetallicityHalfRadtot;
      double gasMassMetalsHalfRadtot[GFM_N_CHEM_ELEMENTS];
      MPI_Allreduce(&gasMassMetallicityHalfRad, &gasMassMetallicityHalfRadtot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      gasMassMetallicityHalfRad = gasMassMetallicityHalfRadtot;
      MPI_Allreduce(gasMassMetalsHalfRad, gasMassMetalsHalfRadtot, GFM_N_CHEM_ELEMENTS, MPI_DOUBLE, MPI_SUM, SubComm);
      for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
        gasMassMetalsHalfRad[j] = gasMassMetalsHalfRadtot[j];

      double stellarMassMetallicityHalfRadtot;
      double stellarMassMetalsHalfRadtot[GFM_N_CHEM_ELEMENTS];
      MPI_Allreduce(&stellarMassMetallicityHalfRad, &stellarMassMetallicityHalfRadtot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      stellarMassMetallicityHalfRad = stellarMassMetallicityHalfRadtot;
      MPI_Allreduce(stellarMassMetalsHalfRad, stellarMassMetalsHalfRadtot, GFM_N_CHEM_ELEMENTS, MPI_DOUBLE, MPI_SUM, SubComm);
      for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
        stellarMassMetalsHalfRad[j] = stellarMassMetalsHalfRadtot[j];

      double gasMassMetallicityMaxRadtot;
      double gasMassMetalsMaxRadtot[GFM_N_CHEM_ELEMENTS];
      MPI_Allreduce(&gasMassMetallicityMaxRad, &gasMassMetallicityMaxRadtot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      gasMassMetallicityMaxRad = gasMassMetallicityMaxRadtot;
      MPI_Allreduce(gasMassMetalsMaxRad, gasMassMetalsMaxRadtot, GFM_N_CHEM_ELEMENTS, MPI_DOUBLE, MPI_SUM, SubComm);
      for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
        gasMassMetalsMaxRad[j] = gasMassMetalsMaxRadtot[j];

      double stellarMassMetallicityMaxRadtot;
      double stellarMassMetalsMaxRadtot[GFM_N_CHEM_ELEMENTS];
      MPI_Allreduce(&stellarMassMetallicityMaxRad, &stellarMassMetallicityMaxRadtot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      stellarMassMetallicityMaxRad = stellarMassMetallicityMaxRadtot;
      MPI_Allreduce(stellarMassMetalsMaxRad, stellarMassMetalsMaxRadtot, GFM_N_CHEM_ELEMENTS, MPI_DOUBLE, MPI_SUM, SubComm);
      for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
        stellarMassMetalsMaxRad[j] = stellarMassMetalsMaxRadtot[j];

#ifdef GFM_DUST
      double gasMassDustMetallicitytot;
      MPI_Allreduce(&gasMassDustMetallicity, &gasMassDustMetallicitytot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      gasMassDustMetallicity = gasMassDustMetallicitytot;

      double gasMassDustMetallicityHalfRadtot;
      MPI_Allreduce(&gasMassDustMetallicityHalfRad, &gasMassDustMetallicityHalfRadtot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      gasMassDustMetallicityHalfRad = gasMassDustMetallicityHalfRadtot;

      double gasMassDustMetallicityMaxRadtot;
      MPI_Allreduce(&gasMassDustMetallicityMaxRad, &gasMassDustMetallicityMaxRadtot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      gasMassDustMetallicityMaxRad = gasMassDustMetallicityMaxRadtot;
#endif

#ifdef USE_SFR
      double gasMassMetallicitySfrtot;
      double gasMassMetalsSfrtot[GFM_N_CHEM_ELEMENTS];
      MPI_Allreduce(&gasMassMetallicitySfr, &gasMassMetallicitySfrtot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      gasMassMetallicitySfr = gasMassMetallicitySfrtot;
      MPI_Allreduce(gasMassMetalsSfr, gasMassMetalsSfrtot, GFM_N_CHEM_ELEMENTS, MPI_DOUBLE, MPI_SUM, SubComm);
      for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
        gasMassMetalsSfr[j] = gasMassMetalsSfrtot[j];

      double gasMassMetallicitySfrWeightedtot;
      double gasMassMetalsSfrWeightedtot[GFM_N_CHEM_ELEMENTS];
      MPI_Allreduce(&gasMassMetallicitySfrWeighted, &gasMassMetallicitySfrWeightedtot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      gasMassMetallicitySfrWeighted = gasMassMetallicitySfrWeightedtot;
      MPI_Allreduce(gasMassMetalsSfrWeighted, gasMassMetalsSfrWeightedtot, GFM_N_CHEM_ELEMENTS, MPI_DOUBLE, MPI_SUM, SubComm);
      for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
        gasMassMetalsSfrWeighted[j] = gasMassMetalsSfrWeightedtot[j];

#ifdef GFM_DUST
      double gasMassDustMetallicitySfrtot;
      MPI_Allreduce(&gasMassDustMetallicitySfr, &gasMassDustMetallicitySfrtot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      gasMassDustMetallicitySfr = gasMassDustMetallicitySfrtot;

      double gasMassDustMetallicitySfrWeightedtot;
      MPI_Allreduce(&gasMassDustMetallicitySfrWeighted, &gasMassDustMetallicitySfrWeightedtot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      gasMassDustMetallicitySfrWeighted = gasMassDustMetallicitySfrWeightedtot;
#endif

#endif
#endif
    }

  if(parallel_flag)
    MPI_Allreduce(&num, &totlen, 1, MPI_INT, MPI_SUM, SubComm);
  else
    totlen = num;

#ifdef MHD
  if(bfld_vol_halo > 0.)
    bfld_halo = sqrt(bfld_halo / bfld_vol_halo);
  if(bfld_vol_disk > 0.)
    bfld_disk = sqrt(bfld_disk / bfld_vol_disk);
#endif

#ifdef SUBFIND_EXTENDED_PROPERTIES
  // finish centre of mass of spheres
  for(j = 0; j < 3; j++)
    {
      if(massinrad > 0)
        {
          sinrad[j] /= massinrad;
          sinrad[j] += pos[j];

#ifdef PERIODIC
          while(sinrad[j] < 0)
            sinrad[j] += boxsize;
          while(sinrad[j] >= boxsize)
            sinrad[j] -= boxsize;
#endif

          vinrad[j] /= massinrad;       // this is comoving (as it should be.)
        }

      if(massinhalfrad > 0)
        {
          sinhalfrad[j] /= massinhalfrad;
          sinhalfrad[j] += pos[j];

#ifdef PERIODIC
          while(sinhalfrad[j] < 0)
            sinhalfrad[j] += boxsize;
          while(sinhalfrad[j] >= boxsize)
            sinhalfrad[j] -= boxsize;
#endif

          vinhalfrad[j] /= massinhalfrad;
        }
    }

  Jtot_inHalfRad[0] = Jtot_inHalfRad[1] = Jtot_inHalfRad[2] = 0;
  Jdm_inHalfRad[0] = Jdm_inHalfRad[1] = Jdm_inHalfRad[2] = 0;
  Jgas_inHalfRad[0] = Jgas_inHalfRad[1] = Jgas_inHalfRad[2] = 0;
  Jstars_inHalfRad[0] = Jstars_inHalfRad[1] = Jstars_inHalfRad[2] = 0;
  Jtot_inRad[0] = Jtot_inRad[1] = Jtot_inRad[2] = 0;
  Jdm_inRad[0] = Jdm_inRad[1] = Jdm_inRad[2] = 0;
  Jgas_inRad[0] = Jgas_inRad[1] = Jgas_inRad[2] = 0;
  Jstars_inRad[0] = Jstars_inRad[1] = Jstars_inRad[2] = 0;

  for(i = 0; i < num; i++)
    {
      /* identify particle type */
      p = d[i].index;

      /* calculate particle radius */
      for(j = 0, rr_tmp = 0; j < 3; j++)
        {
          ddxx = GRAVITY_NEAREST_X(P[p].Pos[j] - pos[j]);
          rr_tmp += ddxx * ddxx;
        }
      rr_tmp = sqrt(rr_tmp);

      int ptype = P[p].Type;
#ifdef GFM_WINDS
      /* count wind as gas for mass, but not for LenType since we use this to construct offset tables */
      if(P[p].Type == 4 && STP(p).BirthTime < 0)
        ptype = 0;
#endif

      /* properties inside SUBFIND_GAL_RADIUS_FAC*(stellar half mass radius) */
      if((massinrad > 0) && (rr_tmp < SUBFIND_GAL_RADIUS_FAC * halfmassradtype[4]))
        {
          for(j = 0; j < 3; j++)
            {
              ddxx = GRAVITY_NEAREST_X(P[p].Pos[j] - pos[j]);
              dx[j] = All.cf_atime * ddxx;
              dv[j] = vel_to_phys * (P[p].Vel[j] - vinrad[j]);
              dv[j] += H_of_a * dx[j];
            }

          Jtot_inRad[0] += P[p].Mass * (dx[1] * dv[2] - dx[2] * dv[1]);
          Jtot_inRad[1] += P[p].Mass * (dx[2] * dv[0] - dx[0] * dv[2]);
          Jtot_inRad[2] += P[p].Mass * (dx[0] * dv[1] - dx[1] * dv[0]);

          if(ptype == 1)        // dm illustris
            {
              Jdm_inRad[0] += P[p].Mass * (dx[1] * dv[2] - dx[2] * dv[1]);
              Jdm_inRad[1] += P[p].Mass * (dx[2] * dv[0] - dx[0] * dv[2]);
              Jdm_inRad[2] += P[p].Mass * (dx[0] * dv[1] - dx[1] * dv[0]);
            }
          if(ptype == 0)        // gas
            {
              Jgas_inRad[0] += P[p].Mass * (dx[1] * dv[2] - dx[2] * dv[1]);
              Jgas_inRad[1] += P[p].Mass * (dx[2] * dv[0] - dx[0] * dv[2]);
              Jgas_inRad[2] += P[p].Mass * (dx[0] * dv[1] - dx[1] * dv[0]);
            }
          if(ptype == 4)        // stars
            {
              Jstars_inRad[0] += P[p].Mass * (dx[1] * dv[2] - dx[2] * dv[1]);
              Jstars_inRad[1] += P[p].Mass * (dx[2] * dv[0] - dx[0] * dv[2]);
              Jstars_inRad[2] += P[p].Mass * (dx[0] * dv[1] - dx[1] * dv[0]);
            }
        }


      /* properties inside (stellar half mass radius) */
      if((massinhalfrad > 0) && (rr_tmp < 1.0 * halfmassradtype[4]))
        {
          for(j = 0; j < 3; j++)
            {
              ddxx = GRAVITY_NEAREST_X(P[p].Pos[j] - pos[j]);
              dx[j] = All.cf_atime * ddxx;
              dv[j] = vel_to_phys * (P[p].Vel[j] - vinhalfrad[j]);
              dv[j] += H_of_a * dx[j];
            }

          Jtot_inHalfRad[0] += P[p].Mass * (dx[1] * dv[2] - dx[2] * dv[1]);
          Jtot_inHalfRad[1] += P[p].Mass * (dx[2] * dv[0] - dx[0] * dv[2]);
          Jtot_inHalfRad[2] += P[p].Mass * (dx[0] * dv[1] - dx[1] * dv[0]);

          if(ptype == 1)        // dm illustris
            {
              Jdm_inHalfRad[0] += P[p].Mass * (dx[1] * dv[2] - dx[2] * dv[1]);
              Jdm_inHalfRad[1] += P[p].Mass * (dx[2] * dv[0] - dx[0] * dv[2]);
              Jdm_inHalfRad[2] += P[p].Mass * (dx[0] * dv[1] - dx[1] * dv[0]);
            }
          if(ptype == 0)        // gas
            {
              Jgas_inHalfRad[0] += P[p].Mass * (dx[1] * dv[2] - dx[2] * dv[1]);
              Jgas_inHalfRad[1] += P[p].Mass * (dx[2] * dv[0] - dx[0] * dv[2]);
              Jgas_inHalfRad[2] += P[p].Mass * (dx[0] * dv[1] - dx[1] * dv[0]);
            }
          if(ptype == 4)        // stars
            {
              Jstars_inHalfRad[0] += P[p].Mass * (dx[1] * dv[2] - dx[2] * dv[1]);
              Jstars_inHalfRad[1] += P[p].Mass * (dx[2] * dv[0] - dx[0] * dv[2]);
              Jstars_inHalfRad[2] += P[p].Mass * (dx[0] * dv[1] - dx[1] * dv[0]);
            }
        }
    }


  if(parallel_flag)
    {
      MPI_Allreduce(MPI_IN_PLACE, Jtot_inRad, 3, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, Jdm_inRad, 3, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, Jgas_inRad, 3, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, Jstars_inRad, 3, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, Jtot_inHalfRad, 3, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, Jdm_inHalfRad, 3, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, Jgas_inHalfRad, 3, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, Jstars_inHalfRad, 3, MPI_DOUBLE, MPI_SUM, SubComm);
    }

  // counter rotating mass fractions
  CMFrac_inHalfRad = CMFrac_inRad = 0;
  for(i = 0; i < NTYPES; i++)
    CMFracType_inHalfRad[i] = CMFracType_inRad[i] = 0;

  for(i = 0; i < num; i++)
    {
      /* identify particle type */
      p = d[i].index;

      /* calculate particle radius */
      for(j = 0, rr_tmp = 0; j < 3; j++)
        {
          ddxx = GRAVITY_NEAREST_X(P[p].Pos[j] - pos[j]);       // counter-rotating mass calc with respect to pot min
          rr_tmp += ddxx * ddxx;
        }
      rr_tmp = sqrt(rr_tmp);

      int ptype = P[p].Type;
#ifdef GFM_WINDS
      /* count wind as gas for mass, but not for LenType since we use this to construct offset tables */
      if(P[p].Type == 4 && STP(p).BirthTime < 0)
        ptype = 0;
#endif

      /* properties inside SUBFIND_GAL_RADIUS_FAC*(stellar half mass radius) */
      if((massinrad > 0) && (rr_tmp < SUBFIND_GAL_RADIUS_FAC * halfmassradtype[4]))
        {
          for(j = 0; j < 3; j++)
            {
              ddxx = GRAVITY_NEAREST_X(P[p].Pos[j] - pos[j]);
              dx[j] = All.cf_atime * ddxx;
              dv[j] = vel_to_phys * (P[p].Vel[j] - vinrad[j]);
              dv[j] += H_of_a * dx[j];
            }

          jpart[0] = P[p].Mass * (dx[1] * dv[2] - dx[2] * dv[1]);
          jpart[1] = P[p].Mass * (dx[2] * dv[0] - dx[0] * dv[2]);
          jpart[2] = P[p].Mass * (dx[0] * dv[1] - dx[1] * dv[0]);

          if((Jtot_inRad[0] * jpart[0] + Jtot_inRad[1] * jpart[1] + Jtot_inRad[2] * jpart[2]) < 0.)
            CMFrac_inRad += P[p].Mass / massinrad;

          if(ptype == 1)        // dm illustris
            if((Jdm_inRad[0] * jpart[0] + Jdm_inRad[1] * jpart[1] + Jdm_inRad[2] * jpart[2]) < 0.)
              CMFracType_inRad[1] += P[p].Mass / massinrad_tab[1];
          if(ptype == 0)        // gas (incl. winds!)
            if((Jgas_inRad[0] * jpart[0] + Jgas_inRad[1] * jpart[1] + Jgas_inRad[2] * jpart[2]) < 0.)
              CMFracType_inRad[0] += P[p].Mass / massinrad_tab[0];
          if(ptype == 4)        // stars
            if((Jstars_inRad[0] * jpart[0] + Jstars_inRad[1] * jpart[1] + Jstars_inRad[2] * jpart[2]) < 0.)
              CMFracType_inRad[4] += P[p].Mass / massinrad_tab[4];
        }

      /* properties inside (stellar half mass radius) */
      if((massinhalfrad > 0) && (rr_tmp < 1.0 * halfmassradtype[4]))
        {
          for(j = 0; j < 3; j++)
            {
              ddxx = GRAVITY_NEAREST_X(P[p].Pos[j] - pos[j]);
              dx[j] = All.cf_atime * ddxx;
              dv[j] = vel_to_phys * (P[p].Vel[j] - vinhalfrad[j]);
              dv[j] += H_of_a * dx[j];
            }

          jpart[0] = P[p].Mass * (dx[1] * dv[2] - dx[2] * dv[1]);
          jpart[1] = P[p].Mass * (dx[2] * dv[0] - dx[0] * dv[2]);
          jpart[2] = P[p].Mass * (dx[0] * dv[1] - dx[1] * dv[0]);

          if((Jtot_inHalfRad[0] * jpart[0] + Jtot_inHalfRad[1] * jpart[1] + Jtot_inHalfRad[2] * jpart[2]) < 0.)
            CMFrac_inHalfRad += P[p].Mass / massinhalfrad;

          if(ptype == 1)        // dm illustris
            if((Jdm_inHalfRad[0] * jpart[0] + Jdm_inHalfRad[1] * jpart[1] + Jdm_inHalfRad[2] * jpart[2]) < 0.)
              CMFracType_inHalfRad[1] += P[p].Mass / massinhalfrad_tab[1];
          if(ptype == 0)        // gas (incl. winds!)
            if((Jgas_inHalfRad[0] * jpart[0] + Jgas_inHalfRad[1] * jpart[1] + Jgas_inHalfRad[2] * jpart[2]) < 0.)
              CMFracType_inHalfRad[0] += P[p].Mass / massinhalfrad_tab[0];
          if(ptype == 4)        // stars
            if((Jstars_inHalfRad[0] * jpart[0] + Jstars_inHalfRad[1] * jpart[1] + Jstars_inHalfRad[2] * jpart[2]) < 0.)
              CMFracType_inHalfRad[4] += P[p].Mass / massinhalfrad_tab[4];
        }
    }

  if(parallel_flag)
    {
      MPI_Allreduce(MPI_IN_PLACE, &CMFrac_inRad, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, &CMFrac_inHalfRad, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, CMFracType_inRad, NTYPES, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(MPI_IN_PLACE, CMFracType_inHalfRad, NTYPES, MPI_DOUBLE, MPI_SUM, SubComm);
    }
#endif



#ifdef SUBFIND_MEASURE_H2MASS
  double H2_Mass = 0;

  for(i = 0; i < num; i++)
    {
      p = d[i].index;
      if(P[p].Type == 0)
        {
          if(SphP[PS[p].OldIndex].Sfr > 0)
            {
              double MeanWeight = 4.0 / (3 * HYDROGEN_MASSFRAC + 1 + 4 * HYDROGEN_MASSFRAC * SphP[PS[p].OldIndex].Ne) * PROTONMASS;
              double Temp = MeanWeight / BOLTZMANN * GAMMA_MINUS1 * SphP[PS[p].OldIndex].Utherm * All.UnitEnergy_in_cgs / All.UnitMass_in_g;
              double nh = SphP[PS[p].OldIndex].Density * ((All.UnitMass_in_g / All.HubbleParam) / pow(All.Time, 3) / pow(All.UnitLength_in_cm / All.HubbleParam, 3) * HYDROGEN_MASSFRAC / PROTONMASS);  //FIXME?
              if(Temp < pow(10.0, 4.2) * (nh / 0.1))
                {
                  double Rsurf = pow(0.5 * nh * Temp / pow(10.0, 4.3), 0.8);
                  double fac = Rsurf / (1 + Rsurf);

                  H2_Mass += fac * HYDROGEN_MASSFRAC * P[p].Mass;
                }
            }
        }
    }

  if(parallel_flag)
    {
      double H2_Mass_tot;
      MPI_Allreduce(&H2_Mass, &H2_Mass_tot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      H2_Mass = H2_Mass_tot;
    }
#endif

#ifdef GFM_STELLAR_PHOTOMETRICS
  double Magnitude_U, Magnitude_B, Magnitude_V, Magnitude_K;
  double Magnitude_g, Magnitude_r, Magnitude_i, Magnitude_z;
  double Luminosity_U = 0, Luminosity_B = 0, Luminosity_V = 0, Luminosity_K = 0;
  double Luminosity_g = 0, Luminosity_r = 0, Luminosity_i = 0, Luminosity_z = 0;
  stellar_photometrics st_photo;

  for(i = 0; i < num; i++)
    {
      p = d[i].index;
      if(P[p].Type == 4 && STP(p).BirthTime >= 0)
        {
          assign_stellar_photometrics(p, &st_photo);
          Luminosity_U += pow(10.0, -0.4 * st_photo.Magnitude_U);
          Luminosity_B += pow(10.0, -0.4 * st_photo.Magnitude_B);
          Luminosity_V += pow(10.0, -0.4 * st_photo.Magnitude_V);
          Luminosity_K += pow(10.0, -0.4 * st_photo.Magnitude_K);
          Luminosity_g += pow(10.0, -0.4 * st_photo.Magnitude_g);
          Luminosity_r += pow(10.0, -0.4 * st_photo.Magnitude_r);
          Luminosity_i += pow(10.0, -0.4 * st_photo.Magnitude_i);
          Luminosity_z += pow(10.0, -0.4 * st_photo.Magnitude_z);
        }
    }

  if(parallel_flag)
    {
      double Luminosity_U_tot, Luminosity_B_tot, Luminosity_V_tot, Luminosity_K_tot;
      double Luminosity_g_tot, Luminosity_r_tot, Luminosity_i_tot, Luminosity_z_tot;
      MPI_Allreduce(&Luminosity_U, &Luminosity_U_tot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(&Luminosity_B, &Luminosity_B_tot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(&Luminosity_V, &Luminosity_V_tot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(&Luminosity_K, &Luminosity_K_tot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(&Luminosity_g, &Luminosity_g_tot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(&Luminosity_r, &Luminosity_r_tot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(&Luminosity_i, &Luminosity_i_tot, 1, MPI_DOUBLE, MPI_SUM, SubComm);
      MPI_Allreduce(&Luminosity_z, &Luminosity_z_tot, 1, MPI_DOUBLE, MPI_SUM, SubComm);


      Luminosity_U = Luminosity_U_tot;
      Luminosity_B = Luminosity_B_tot;
      Luminosity_V = Luminosity_V_tot;
      Luminosity_K = Luminosity_K_tot;
      Luminosity_g = Luminosity_g_tot;
      Luminosity_r = Luminosity_r_tot;
      Luminosity_i = Luminosity_i_tot;
      Luminosity_z = Luminosity_z_tot;

    }

  if(Luminosity_U > 0 && Luminosity_B > 0 && Luminosity_V > 0 && Luminosity_K > 0 && Luminosity_g > 0 && Luminosity_r > 0 && Luminosity_i > 0 && Luminosity_z)
    {
      Magnitude_U = -2.5 * log10(Luminosity_U);
      Magnitude_B = -2.5 * log10(Luminosity_B);
      Magnitude_V = -2.5 * log10(Luminosity_V);
      Magnitude_K = -2.5 * log10(Luminosity_K);
      Magnitude_g = -2.5 * log10(Luminosity_g);
      Magnitude_r = -2.5 * log10(Luminosity_r);
      Magnitude_i = -2.5 * log10(Luminosity_i);
      Magnitude_z = -2.5 * log10(Luminosity_z);
    }
  else
    {
      Magnitude_U = MAX_FLOAT_NUMBER;
      Magnitude_B = MAX_FLOAT_NUMBER;
      Magnitude_V = MAX_FLOAT_NUMBER;
      Magnitude_K = MAX_FLOAT_NUMBER;
      Magnitude_g = MAX_FLOAT_NUMBER;
      Magnitude_r = MAX_FLOAT_NUMBER;
      Magnitude_i = MAX_FLOAT_NUMBER;
      Magnitude_z = MAX_FLOAT_NUMBER;
    }

  /* calculate surface brightness profile and find the maximum radius that is still above the threshold */
  int idir, irad;
  double rotxy[2], rxy;
  double sintheta, costheta, sinphi, cosphi;
  double radii_grid[GFM_STELLAR_PHOTOMETRICS_RADII];
  MyFloat luminosity_grid[GFM_STELLAR_PHOTOMETRICS_DIRECTIONS][GFM_STELLAR_PHOTOMETRICS_RADII - 1];
  double mass_grid[GFM_STELLAR_PHOTOMETRICS_DIRECTIONS][GFM_STELLAR_PHOTOMETRICS_RADII - 1];
  MyFloat tmp_surface_brightness;
  sort_r2list brightness_limit_list[GFM_STELLAR_PHOTOMETRICS_DIRECTIONS];

  memset(luminosity_grid, 0, GFM_STELLAR_PHOTOMETRICS_DIRECTIONS * (GFM_STELLAR_PHOTOMETRICS_RADII - 1) * sizeof(MyFloat));
  memset(mass_grid, 0, GFM_STELLAR_PHOTOMETRICS_DIRECTIONS * (GFM_STELLAR_PHOTOMETRICS_RADII - 1) * sizeof(double));
  memset(brightness_limit_list, 0, GFM_STELLAR_PHOTOMETRICS_DIRECTIONS * sizeof(sort_r2list));

  if(mass_tab[4] > 0)
    {
      radii_grid[0] = 0;
      for(irad = 1; irad < GFM_STELLAR_PHOTOMETRICS_RADII; irad++)
        {
          if(get_default_softening_of_particletype(4) < max_stellar_rad)
            radii_grid[irad] = get_default_softening_of_particletype(4) * pow(max_stellar_rad / get_default_softening_of_particletype(4), (double) (irad - 1) / (GFM_STELLAR_PHOTOMETRICS_RADII - 2));
          else
            radii_grid[irad] = max_stellar_rad;
        }

      for(idir = 0; idir < GFM_STELLAR_PHOTOMETRICS_DIRECTIONS; idir++)
        {
          sintheta = sin(StellarPhotometricsRandomAngles[idir][0]);
          costheta = cos(StellarPhotometricsRandomAngles[idir][0]);
          sinphi = sin(StellarPhotometricsRandomAngles[idir][1]);
          cosphi = cos(StellarPhotometricsRandomAngles[idir][1]);

          for(i = 0; i < num; i++)
            {
              p = d[i].index;
              if(P[p].Type == 4 && STP(p).BirthTime >= 0)
                {
                  rotxy[0] = GRAVITY_NEAREST_X(P[p].Pos[0] - pos[0]) * cosphi + GRAVITY_NEAREST_Y(P[p].Pos[1] - pos[1]) * sinphi;
                  rotxy[1] =
                    (-1) * GRAVITY_NEAREST_X(P[p].Pos[0] - pos[0]) * costheta * sinphi + GRAVITY_NEAREST_Y(P[p].Pos[1] - pos[1]) * costheta * cosphi +
                    GRAVITY_NEAREST_Z(P[p].Pos[2] - pos[2]) * sintheta;
                  rxy = sqrt(rotxy[0] * rotxy[0] + rotxy[1] * rotxy[1]);
                  if(rxy > 0 && rxy >= max_stellar_rad)
                    terminate("SUBFIND GFM_STELLAR_PHOTOMETRICS: rxy=%g >= max_stellar_rad=%g idir=%d orig-pos=[%g %g %g]", rxy, max_stellar_rad,
                              idir, GRAVITY_NEAREST_X(P[p].Pos[0] - pos[0]), GRAVITY_NEAREST_Y(P[p].Pos[1] - pos[1]), GRAVITY_NEAREST_Z(P[p].Pos[2] - pos[2]));

                  for(irad = 0; irad < GFM_STELLAR_PHOTOMETRICS_RADII - 1; irad++)
                    {
                      if((rxy == 0 && irad == 0) || (rxy > radii_grid[irad] && rxy <= radii_grid[irad + 1]))
                        {
                          assign_stellar_photometrics(p, &st_photo);
                          luminosity_grid[idir][irad] += pow(10.0, -0.4 * st_photo.Magnitude_K);
                          mass_grid[idir][irad] += P[p].Mass;
                          break;
                        }
                    }
                }
            }
        }

      if(parallel_flag)
        {
          MyFloat luminosity_grid_tot[GFM_STELLAR_PHOTOMETRICS_DIRECTIONS][GFM_STELLAR_PHOTOMETRICS_RADII - 1];
          MPI_Allreduce(&luminosity_grid[0][0], &luminosity_grid_tot[0][0], GFM_STELLAR_PHOTOMETRICS_DIRECTIONS * (GFM_STELLAR_PHOTOMETRICS_RADII - 1), MPI_MYFLOAT, MPI_SUM, SubComm);
          double mass_grid_tot[GFM_STELLAR_PHOTOMETRICS_DIRECTIONS][GFM_STELLAR_PHOTOMETRICS_RADII - 1];
          MPI_Allreduce(&mass_grid[0][0], &mass_grid_tot[0][0], GFM_STELLAR_PHOTOMETRICS_DIRECTIONS * (GFM_STELLAR_PHOTOMETRICS_RADII - 1), MPI_DOUBLE, MPI_SUM, SubComm);
          for(idir = 0; idir < GFM_STELLAR_PHOTOMETRICS_DIRECTIONS; idir++)
            for(irad = 0; irad < GFM_STELLAR_PHOTOMETRICS_RADII - 1; irad++)
              {
                luminosity_grid[idir][irad] = luminosity_grid_tot[idir][irad];
                mass_grid[idir][irad] = mass_grid_tot[idir][irad];
              }
        }

      for(idir = 0; idir < GFM_STELLAR_PHOTOMETRICS_DIRECTIONS; idir++)
        {
          brightness_limit_list[idir].r = radii_grid[1];
          for(irad = GFM_STELLAR_PHOTOMETRICS_RADII - 2; irad > 0; irad--)
            {
              /* convert luminosity/pc^2 to mag/arcsec^2, including cosmological surface brightness dimming */
              tmp_surface_brightness = 5 * log10(3600 * 360 / 10 / 2 / M_PI) - 2.5 * log10(luminosity_grid[idir][irad])
                + 2.5 * log10(M_PI * (pow(radii_grid[irad + 1], 2) - pow(radii_grid[irad], 2)))
                + 2.5 * log10(pow(All.cf_atime * All.UnitLength_in_cm / PARSEC / All.HubbleParam, 2)) - 2.5 * log10(pow(All.cf_atime, 4));
              if(luminosity_grid[idir][irad] > 0 && GFM_STELLAR_PHOTOMETRICS_K_LIMIT > tmp_surface_brightness)
                {
                  brightness_limit_list[idir].r = radii_grid[irad + 1];
                  break;
                }
            }
          for(; irad >= 0; irad--)
            brightness_limit_list[idir].mass += mass_grid[idir][irad];
        }

      mysort(brightness_limit_list, GFM_STELLAR_PHOTOMETRICS_DIRECTIONS, sizeof(sort_r2list), subfind_compare_dist_rotcurve);

      brightness_limit_rad = brightness_limit_list[(int) (GFM_STELLAR_PHOTOMETRICS_DIRECTIONS / 2)].r;
      stellar_mass_in_phot_rad = brightness_limit_list[(int) (GFM_STELLAR_PHOTOMETRICS_DIRECTIONS / 2)].mass;
    }
#endif


  /* now store the calculated properties in the subgroup structure */
  if(parallel_flag == 0 || SubThisTask == 0)
    {
      subgroup->Len = totlen;
      subgroup->Mass = mass;
      subgroup->SubMassInRad = massinrad;
      subgroup->SubMassInHalfRad = massinhalfrad;
      subgroup->SubMassInMaxRad = massinmaxrad;
#ifdef SUBFIND_EXTENDED_PROPERTIES
      subgroup->Ekin = Ekin;
      subgroup->Epot = Epot;
      subgroup->Ethr = Ethr;
      subgroup->CMFrac = CMFrac;
      subgroup->CMFrac_inHalfRad = CMFrac_inHalfRad;
      subgroup->CMFrac_inRad = CMFrac_inRad;
#endif

#ifdef MHD
      subgroup->Bfld_Halo = bfld_halo;
      subgroup->Bfld_Disk = bfld_disk;
#endif

      for(j = 0; j < 6; j++)
        {
          subgroup->MassType[j] = mass_tab[j];
          subgroup->LenType[j] = len_type[j];
          subgroup->SubHalfMassRadType[j] = halfmassradtype[j];
          subgroup->SubMassInRadType[j] = massinrad_tab[j];
          subgroup->SubMassInHalfRadType[j] = massinhalfrad_tab[j];
          subgroup->SubMassInMaxRadType[j] = massinmaxrad_tab[j];
#ifdef SUBFIND_EXTENDED_PROPERTIES
          subgroup->CMFracType[j] = CMFracType[j];
          subgroup->CMFracType_inHalfRad[j] = CMFracType_inHalfRad[j];
          subgroup->CMFracType_inRad[j] = CMFracType_inRad[j];
#endif
        }
      for(j = 0; j < 3; j++)
        {
          subgroup->Pos[j] = pos[j];
          subgroup->Vel[j] = vel[j];
          subgroup->CM[j] = cm[j];
          subgroup->Spin[j] = spin[j];
#ifdef SUBFIND_EXTENDED_PROPERTIES
          subgroup->J[j] = Jtot[j];
          subgroup->Jdm[j] = Jdm[j];
          subgroup->Jgas[j] = Jgas[j];
          subgroup->Jstars[j] = Jstars[j];
          subgroup->J_inHalfRad[j] = Jtot_inHalfRad[j];
          subgroup->Jdm_inHalfRad[j] = Jdm_inHalfRad[j];
          subgroup->Jgas_inHalfRad[j] = Jgas_inHalfRad[j];
          subgroup->Jstars_inHalfRad[j] = Jstars_inHalfRad[j];
          subgroup->J_inRad[j] = Jtot_inRad[j];
          subgroup->Jdm_inRad[j] = Jdm_inRad[j];
          subgroup->Jgas_inRad[j] = Jgas_inRad[j];
          subgroup->Jstars_inRad[j] = Jstars_inRad[j];
#endif
        }

      subgroup->SubMostBoundID = mostboundid;
      subgroup->SubVelDisp = veldisp;
      subgroup->SubVmax = vmax;
      subgroup->SubVmaxRad = vmaxrad;
      subgroup->SubHalfMassRad = halfmassrad;

#ifdef USE_SFR
      subgroup->Sfr = sfr;
      subgroup->SfrInRad = sfrinrad;
      subgroup->SfrInHalfRad = sfrinhalfrad;
      subgroup->SfrInMaxRad = sfrinmaxrad;
      subgroup->GasMassSfr = gasMassSfr;
#endif
#ifdef GFM_STELLAR_EVOLUTION
      subgroup->GasMassMetallicity = gasMassMetallicity;
      subgroup->StellarMassMetallicity = stellarMassMetallicity;
      subgroup->GasMassMetallicityHalfRad = gasMassMetallicityHalfRad;
      subgroup->StellarMassMetallicityHalfRad = stellarMassMetallicityHalfRad;
      subgroup->GasMassMetallicityMaxRad = gasMassMetallicityMaxRad;
      subgroup->StellarMassMetallicityMaxRad = stellarMassMetallicityMaxRad;
#ifdef GFM_DUST
      subgroup->GasMassDustMetallicity = gasMassDustMetallicity;
      subgroup->GasMassDustMetallicityHalfRad = gasMassDustMetallicityHalfRad;
      subgroup->GasMassDustMetallicityMaxRad = gasMassDustMetallicityMaxRad;
#endif
#ifdef USE_SFR
      subgroup->GasMassMetallicitySfr = gasMassMetallicitySfr;
      subgroup->GasMassMetallicitySfrWeighted = gasMassMetallicitySfrWeighted;
#ifdef GFM_DUST
      subgroup->GasMassDustMetallicitySfr = gasMassDustMetallicitySfr;
      subgroup->GasMassDustMetallicitySfrWeighted = gasMassDustMetallicitySfrWeighted;
#endif
#endif
      for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
        {
          subgroup->GasMassMetals[j] = gasMassMetals[j];
          subgroup->StellarMassMetals[j] = stellarMassMetals[j];
          subgroup->GasMassMetalsHalfRad[j] = gasMassMetalsHalfRad[j];
          subgroup->StellarMassMetalsHalfRad[j] = stellarMassMetalsHalfRad[j];
          subgroup->GasMassMetalsMaxRad[j] = gasMassMetalsMaxRad[j];
          subgroup->StellarMassMetalsMaxRad[j] = stellarMassMetalsMaxRad[j];
#ifdef USE_SFR
          subgroup->GasMassMetalsSfr[j] = gasMassMetalsSfr[j];
          subgroup->GasMassMetalsSfrWeighted[j] = gasMassMetalsSfrWeighted[j];
#endif
        }
#endif
#ifdef BLACK_HOLES
      subgroup->BH_Mass = bh_Mass;
      subgroup->BH_Mdot = bh_Mdot;
#endif
#ifdef GFM_WINDS
      subgroup->WindMass = windMass;
#endif
#ifdef SUBFIND_MEASURE_H2MASS
      subgroup->H2_Mass = H2_Mass;
#endif
#ifdef GFM_STELLAR_PHOTOMETRICS
      subgroup->Magnitude_U = Magnitude_U;
      subgroup->Magnitude_B = Magnitude_B;
      subgroup->Magnitude_V = Magnitude_V;
      subgroup->Magnitude_K = Magnitude_K;
      subgroup->Magnitude_g = Magnitude_g;
      subgroup->Magnitude_r = Magnitude_r;
      subgroup->Magnitude_i = Magnitude_i;
      subgroup->Magnitude_z = Magnitude_z;
      subgroup->SurfaceBrightnessLimitRad = brightness_limit_rad;
      subgroup->SubMassInPhotRad = stellar_mass_in_phot_rad;
#endif
    }
}




#endif
