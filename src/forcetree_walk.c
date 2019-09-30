/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/forcetree_walk.c
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

#include "allvars.h"
#include "proto.h"
#if defined(RADCOOL) && !defined(GFM)
#error "Need to compile with GFM if RADCOOL option is chosen"
#endif

#ifdef RADCOOL_HOTHALO
static double effectivemass;
#endif

/*! \file forcetree_walk.c
 *  \brief Gravitational tree walk code
 *
 *  This file contains the various gravitational tree walks.
 */

/*! \brief variable for short-range lookup table
 *
 *  contains the factor needed for the short range
 *  contribution of the tree to the gravity force
 */


static float shortrange_table[NTAB + 1];

/*! \brief variable for short-range lookup table
 *
 *  contains the factor needed for the short range
 *  contribution of the tree to the potential energy
 */
static float shortrange_table_potential[NTAB + 1];

/*! \brief initializes the short range table
 *
 * The short range table contains the complementary error function
 * needed for the computation of the short range part of the gravity
 * force/potential in case of the TreePM algorithm.
 */
void force_short_range_init(void)
{
  for(int i = 0; i <= NTAB; i++)
    {
      double u = ((RCUT / 2.0) / NTAB) * i;

      shortrange_table_potential[i] = -erfc(u);   /* -r * g(r) */

      if(u > 0)
        shortrange_table[i] = (erfc(u) + 2.0 * u / sqrt(M_PI) * exp(-u * u) - 1.0) / (u * u);   /* -g'(r) - 1/r^2 */
      else
        shortrange_table[i] = 0;
    }
}

/*! \brief This routine calculates the (short range) force contribution
 *   for a given particle in case the Tree(PM) algorithm is used.
 *
 *  In the TreePM algorithm, the tree is walked only locally around the
 *  target coordinate.  Tree nodes that fall outside a box of half
 *  side-length Rcut= RCUT*ASMTH*MeshSize can be discarded. The short-range
 *  potential is modified by a complementary error function, multiplied
 *  with the Newtonian form. The resulting short-range suppression compared
 *  to the Newtonian force is tabulated, because looking up from this table
 *  is faster than recomputing the corresponding factor, despite the
 *  memory-access penalty (which reduces cache performance) incurred by the
 *  table.
 *
 *  Depending on the value of TypeOfOpeningCriterion, either the geometrical BH
 *  cell-opening criterion, or the `relative' opening criterion is used.
 *
 *  \param i index of the particle to be processed
 *  \param mode 0: process local particle (phase 1), 1: process imported particle (phase 2)
 *  \param thread_id: id of this thread
 *  \param measure_cost_flag whether the cost of the tree walk should be measured
 *  \return number of interactions processed for particle i
 */
int force_treeevaluate(gravdata_in * in, gravdata_out * out, int target, int mode, int thread_id, int numnodes, int *firstnode, int measure_cost_flag)
{
  struct NODE *nop = NULL;
#ifdef MULTIPLE_NODE_SOFTENING
  struct ExtNODE *extnop = 0;
#endif
#if (defined(PERIODIC) && !defined(GRAVITY_NOT_PERIODIC)) || defined(GRAVITY_TALLBOX)
  double xtmp, ytmp, ztmp;
#endif

  double acc_x = 0;
  double acc_y = 0;
  double acc_z = 0;
#ifdef EVALPOTENTIAL
  double pot = 0.0;
#endif


#ifdef TREE_RAD
  double Projection[NPIX];
  double ProjectionH2[NPIX];
  double ProjectionCO[NPIX];
  double pi, pisqrt, mass_rad, h2_mass_rad, co_mass_rad, area, angular_radius;
  double h2_mass_fraction, co_mass_fraction;
  double disti, ang_dist_try;
  double dx_rad, dy_rad, dz_rad;
  double dx_heal, dy_heal, dz_heal;
  double column, columnH2, columnCO;
  int i_projection_h2, i_projection_co;
  int j, particleid, inode, no_store;
#endif

#ifdef TREE_RAD_VEL
  double vx_p, vy_p, vz_p; /* velocities of particles */
  double vx_n, vy_n, vz_n; /* velocities of nodes */
  double v_diff2, v_th2, v_th2_compare; /* squared velocity difference between particle and node, thermal gas velocities for H - has to be rescaled for H2 and CO */
  double v_th2_compare_H2, v_th2_compare_CO;
  double f_overl; /* overlap fraction */
  v_th2 = 0;
#endif

#ifdef RADCOOL
  double youngstellarmass, oldstellarmass;
  double walkBirthTime, r2_ys, r2_os, r2_invfac_ys, r2_invfac_os;
  double dx_ys, dy_ys, dz_ys, dx_os, dy_os, dz_os;
  double walk_Phios, walk_Phins;
  walk_Phios = 0.0;
  walk_Phins = 0.0;
#ifdef RADCOOL_HOTHALO
  double T6gasmass, T7gasmass, T8gasmass;
  double log10temp, r2_T6, r2_T7, r2_T8, r2_invfac_T6, r2_invfac_T7, r2_invfac_T8;
  double dx_T6, dy_T6, dz_T6, dx_T7, dy_T7, dz_T7, dx_T8, dy_T8, dz_T8;
  double walk_PhiT6, walk_PhiT7, walk_PhiT8;
  walk_PhiT6 = 0.0;
  walk_PhiT7 = 0.0;
  walk_PhiT8 = 0.0;
#endif
#endif

#ifdef TREE_RAD
  pi = 3.14159265358979323846;
  pisqrt = 1.77245385090552;
  mass_rad = h2_mass_rad = co_mass_rad = 0.;

  for(j = 0; j < NPIX; j++)
    {
      Projection[j] = 0.0;
#ifdef TREE_RAD_H2
      ProjectionH2[j] = 0.0;
#endif
#ifdef TREE_RAD_CO
      ProjectionCO[j] = 0.0;
#endif
    }
#endif
#ifdef MODGRAV_EFF_MASS
  double modgrav_acc_x = 0;
  double modgrav_acc_y = 0;
  double modgrav_acc_z = 0;
#endif

  int ninteractions = 0;


  double pos_x = in->Pos[0];
  double pos_y = in->Pos[1];
  double pos_z = in->Pos[2];
  double aold = All.ErrTolForceAcc * in->OldAcc;
  double h_i = All.ForceSoftening[in->SofteningType];

#ifdef TREE_RAD_VEL
  vx_p  = in->Vel[0];
  vy_p  = in->Vel[1];
  vz_p  = in->Vel[2];
  v_th2 = in->Vth2;
#endif

#ifdef PMGRID
  double rcut = All.Rcut[0];
  double asmth = All.Asmth[0];
#ifdef PLACEHIGHRESREGION
  if(pmforce_is_particle_high_res(in->Type, in->Pos))
    {
      rcut = All.Rcut[1];
      asmth = All.Asmth[1];
    }
#endif

  double rcut2 = rcut * rcut;
  double asmthinv = 0.5 / asmth;
  double asmthinv2 = asmthinv * asmthinv;
  double asmthfac = asmthinv * (NTAB / (RCUT / 2.0));
#endif


  for(int k = 0; k < numnodes; k++)
    {
      int no;

      if(mode == 0)
        no = Tree_MaxPart;      /* root node */
      else
        {
          no = firstnode[k];
          no = Nodes[no].u.d.nextnode;  /* open it */
        }

      while(no >= 0)
        {
          double dx, dy, dz, r2, mass, hmax;

#ifdef MULTIPLE_NODE_SOFTENING
          int indi_flag1 = -1, indi_flag2 = 0;
#endif

          if(no < Tree_MaxPart) /* single particle */
            {
              dx = GRAVITY_NEAREST_X(Tree_Pos_list[3 * no + 0] - pos_x);
              dy = GRAVITY_NEAREST_Y(Tree_Pos_list[3 * no + 1] - pos_y);
              dz = GRAVITY_NEAREST_Z(Tree_Pos_list[3 * no + 2] - pos_z);
#ifdef TREE_RAD_VEL
	      vx_n = P[no].Vel[0];
	      vy_n = P[no].Vel[1];
	      vz_n = P[no].Vel[2];
#endif
              r2 = dx * dx + dy * dy + dz * dz;

              mass = P[no].Mass;

#ifdef RADCOOL
              if(P[no].Type == 4)
                {
                  r2_ys = r2;
                  r2_os = r2;
                  walkBirthTime = get_time_difference_in_Gyr(STP(np).BirthTime, All.Time);
                  if(walkBirthTime <= TIMEON_NEWSTARS)
                    {
                      youngstellarmass = P[no].Mass;
                      oldstellarmass = 0.0;
                    }
                  else if(walkBirthTime >= TIMEON_OLDSTARS)
                    {
                      oldstellarmass = P[no].Mass;
                      youngstellarmass = 0.0;
                    }
                  else
                    {
                      oldstellarmass = 0.0;
                      youngstellarmass = 0.0;
                    }
                }
              else
                {
                  oldstellarmass = 0.0;
                  youngstellarmass = 0.0;
                  r2_ys = r2;
                  r2_os = r2;
                }
#ifdef RADCOOL_HOTHALO
              if(P[no].Type == 0)
                {
                  r2_T6 = r2;
                  r2_T7 = r2;
                  r2_T8 = r2;
                  log10temp = log10(calculate_HH_temperature(SphP[no].Utherm
#ifdef COOLING
                                                             , SphP[no].Ne
#endif
                  ));
                  if((log10temp > TLOGMIN6) && (log10temp <= TLOGMAX6))
                    {
                      T6gasmass = P[no].Mass * SphP[no].Density
#ifdef RADCOOL_HOTHALO_METAl_BOOST
                          * (1.0 + T6METALBOOSTFACTOR * SphP[no].Metallicity / GFM_SOLAR_METALLICITY)
#endif
                          ;
                      T7gasmass = 0.0;
                      T8gasmass = 0.0;
                    }
                  else if((log10temp > TLOGMIN7) && (log10temp <= TLOGMAX7))
                    {
                      T7gasmass = P[no].Mass * SphP[no].Density
#ifdef RADCOOL_HOTHALO_METAl_BOOST
                          * (1.0 + T7METALBOOSTFACTOR * SphP[no].Metallicity / GFM_SOLAR_METALLICITY)
#endif
                          ;
                      T8gasmass = 0.0;
                      T6gasmass = 0.0;
                    }
                  else if((log10temp > TLOGMIN8) && (log10temp <= TLOGMAX8))
                    {
                      T8gasmass = P[no].Mass * SphP[no].Density
#ifdef RADCOOL_HOTHALO_METAl_BOOST
                          * (1.0 + T8METALBOOSTFACTOR * SphP[no].Metallicity / GFM_SOLAR_METALLICITY)
#endif
                          ;
                      T6gasmass = 0.0;
                      T7gasmass = 0.0;
                    }
                  else
                    {
                      T6gasmass = 0.0;
                      T7gasmass = 0.0;
                      T8gasmass = 0.0;
                    }
                }
              else
                {
                  T6gasmass = 0.0;
                  T7gasmass = 0.0;
                  T8gasmass = 0.0;
                  r2_T6 = r2;
                  r2_T7 = r2;
                  r2_T8 = r2;
                }
#endif
#endif

#ifdef TREE_RAD
              dx_rad = dx;
              dy_rad = dy;
              dz_rad = dz;
              inode = 0;
              no_store = no;
              if(P[no].Type == 0)
                {
                  mass_rad = P[no].Mass;
#ifdef TREE_RAD_H2
                  h2_mass_fraction = 2.0 * SphP[no].TracAbund[IH2] / (1.0 + 4.0 * ABHE);
                  h2_mass_rad = P[no].Mass * h2_mass_fraction;
#endif
#ifdef TREE_RAD_CO
                  co_mass_fraction = 28.0 * SphP[no].TracAbund[ICO] / (1.0 + 4.0 * ABHE);
                  co_mass_rad = P[no].Mass * co_mass_fraction;
#endif
                  angular_radius = 0.5 * pow(SphP[no].Volume, 1. / 3.); /* Will be converted to angle later */
                  area = 3.14159265358979 * angular_radius * angular_radius;
                }
              else
                {
                  mass_rad = 0.;
                  h2_mass_rad = 0.;
                  co_mass_rad = 0.;
                  area = 1.;
                  angular_radius = 0.;
                }
#endif
              if(measure_cost_flag)
                Thread[thread_id].P_CostCount[no]++;

              double h_j = All.ForceSoftening[P[no].SofteningType];

              hmax = (h_j > h_i) ? h_j : h_i;

              no = Nextnode[no];
            }
          else if(no < Tree_MaxPart + Tree_MaxNodes)    /* we have an  internal node */
            {
              if(mode == 1)
                {
                  if(no < Tree_FirstNonTopLevelNode)    /* we reached a top-level node again, which means that we are done with the branch */
                    {
                      no = -1;
                      continue;
                    }
                }

              nop = &Nodes[no];

#ifdef MODGRAV_EFF_MASS
              if(All.PM_Ti_endstep == All.Ti_Current && nop->mg_bitflag.amr_node == 1)
                {
                  //nop->GravCost += 1;
                  modgrav_tree_fifth_force(pos_x, pos_y, pos_z, no, asmthfac, h_i, shortrange_table, &modgrav_acc_x, &modgrav_acc_y, &modgrav_acc_z
#ifdef EVALPOTENTIAL
                              , &pot
#endif
                  );
                }
#endif

              mass = nop->u.d.mass;
              dx = GRAVITY_NEAREST_X(nop->u.d.s[0] - pos_x);
              dy = GRAVITY_NEAREST_Y(nop->u.d.s[1] - pos_y);
              dz = GRAVITY_NEAREST_Z(nop->u.d.s[2] - pos_z);

              r2 = dx * dx + dy * dy + dz * dz;

#if defined(PMGRID) && !defined(TREE_RAD)
              if(r2 > rcut2)
                {
                  /* check whether we can stop walking along this branch */
                  double eff_dist = rcut + 0.5 * nop->len;

                  double dist = GRAVITY_NEAREST_X(nop->center[0] - pos_x);
                  if(dist < -eff_dist || dist > eff_dist)
                    {
                      no = nop->u.d.sibling;
                      continue;
                    }

                  dist = GRAVITY_NEAREST_Y(nop->center[1] - pos_y);
                  if(dist < -eff_dist || dist > eff_dist)
                    {
                      no = nop->u.d.sibling;
                      continue;
                    }

                  dist = GRAVITY_NEAREST_Z(nop->center[2] - pos_z);
                  if(dist < -eff_dist || dist > eff_dist)
                    {
                      no = nop->u.d.sibling;
                      continue;
                    }
                }
#endif

#ifdef RADCOOL

              dx_ys = GRAVITY_NEAREST_X(nop->young_stellar_s[0] - pos_x);
              dy_ys = GRAVITY_NEAREST_X(nop->young_stellar_s[1] - pos_y);
              dz_ys = GRAVITY_NEAREST_X(nop->young_stellar_s[2] - pos_z);

              dx_os = GRAVITY_NEAREST_X(nop->old_stellar_s[0] - pos_x);
              dy_os = GRAVITY_NEAREST_X(nop->old_stellar_s[1] - pos_y);
              dz_os = GRAVITY_NEAREST_X(nop->old_stellar_s[2] - pos_z);

              r2_ys = dx_ys * dx_ys + dy_ys * dy_ys + dz_ys * dz_ys;
              r2_os = dx_os * dx_os + dy_os * dy_os + dz_os * dz_os;

              youngstellarmass = nop->young_stellar_mass;
              oldstellarmass = nop->old_stellar_mass;
#ifdef RADCOOL_HOTHALO
              dx_T6 = GRAVITY_NEAREST_X(nop->T6_gas_s[0] - pos_x);
              dy_T6 = GRAVITY_NEAREST_Y(nop->T6_gas_s[1] - pos_y);
              dz_T6 = GRAVITY_NEAREST_Z(nop->T6_gas_s[2] - pos_z);

              dx_T7 = GRAVITY_NEAREST_X(nop->T7_gas_s[0] - pos_x);
              dy_T7 = GRAVITY_NEAREST_Y(nop->T7_gas_s[1] - pos_y);
              dz_T7 = GRAVITY_NEAREST_Z(nop->T7_gas_s[2] - pos_z);

              dx_T8 = GRAVITY_NEAREST_X(nop->T8_gas_s[0] - pos_x);
              dy_T8 = GRAVITY_NEAREST_Y(nop->T8_gas_s[1] - pos_y);
              dz_T8 = GRAVITY_NEAREST_Z(nop->T8_gas_s[2] - pos_z);

              r2_T6 = dx_T6 * dx_T6 + dy_T6 * dy_T6 + dz_T6 * dz_T6;
              r2_T7 = dx_T7 * dx_T7 + dy_T7 * dy_T7 + dz_T7 * dz_T7;
              r2_T8 = dx_T8 * dx_T8 + dy_T8 * dy_T8 + dz_T8 * dz_T8;

              T6gasmass = nop->T6_gas_mass;
              T7gasmass = nop->T7_gas_mass;
              T8gasmass = nop->T8_gas_mass;
#endif
#endif

#ifdef TREE_RAD
              dx_rad = dx;
              dy_rad = dy;
              dz_rad = dz;
#ifdef TREE_RAD_VEL
              vx_n = nop->Vel[0];
              vy_n = nop->Vel[1];
              vz_n = nop->Vel[2];
#endif
              inode = 1;
              mass_rad = mass;
              if(isnan(mass_rad))
                printf("Node mass is NaN %d\n", no);
#ifdef TREE_RAD_H2
              h2_mass_rad = nop->u.d.h2mass;
#endif
#ifdef TREE_RAD_CO
              co_mass_rad = nop->u.d.comass;
#endif
              area = nop->len * nop->len;
              angular_radius = 0.5 * nop->len;
#endif


              if(All.ErrTolTheta)       /* check Barnes-Hut opening criterion */
                {
                  if(nop->len * nop->len > r2 * All.ErrTolTheta * All.ErrTolTheta)
                    {
                      /* open cell */
                      no = nop->u.d.nextnode;
                      continue;
                    }
                }
              else              /* check relative opening criterion */
                {
                  double len2 = nop->len * nop->len;

                  if(len2 > r2 * (1.2 * 1.2))  /* add a worst case protection */
                    {
                      /* open cell */
                      no = nop->u.d.nextnode;
                      continue;
                    }

#ifdef ACTIVATE_MINIMUM_OPENING_ANGLE
                  if(mass * len2 > r2 * r2 * aold && len2 > r2 * (0.4 * 0.4))
#else
                    if(mass * len2 > r2 * r2 * aold)
#endif
                      {
                        /* open cell */
                        no = nop->u.d.nextnode;
                        continue;
                      }

                  /* check in addition whether we lie inside the cell */

                  if(fabs(nop->center[0] - pos_x) < 0.60 * nop->len)
                    {
                      if(fabs(nop->center[1] - pos_y) < 0.60 * nop->len)
                        {
                          if(fabs(nop->center[2] - pos_z) < 0.60 * nop->len)
                            {
                              no = nop->u.d.nextnode;
                              continue;
                            }
                        }
                    }
                }

              double h_j = All.ForceSoftening[nop->u.d.maxsofttype];

              if(h_j > h_i)
                {
#ifdef MULTIPLE_NODE_SOFTENING
#ifdef ADAPTIVE_HYDRO_SOFTENING
                  if(nop->u.d.maxhydrosofttype != nop->u.d.minhydrosofttype)
                    if(ExtNodes[no].mass_per_type[0] > 0)
                      if(r2 < All.ForceSoftening[nop->u.d.maxhydrosofttype] * All.ForceSoftening[nop->u.d.maxhydrosofttype])
                        {
                          /* open cell */
                          no = nop->u.d.nextnode;
                          continue;
                        }
#endif
                  indi_flag1 = 0;
                  indi_flag2 = NSOFTTYPES;
#else
                  if(r2 < h_j * h_j)
                    {
                      /* open cell */
                      no = nop->u.d.nextnode;
                      continue;
                    }
#endif
                  hmax = h_j;
                }
              else
                hmax = h_i;

              /* ok, node can be used */
#ifdef MULTIPLE_NODE_SOFTENING
              extnop = &ExtNodes[no];
#endif
              if(measure_cost_flag && mass)
                Thread[thread_id].Node_CostCount[no]++;

#ifdef MODGRAV_EFF_MASS
              if(All.PM_Ti_endstep == All.Ti_Current && nop->mg_bitflag.amr_node == 0)
                {
                  //nop->GravCost += 1;
                  modgrav_tree_fifth_force(pos_x, pos_y, pos_z, no, asmthfac, h_i, shortrange_table,  &modgrav_acc_x, &modgrav_acc_y, &modgrav_acc_z
#ifdef EVALPOTENTIAL
                              , &pot
#endif
                  );
                }
#endif

              no = nop->u.d.sibling;
            }
          else if(no >= Tree_ImportedNodeOffset)        /* point from imported nodelist */
            {
              int n = no - Tree_ImportedNodeOffset;

              dx = GRAVITY_NEAREST_X(Tree_Points[n].Pos[0] - pos_x);
              dy = GRAVITY_NEAREST_Y(Tree_Points[n].Pos[1] - pos_y);
              dz = GRAVITY_NEAREST_Z(Tree_Points[n].Pos[2] - pos_z);

              r2 = dx * dx + dy * dy + dz * dz;

              mass = Tree_Points[n].Mass;

              if(measure_cost_flag)
                Thread[thread_id].TreePoints_CostCount[n]++;

#ifdef RADCOOL
              if((Tree_Points[n].Type) == 4)
                {
                  r2_ys = r2;
                  r2_os = r2;
                  walkBirthTime = get_time_difference_in_Gyr(Tree_Points[n].BirthTime, All.Time);
                  if(walkBirthTime <= TIMEON_NEWSTARS)
                    {
                      youngstellarmass = Tree_Points[n].Mass;
                      oldstellarmass = 0.0;
                    }
                  else if(walkBirthTime >= TIMEON_OLDSTARS)
                    {
                      oldstellarmass = Tree_Points[n].Mass;
                      youngstellarmass = 0.0;
                    }
                  else
                    {
                      youngstellarmass = 0.0;
                      oldstellarmass = 0.0;
                    }
                }
              else
                {
                  youngstellarmass = 0.0;
                  oldstellarmass = 0.0;
                  r2_ys = r2;
                  r2_os = r2;
                }
#ifdef RADCOOL_HOTHALO
              if((Tree_Points[n].Type) == 0)
                {
                  r2_T6 = r2;
                  r2_T7 = r2;
                  r2_T8 = r2;
                  log10temp = log10(calculate_HH_temperature(Tree_Points[n].Utherm
#ifdef COOLING
                                                             , Tree_Points[n].Ne
#endif
                  ));
                  if((log10temp > TLOGMIN6) && (log10temp <= TLOGMAX6))
                    {
                      T6gasmass = Tree_Points[n].Mass * Tree_Points[n].Density
#ifdef RADCOOL_HOTHALO_METAl_BOOST
                          * (1.0 + T6METALBOOSTFACTOR * Tree_Points[n].Metallicity / GFM_SOLAR_METALLICITY)
#endif
                          ;
                      T7gasmass = 0.0;
                      T8gasmass = 0.0;
                    }
                  else if((log10temp > TLOGMIN7) && (log10temp <= TLOGMAX7))
                    {
                      T7gasmass = Tree_Points[n].Mass * Tree_Points[n].Density
#ifdef RADCOOL_HOTHALO_METAl_BOOST
                          * (1.0 + T7METALBOOSTFACTOR * Tree_Points[n].Metallicity / GFM_SOLAR_METALLICITY)
#endif
                          ;
                      T8gasmass = 0.0;
                      T6gasmass = 0.0;
                    }
                  else if((log10temp > TLOGMIN8) && (log10temp <= TLOGMAX8))
                    {
                      T8gasmass = Tree_Points[n].Mass * Tree_Points[n].Density
#ifdef RADCOOL_HOTHALO_METAl_BOOST
                          * (1.0 + T8METALBOOSTFACTOR * Tree_Points[n].Metallicity / GFM_SOLAR_METALLICITY)
#endif
                          ;
                      T6gasmass = 0.0;
                      T7gasmass = 0.0;
                    }
                  else
                    {
                      T6gasmass = 0.0;
                      T7gasmass = 0.0;
                      T8gasmass = 0.0;
                    }
                }
              else
                {
                  T6gasmass = 0.0;
                  T7gasmass = 0.0;
                  T8gasmass = 0.0;
                  r2_T6 = r2;
                  r2_T7 = r2;
                  r2_T8 = r2;
                }
#endif
#endif

#ifdef TREE_RAD
              dx_rad = dx;
              dy_rad = dy;
              dz_rad = dz;
#ifdef TREE_RAD_VEL
              vx_n = Tree_Points[n].Vel[0];
              vy_n = Tree_Points[n].Vel[1];
              vz_n = Tree_Points[n].Vel[2];
#endif
              inode = 2;


              if((Tree_Points[n].Type) == 0)
                {
                  mass_rad = mass;
#ifdef TREE_RAD_H2
                  h2_mass_rad = mass * Tree_Points[n].TracAbund[IH2];
#endif
#ifdef TREE_RAD_CO
                  co_mass_rad = mass * Tree_Points[n].TracAbund[ICO];
#endif
                  area = 3.14159265358979 * Tree_Points[n].Hsml * Tree_Points[n].Hsml;
                  angular_radius = Tree_Points[n].Hsml;
                }
              else
                {
                  mass_rad = 0;
                  h2_mass_rad = 0;
                  co_mass_rad = 0;
                  area = 0;
                  angular_radius = 0;
                }
#endif

              double h_j = All.ForceSoftening[Tree_Points[n].SofteningType];

              hmax = (h_j > h_i) ? h_j : h_i;

              no = Nextnode[no - Tree_MaxNodes];
            }
          else                  /* pseudo particle */
            {
              if(mode == 0)
                {
                  tree_treefind_export_node_threads(no, target, thread_id);
                }

              no = Nextnode[no - Tree_MaxNodes];
              continue;
            }


          /* now calculate the shielding */
#ifdef TREE_RAD
#ifdef TREE_RAD_VEL
	  i_projection_h2 = i_projection_co = 1; /* Default is to project everything */

	  /* Compute squared velocity difference between particle and node. In comoving runs, we assume that the peculiar velocity dominates over the effects
	   * of the Hubble flow; this should generally be a good approximation for the gas that dominates the local shielding
	   */
	  v_diff2 = (vx_p-vx_n)*(vx_p-vx_n) + (vy_p-vy_n)*(vy_p-vy_n) + (vz_p-vz_n)*(vz_p-vz_n);
	  v_diff2 *= All.cf_a2inv; /* Convert to physical velocity; note that in non-comoving runs, All.cf_a2inv=1 */
	  f_overl = All.FracOverlap;
	  v_th2_compare = v_th2 * f_overl * f_overl;

	  /* Now need to rescale thermal velocity for H2 or CO */
	  v_th2_compare_H2 = v_th2_compare / 2.0;
	  v_th2_compare_CO = v_th2_compare / 28.0;
	  if(v_diff2 >  v_th2_compare_H2) /* Add matter only if relative velocities are small enough */
	    {
	      i_projection_h2 = 0; /* Do not project */
	      i_projection_co = 0; /* Thermal velocity of CO < thermal velocity of H2  */
	    }
	  else if (v_diff2 > v_th2_compare_CO)
	    {
	      i_projection_co = 0;
	    }
#else
	  /* Always project when not using TREE_RAD_VEL option */
	  i_projection_h2 = i_projection_co = 1;
#endif
          if(mass)
            {
#ifdef MULTIPLE_NODE_SOFTENING
              int type;
              for(type = indi_flag1; type < indi_flag2; type++)
                {
                  if(type >= 0)
                    {
                      mass = extnop->mass_per_type[type];
                    }

                  if(mass)
                    {
#endif
                      /*compute column densities */
#ifdef MULTIPLE_NODE_SOFTENING
                      if(type == indi_flag1)
#endif
                        {
                          // Do TreeCol projection
                          pi = 3.14159265358979323846;
                          if(mass_rad > 0. || h2_mass_rad > 0.)
                            {
                              if(area > 0)
                                {
                                  double dist = sqrt(dx_rad * dx_rad + dy_rad * dy_rad + dz_rad * dz_rad + 1e-99);

				  double ascale = 1.0;
				  if (All.ComovingIntegrationOn) {
				    ascale = All.Time;
				  }

                                  if(dist < All.TreeColMaxDistance * All.HubbleParam / ascale)
                                    {
                                      /* Don't include contributions further away than our specified maximum distance. This is useful for e.g.
                                       * galactic simulations, where we expect the typical separation between UV sources to be much less than
                                       * the size of the simulation volume.
                                       */
                                      disti = 1.0 / dist;
				      /* Note that it doesn't matter whether the distances here are comoving or real, as we're only dealing with the ratios. */
                                      dx_heal = dx_rad * disti;
                                      dy_heal = dy_rad * disti;
                                      dz_heal = dz_rad * disti;
                                      ang_dist_try = angular_radius * disti; /* "angular_radius" is still a physical distance here, but now gets converted */
                                      /* Use small angle formula only when angles are SMALL! */
                                      if(ang_dist_try < 0.1)
                                        angular_radius = ang_dist_try;
                                      else
                                        angular_radius = atan2(angular_radius, dist);

				      /* Note: when using comoving coordinates, these are the comoving column densities. They are converted to physical units
				       * inside sgchem.c.
				       */
                                      column = mass_rad / area;
#ifdef TREE_RAD_H2
				      if (i_projection_h2)
					{
                                          columnH2 = h2_mass_rad / area;
					}
				      else
					{
					  columnH2 = 0.0;
					}
#endif
#ifdef TREE_RAD_CO
				      if (i_projection_co)
					{
                                          columnCO = co_mass_rad / area;
					}
				      else
					{
					  columnCO = 0.0;
					}
#endif
                                      if(column > 0 || columnH2 > 0 || columnCO > 0)
                                        {
                                          PROJECT_COLUMN(&column, &columnH2, &columnCO, &dx_heal, &dy_heal, &dz_heal, &angular_radius, Projection, ProjectionH2, ProjectionCO, &particleid);
                                        }
                                    }
                                }

                            }
                        }
#ifdef MULTIPLE_NODE_SOFTENING
                    }
                }
#endif
            }
#endif

          /* now evaluate the multipole moment */
          if(mass)
            {
              double r = sqrt(r2);

#ifdef PMGRID
              double tabentry = asmthfac * r;
              int tabindex = (int) tabentry;

              if(tabindex < NTAB)
                {
                  double tabweight = tabentry - tabindex;
                  double factor_force = (1.0 - tabweight) * shortrange_table[tabindex] + tabweight * shortrange_table[tabindex + 1];
#ifdef EVALPOTENTIAL
                  double factor_pot = (1.0 - tabweight) * shortrange_table_potential[tabindex] + tabweight * shortrange_table_potential[tabindex + 1];
#endif
#endif

#ifdef MULTIPLE_NODE_SOFTENING
                  for(int type = indi_flag1; type < indi_flag2; type++)
                    {
                      if(type >= 0)
                        {
                          mass = extnop->mass_per_type[type];
                          double h_j;
#ifdef ADAPTIVE_HYDRO_SOFTENING
                          if(type == 0)
                            h_j = All.ForceSoftening[nop->u.d.maxhydrosofttype];
                          else
#endif
                            h_j = All.ForceSoftening[type];

                          hmax = (h_j > h_i) ? h_j : h_i;
                        }

                      if(mass)
                        {
#endif
                          double fac;
#ifdef EVALPOTENTIAL
                          double wp;
#endif

                          if(r >= hmax)
                            {
                              double rinv = 1.0 / r;
                              double rinv3 = rinv * rinv * rinv;
#ifdef PMGRID
                              fac =  rinv3 + rinv * factor_force * asmthinv2; /* fac  = -g'(r)/r */
#ifdef EVALPOTENTIAL
                              wp = rinv * factor_pot; /* wp   = -g(r)    */
#endif
#else
                              fac = rinv3;
#ifdef EVALPOTENTIAL
                              wp = -rinv;
#endif
#endif
                            }
                          else
                            {
                              double h_inv = 1.0 / hmax;
                              double h3_inv = h_inv * h_inv * h_inv;
                              double u = r * h_inv;

                              if(u < 0.5)
                                {
                                  double u2 = u * u;
                                  fac = h3_inv * (SOFTFAC1 + u2 * (SOFTFAC2 * u + SOFTFAC3));
#ifdef EVALPOTENTIAL
                                  wp = h_inv * (SOFTFAC4 + u2 * (SOFTFAC5 + u2 * (SOFTFAC6 * u + SOFTFAC7)));
#endif
                                }
                              else
                                {
                                  double u2 = u * u;
                                  double u3 = u2 * u;
                                  fac = h3_inv * (SOFTFAC8 + SOFTFAC9 * u + SOFTFAC10 * u2 + SOFTFAC11 * u3 + SOFTFAC12 / u3);
#ifdef EVALPOTENTIAL
                                  wp = h_inv * (SOFTFAC13 + SOFTFAC14 / u + u2 * (SOFTFAC1 + u * (SOFTFAC15 + u * (SOFTFAC16 + SOFTFAC17 * u))));
#endif
                                }

#ifdef PMGRID
                              if(r > 0)
                                {
                                  double rinv = 1.0 / r;
                                  fac += rinv * factor_force * asmthinv2; /* fac  = -g'(r)/r */
#ifdef EVALPOTENTIAL
                                  wp += rinv * (factor_pot + 1.0);     /* wp   = -g(r)    */
#endif
                                }
#endif
                            }


#ifdef RADCOOL
#ifdef MULTIPLE_NODE_SOFTENING
                          if(type == indi_flag1)
#endif
                            {
                              r2_invfac_ys = 1.0 / (r2_ys + hmax * hmax);
                              r2_invfac_os = 1.0 / (r2_os + hmax * hmax);
                              walk_Phins += youngstellarmass * r2_invfac_ys;
                              walk_Phios += oldstellarmass * r2_invfac_os;
#ifdef RADCOOL_HOTHALO
                              r2_invfac_T6 = 1.0 / (r2_T6 + hmax * hmax);
                              r2_invfac_T7 = 1.0 / (r2_T7 + hmax * hmax);
                              r2_invfac_T8 = 1.0 / (r2_T8 + hmax * hmax);
                              walk_PhiT6 += T6gasmass * r2_invfac_T6;
                              walk_PhiT7 += T6gasmass * r2_invfac_T7;
                              walk_PhiT8 += T6gasmass * r2_invfac_T8;
#endif
                            }
#endif

#ifdef EVALPOTENTIAL
                          pot += mass * wp;
#endif
                          fac *= mass;

                          acc_x += dx * fac;
                          acc_y += dy * fac;
                          acc_z += dz * fac;


#if !defined(PMGRID) && defined(SELFGRAVITY) && defined(PERIODIC) && !defined(GRAVITY_NOT_PERIODIC) && !defined(ONEDIMS_SPHERICAL)
                          double fcorr[3];
                          ewald_corr(dx, dy, dz, fcorr);
                          acc_x += mass * fcorr[0];
                          acc_y += mass * fcorr[1];
                          acc_z += mass * fcorr[2];
#ifdef EVALPOTENTIAL
                          pot += mass * ewald_pot_corr(dx, dy, dz);
#endif
#endif

#ifdef MULTIPLE_NODE_SOFTENING
                        }
                    }
#endif
                  ninteractions++;
#ifdef PMGRID
                }
#endif
            }
        }
    }


  out->Acc[0] = acc_x;
  out->Acc[1] = acc_y;
  out->Acc[2] = acc_z;
#ifdef EVALPOTENTIAL
  out->Potential = pot;
#endif
#ifdef NO_GRAVITY_TYPE
  if (in->Type == NO_GRAVITY_TYPE)
    {
      out->Acc[0] = 0.0;
      out->Acc[1] = 0.0;
      out->Acc[2] = 0.0;
#ifdef EVALPOTENTIAL
      out->Potential = 0.0;
#endif
    }
#endif
#ifdef OUTPUTGRAVINTERACTIONS
  out->GravInteractions = ninteractions;
#endif
#ifdef RADCOOL
  out->Phios = walk_Phios;
  out->Phins = walk_Phins;
#ifdef RADCOOL_HOTHALO
  out->PhiT6 = walk_PhiT6;
  out->PhiT7 = walk_PhiT7;
  out->PhiT8 = walk_PhiT8;
#endif
#endif
#ifdef TREE_RAD
  for(j = 0; j < NPIX; j++)
    {
      out->Projection[j] = Projection[j];
#ifdef TREE_RAD_H2
      out->ProjectionH2[j] = ProjectionH2[j];
#endif
#ifdef TREE_RAD_CO
      out->ProjectionCO[j] = ProjectionCO[j];
#endif
    }
#endif
#ifdef MODGRAV
  out->ModgravAcc[0] = modgrav_acc_x;
  out->ModgravAcc[1] = modgrav_acc_y;
  out->ModgravAcc[2] = modgrav_acc_z;
#endif

  return ninteractions;
}




int tree_treefind_export_node_threads(int no, int i, int thread_id)
{
  /* The task indicated by the pseudoparticle node */
  int task = DomainNewTask[no - (Tree_MaxPart + Tree_MaxNodes)];

  if(Thread[thread_id].Exportflag[task] != i)
    {
      Thread[thread_id].Exportflag[task] = i;
      int nexp = Thread[thread_id].Nexport++;
      Thread[thread_id].PartList[nexp].Task = task;
      Thread[thread_id].PartList[nexp].Index = i;
      Thread[thread_id].ExportSpace -= Thread[thread_id].ItemSize;
    }

  int nexp = Thread[thread_id].NexportNodes++;
  nexp = -1 - nexp;
  struct datanodelist *nodelist = (struct datanodelist *) (((char *) Thread[thread_id].PartList) + Thread[thread_id].InitialSpace);
  nodelist[nexp].Task = task;
  nodelist[nexp].Index = i;
  nodelist[nexp].Node = DomainNodeIndex[no - (Tree_MaxPart + Tree_MaxNodes)];
  Thread[thread_id].ExportSpace -= sizeof(struct datanodelist) + sizeof(int);
  return 0;
}




#ifdef ALLOW_DIRECT_SUMMATION
void force_evaluate_direct(int target, int result_idx, int nimport)
{
#if defined(PERIODIC) && !defined(GRAVITY_NOT_PERIODIC)
  double xtmp, ytmp, ztmp;
#endif

  double acc_x = 0;
  double acc_y = 0;
  double acc_z = 0;
#ifdef EVALPOTENTIAL
  double pot = 0.0;
#endif

  double pos_x = DirectDataAll[target].Pos[0];
  double pos_y = DirectDataAll[target].Pos[1];
  double pos_z = DirectDataAll[target].Pos[2];
  double h_i = All.ForceSoftening[DirectDataAll[target].SofteningType];

#ifdef PMGRID
  double asmth = All.Asmth[0];
#if defined(PLACEHIGHRESREGION)
  int ptype_i = DirectDataAll[target].Type;
  if(pmforce_is_particle_high_res(ptype_i, DirectDataAll[target].Pos))
    asmth = All.Asmth[1];
#endif
  double asmthinv = 0.5 / asmth;
  double asmthinv2 = asmthinv * asmthinv;
  double asmthfac = asmthinv * (NTAB / (RCUT / 2.0));
#endif

  for(int j = 0; j < nimport; j++)
    {
      double h_j = All.ForceSoftening[DirectDataAll[j].SofteningType];

      double hmax = (h_j > h_i) ? h_j : h_i;

      double dx = GRAVITY_NEAREST_X(DirectDataAll[j].Pos[0] - pos_x);
      double dy = GRAVITY_NEAREST_Y(DirectDataAll[j].Pos[1] - pos_y);
      double dz = GRAVITY_NEAREST_Z(DirectDataAll[j].Pos[2] - pos_z);

      double r2 = dx * dx + dy * dy + dz * dz;

      double mass = DirectDataAll[j].Mass;

      /* now evaluate the force component */

      double r = sqrt(r2);

#ifdef PMGRID
      double tabentry = asmthfac * r;
      int tabindex = (int) tabentry;

      if(tabindex < NTAB)
        {
          double tabweight = tabentry - tabindex;
          double factor_force = (1.0 - tabweight) * shortrange_table[tabindex] + tabweight * shortrange_table[tabindex + 1];
#ifdef EVALPOTENTIAL
          double factor_pot = (1.0 - tabweight) * shortrange_table_potential[tabindex] + tabweight * shortrange_table_potential[tabindex + 1];
#endif
#endif

          double fac;
#ifdef EVALPOTENTIAL
          double wp;
#endif

          if(r >= hmax)
            {
              double rinv = 1.0 / r;
              double rinv3 = rinv * rinv * rinv;
#ifdef PMGRID
              fac =  rinv3 + rinv * factor_force * asmthinv2;   /* fac  = -g'(r)/r */
#ifdef EVALPOTENTIAL
              wp = rinv * factor_pot; /* wp   = -g(r)    */
#endif
#else
              fac = rinv3;
#ifdef EVALPOTENTIAL
              wp = -rinv;
#endif
#endif
            }
          else
            {
              double h_inv = 1.0 / hmax;
              double h3_inv = h_inv * h_inv * h_inv;
              double u = r * h_inv;

              if(u < 0.5)
                {
                  double u2 = u * u;
                  fac = h3_inv * (SOFTFAC1 + u2 * (SOFTFAC2 * u + SOFTFAC3));
#ifdef EVALPOTENTIAL
                  wp = h_inv * (SOFTFAC4 + u2 * (SOFTFAC5 + u2 * (SOFTFAC6 * u + SOFTFAC7)));
#endif
                }
              else
                {
                  double u2 = u * u;
                  double u3 = u2 * u;
                  fac = h3_inv * (SOFTFAC8 + SOFTFAC9 * u + SOFTFAC10 * u2 + SOFTFAC11 * u3 + SOFTFAC12 / u3);
#ifdef EVALPOTENTIAL
                  wp = h_inv * (SOFTFAC13 + SOFTFAC14 / u + u2 * (SOFTFAC1 + u * (SOFTFAC15 + u * (SOFTFAC16 + SOFTFAC17 * u))));
#endif
                }
#ifdef PMGRID
              if(r > 0)
                {
                  double rinv = 1.0 / r;
                  fac += rinv * factor_force * asmthinv2; /* fac  = -g'(r)/r */
#ifdef EVALPOTENTIAL
                  wp += rinv * (factor_pot + 1.0); /* wp   = -g(r)    */
#endif
                }
#endif
            }

#ifdef EVALPOTENTIAL
          pot += mass * wp;
#endif
          fac *= mass;

          acc_x += dx * fac;
          acc_y += dy * fac;
          acc_z += dz * fac;


#if !defined(PMGRID) && defined(SELFGRAVITY) && !defined(GRAVITY_NOT_PERIODIC) && defined(PERIODIC) && !defined(ONEDIMS_SPHERICAL)
          {
            double fcorr[3];
            ewald_corr(dx, dy, dz, fcorr);
            acc_x += mass * fcorr[0];
            acc_y += mass * fcorr[1];
            acc_z += mass * fcorr[2];
#if defined(EVALPOTENTIAL)
            pot += mass * ewald_pot_corr(dx, dy, dz);
#endif
          }
#endif


#ifdef PMGRID
        }
#endif

    }

  DirectAccOut[result_idx].Acc[0] = acc_x;
  DirectAccOut[result_idx].Acc[1] = acc_y;
  DirectAccOut[result_idx].Acc[2] = acc_z;
#ifdef EVALPOTENTIAL
  DirectAccOut[result_idx].Potential = pot;
#endif
}
#endif
