/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/criterion_refinement.c
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

#if defined(REFINEMENT_SPLIT_CELLS) && !defined(ONEDIMS)

static int refine_criterion_anna(int i);
static int refine_criterion_gmcs(int i);
static int anna_refinement_criteria(int i);
static int jeans_refinement_criteria(int i);

#ifdef ROTATING_HIGHRES_REGION
static int refine_rotating_highres_region(int i);
#endif

#if defined(SPIRAL) && defined(RAMP_REFINE)
static int refine_disczoom(int i);
#endif

#ifdef REFINEMENT_MERGE_CELLS
char *FlagDoNotRefine;
#endif


int should_this_cell_be_split(int i)
{
#ifdef REFINEMENT_MERGE_CELLS
  if(FlagDoNotRefine[i])
    return 0;
#endif

  if(P[i].Mass == 0 && P[i].ID == 0)    /* skip cells that have been swallowed or dissolved */
    return 0;

#ifdef BOUNDARY_INFLOWOUTFLOW_MINID
  if(P[i].ID >= BOUNDARY_INFLOWOUTFLOW_MINID && P[i].ID < BOUNDARY_INFLOWOUTFLOW_MAXID)
    return 0;
#endif

#if defined(BOUNDARY_STICKY_MINID) && defined(BOUNDARY_STICKY_MAXID)
  if(P[i].ID >= BOUNDARY_STICKY_MINID && P[i].ID < BOUNDARY_STICKY_MAXID)
    return 0;
#endif

#ifdef STICKYFLAGS
  if(SphP[i].StickyFlag > 0)
    return 0;
#endif

#if defined(BOUNDARY_REFL_FLUIDSIDE_MINID) && defined(BOUNDARY_REFL_FLUIDSIDE_MAXID)
  if(P[i].ID >= BOUNDARY_REFL_FLUIDSIDE_MINID && P[i].ID < BOUNDARY_REFL_FLUIDSIDE_MAXID)
    return 0;
#endif

#if defined(BOUNDARY_REFL_SOLIDSIDE_MINID) && defined(BOUNDARY_REFL_SOLIDSIDE_MAXID)
  if(P[i].ID >= BOUNDARY_REFL_SOLIDSIDE_MINID && P[i].ID < BOUNDARY_REFL_SOLIDSIDE_MAXID)
    return 0;
#endif

#ifdef SPECIAL_BOUNDARY
  if(P[i].ID == -1 || P[i].ID == -2)
    return 0;
#endif

#if defined(REFINEMENT_VOLUME_LIMIT) && !defined(WINDTUNNEL_REFINEMENT_VOLUME_LIMIT)
  double maxvolume = All.MaxVolume;
  double minvolume = All.MinVolume;
#ifdef REFINE_ONLY_WITH_TRACER
  if(P[i].TracerHead != -1)
    {
      maxvolume = All.MaxTracerVolume;
      minvolume = All.MinTracerVolume;
    }
#endif

  if(SphP[i].Volume > 2. * maxvolume)
    if(can_this_cell_be_split(i))
      return 1;
  
  if(SphP[i].Volume < 2. * minvolume)
    return 0;

#ifdef REFINEMENT_VOLUME_LIMIT_MASS_LIMIT
  if(P[i].Mass < 0.1 * All.TargetGasMass)
    return 0;
#endif

  if(refine_criterion_volume(i))
    if(can_this_cell_be_split(i))
      	return 1;
#endif

#ifdef REFINEMENT_AROUND_BH
  if(SphP[i].RefBHFlag)
    {
      if(SphP[i].RefBHMaxRad < get_cell_radius(i))
        {
          if(P[i].Mass < All.RefBHMinCellMass)
            {
              return 0;
            }
#if (REFINEMENT_AROUND_BH==0)
          if(can_this_cell_be_split(i))
            return 1;
#endif
#if (REFINEMENT_AROUND_BH==1)
          return 1;
#endif
        }
    }
#endif

#ifdef REFINEMENT_AROUND_DM
  int j;
  for(j = 0; j < All.TotPartDM; j++)
    {
      double dx = P[i].Pos[0] - DMPartListGlobal[j].pos[0];
      double dy = P[i].Pos[1] - DMPartListGlobal[j].pos[1];
      double dz = P[i].Pos[2] - DMPartListGlobal[j].pos[2];
      double dist = sqrt(dx * dx + dy * dy + dz * dz);
      if(dist < 1.0 * DMPartListGlobal[j].softening)
        if(get_cell_radius(i) > DMPartListGlobal[j].softening / All.RefinementCellsPerSoftening && can_this_cell_be_split(i))
          {
#ifdef DEBUG_REFINE
            printf("Refining particle ID %d with radius %e around DM particle %d with softening length %e\n", P[i].ID, get_cell_radius(i), DMPartListGlobal[j].ID, DMPartListGlobal[j].softening);
#endif
            return 1;
          }
    }
#endif


#if defined(CIRCUMSTELLAR) && defined(CIRCUMSTELLAR_REFINEMENTS)
  double StellarDistance;
  StellarDistance = get_circumstellar_distance(i);
  if((StellarDistance < All.CircumstellarDerefinementDistance) && (SphP[i].Volume < 0.75 / M_PI / pow(All.CircumstellarDerefinementDistance, 3)))
    return 0;
#endif

#if defined SINK_PARTICLES_REFINEMENT_LIMIT
  int isink;
#ifdef JEANS_REFINEMENT
  double ncells_per_sink = JEANS_REFINEMENT;
#else
  double ncells_per_sink = 4;
#endif
  double dx, dy, dz, dist, dist_between_sink_cell_surface;
  double cell_rad = get_cell_radius(i);
  int number_of_sink_radii;
  if(can_this_cell_be_split(i))
    {
      if (NSinksAllTasks > 0)
        {
          for(isink = 0; isink < NSinksAllTasks; isink++)
            {
              /* Ensure that the regions around the sink particles are well resolved, and that the
                 cells refine smoothly around the sinks. 
              */
              dx = P[i].Pos[0] - SinkP[isink].Pos[0];
              dy = P[i].Pos[1] - SinkP[isink].Pos[1];
              dz = P[i].Pos[2] - SinkP[isink].Pos[2];	
              dist = sqrt(dx * dx  +  dy * dy  +  dz * dz);
              dist_between_sink_cell_surface = dist - cell_rad - SinkAccretionRadius;
              if (dist_between_sink_cell_surface < 0) 
                dist_between_sink_cell_surface = 0;
              number_of_sink_radii = dist_between_sink_cell_surface / SinkAccretionRadius;
              if ((number_of_sink_radii + 1)*dist_between_sink_cell_surface/ncells_per_sink > cell_rad)
                return 1;
            }
        }
    }
#endif

  switch (All.RefinementCriterion)      /* select the function that evaluates the refinement criterion */
    {
    case 0:
      return 0;
      break;

    case 1:
      return refine_criterion_default(i);
      break;

    case 2:
      return refine_criterion_jeans_ref(i);
      break;

    case 3:
      return refine_criterion_windtunnel(i);
      break;

    case 4:
      return refine_criterion_special_boundary(i);
      break;

    case 5:
      return refine_criterion_gmcs(i);
      break;

#ifdef RAMP_REFINE
    case 6:
      return refine_disczoom(i);
      break;
#endif

#ifdef SGCHEM
    case 7:
      return refine_criterion_anna(i);
      break;
#endif

#ifdef DISC_REFINE_ONLY
    case 8:
      return refine_disconly(i);
      break;
#endif

#ifdef ROTATING_HIGHRES_REGION
    case 9:
      return refine_rotating_highres_region(i);
#endif

    default:
      terminate("invalid refinement criterion specified");
      break;
    }

  return 0;
}


int can_this_cell_be_split(int i)
{
#ifdef REGULARIZE_MESH_FACE_ANGLE
  if(SphP[i].MaxFaceAngle < 1.5 * All.CellMaxAngleFactor)
#else
  double dx = nearest_x(P[i].Pos[0] - SphP[i].Center[0]);
  double dy = nearest_y(P[i].Pos[1] - SphP[i].Center[1]);
  double dz = nearest_z(P[i].Pos[2] - SphP[i].Center[2]);
  double d = sqrt(dx * dx + dy * dy + dz * dz);
  double cellrad = get_cell_radius(i);

  if(d < 2.0 * All.CellShapingFactor * cellrad) /* only refine cells which are reasonably 'round' */
#endif
    return 1;

  return 0;
}


int refine_criterion_default(int i)
{
#ifdef REFINEMENT_HIGH_RES_GAS
  if(SphP[i].AllowRefinement != 0)
#ifndef TGSET
    if(SphP[i].HighResMass > HIGHRESMASSFAC * P[i].Mass)
#endif
#endif
      if(can_this_cell_be_split(i) && P[i].Mass > 2.0 * All.TargetGasMass)
	return 1;
  return 0;                     /* default is not to refine */
}


int refine_criterion_jeans_ref(int i)
{
#ifdef REFINEMENT_HIGH_RES_GAS
  if(SphP[i].AllowRefinement != 0)
#ifndef TGSET
    if(SphP[i].HighResMass > HIGHRESMASSFAC * P[i].Mass)
#endif
#endif
      if(can_this_cell_be_split(i))
        {
          if(P[i].Mass > 2.0 * All.TargetGasMass)
            return 1;
#ifdef TGSET
          return tgset_jeans_ref(0, i);
#else
#ifdef JEANS_REFINEMENT
          return jeans_refinement_criteria(i);
#else
          return 0;
#endif
#endif
        }

  return 0;
}

int refine_criterion_gmcs(int i)
{
  if(can_this_cell_be_split(i))
    {
#ifdef GMC_REFINEMENT
      return gmc_refinement_criteria(i);        // currently only jeans criterion here */
#else
      return 0;
#endif
    }
  return 0;
}

int refine_criterion_anna(int i)
{
  if(can_this_cell_be_split(i))
    {
#ifdef SGCHEM
      return anna_refinement_criteria(i);        // jeans criterion here, applied to high densities */
#else
      return 0;
#endif
    }
  return 0;
}

int jeans_refinement_criteria(int i)
{
#ifdef SINK_PARTICLES
  if(SphP[i].Density > 0.1 * SinkCreationDensityCodeUnits)
    return 0;
#endif
#ifdef JEANS_REFINEMENT
  if(can_this_cell_be_split(i))
    {
      double jeans_number, jeans_length, sound_speed, dx;
      sound_speed = sqrt(GAMMA * SphP[i].Pressure / SphP[i].Density);
      //sound_speed = get_sound_speed(i);
      jeans_length = sqrt(M_PI / All.G / SphP[i].Density) * sound_speed;
      dx = 2.0 * get_cell_radius(i);
      jeans_number = jeans_length / dx;
#ifdef EXTERNALSHEARBOX
      double bscale0 = 61.0 * (All.ShearBoxFg / 0.1 / All.ShearBoxMu) / (All.ShearBoxSigma0 / 10.0);
      double rho_thresh = All.ShearBoxSigma0 / 2.0 / bscale0;
      if(SphP[i].Density < rho_thresh)
        return 0;
#endif
      if(jeans_number < JEANS_REFINEMENT)
      {
        return 1;
      }
    }
#endif
  return 0;
}

int anna_refinement_criteria(int i)
{
  if(can_this_cell_be_split(i))
  {
    double n_dens_anna;
    double n_dens_anna_min;

    n_dens_anna = SphP[i].Density * All.UnitDensity_in_cgs * All.cf_a3inv / (All.HubbleParam *All.HubbleParam);
    n_dens_anna = n_dens_anna * 0.81967 / PROTONMASS;
    n_dens_anna_min = 1.0e2;    /* for number densities larger than 100, use jeans refinement! */
#ifdef JEANS_REFINEMENT
    if (n_dens_anna > n_dens_anna_min)
    {
      return jeans_refinement_criteria(i);
    }  
#endif
    return refine_criterion_default(i);
  }
  return 0;
}



int refine_criterion_windtunnel(int i)
{
#ifdef WINDTUNNEL
  double boxsizes[3] = { boxSize_X, boxSize_Y, boxSize_Z };
  if(P[i].Type == 0)
    {
		
#ifdef WINDTUNNEL_REFINEMENT_VOLUME_LIMIT

        if(SphP[i].Volume > 2. * All.MaxVolume)
            if(can_this_cell_be_split(i))
                return 1;

        if(SphP[i].Volume < 2. * All.MinVolume)
            return 0;

        if(P[i].Pos[WINDTUNNEL_COORD] < All.InjectionRegion)//injection region
        { 
            //if(SphP[i].Volume > 1.5 * All.InjectionVolume)//|| P[i].Mass > 2.0 * All.TargetGasMass
            //    if(can_this_cell_be_split(i))
            //        return 1;
        }   
        else if(P[i].Pos[WINDTUNNEL_COORD] > All.InjectionRegion && boxsizes[WINDTUNNEL_COORD] - P[i].Pos[WINDTUNNEL_COORD] > All.InjectionRegion )//the "normal region"
        {

            if(refine_criterion_volume(i))
                if(can_this_cell_be_split(i))
                    return 1;  
        }
        else if (boxsizes[WINDTUNNEL_COORD] - P[i].Pos[WINDTUNNEL_COORD] < All.InjectionRegion)//outflow region
        {
            //if (SphP[i].Volume > 2.0 * All.InjectionVolume)
            //    if(can_this_cell_be_split(i))
            //        return 1;		   
        }
#endif

       if(can_this_cell_be_split(i))  
        {
          if((P[i].Pos[WINDTUNNEL_COORD] < All.InjectionRegion && SphP[i].Volume > 1.5 * All.InjectionVolume) || (SphP[i].Volume > 2.0 * All.InjectionVolume))
            return 1;

          if(P[i].Mass > 2.0 * All.TargetGasMass)
            return 1;
        }        
    }

#endif

  return 0;                     /* default is not to refine */
}


int refine_criterion_special_boundary(int i)
{
#ifdef SPECIAL_BOUNDARY
  if(SphP[i].MinDistBoundaryCell < 4 * All.BoundaryLayerScaleFactor && P[i].ID != -1 && P[i].ID != -2)
    return 0;
  //else
  //  return refine_criterion_default(i);
#endif
  return 0;
}


#ifdef REFINEMENT_VOLUME_LIMIT
int refine_criterion_volume(int i)
{
  if(All.MaxVolumeDiff > 0 && SphP[i].Volume > All.MaxVolumeDiff * SphP[i].MinNgbVolume)
    {
#ifdef REGULARIZE_MESH_FACE_ANGLE
      if(SphP[i].MaxFaceAngle < 1.5 * All.CellMaxAngleFactor)
        return 1;
#else
      double dx = nearest_x(P[i].Pos[0] - SphP[i].Center[0]);
      double dy = nearest_y(P[i].Pos[1] - SphP[i].Center[1]);
      double dz = nearest_z(P[i].Pos[2] - SphP[i].Center[2]);

      double d = sqrt(dx * dx + dy * dy + dz * dz);

      double cellrad = get_cell_radius(i);

      if(d < 2.0 * All.CellShapingFactor * cellrad)     /* only refine cells which are reasonably 'round' */
        return 1;
#endif
    }

  return 0;
}
#endif

#ifdef ROTATING_HIGHRES_REGION
int refine_rotating_highres_region(int i)
{
  double dx = P[i].Pos[0] - boxHalf_X; 
  double dy = P[i].Pos[1] - boxHalf_Y; 
  double dz = P[i].Pos[2] - boxHalf_Z; 
  double erre2 = pow(dx,2) + pow(dy,2) + pow(dz,2);

  double zoomtarget;
  if(erre2 < pow(20. * KILOPARSEC / All.UnitLength_in_cm,2))
    zoomtarget = All.TargetGasMass;
  else if(erre2 < pow(50. * KILOPARSEC / All.UnitLength_in_cm,2))
    zoomtarget = 10. * All.TargetGasMass;
  else
    zoomtarget = 100. * All.TargetGasMass;

  if (fabs(dz) < All.Highres_deltaz)
    {
      double Rref = sqrt(pow((All.Highres_y0 - boxHalf_Y),2) + pow((All.Highres_x0 - boxHalf_X),2)); //radius at which the high res region is centered

      double dR = sqrt(dx*dx+dy*dy);
  
      if (fabs(Rref - dR) < All.Highres_delta)
        {
          double theta_0 = atan2((All.Highres_y0 - boxHalf_Y), (All.Highres_x0 - boxHalf_X)); //initial angle of the high res region
          double wr = All.Highres_vrot / Rref; //angular velocity of the high res region
  
          double thref = theta_0 + wr * (All.Time - All.Highres_t0); //current angle of the high res region
          thref -= (thref < 0 ? -1. : 1.) * floor( (fabs(thref) + M_PI) / (2.*M_PI) ) * 2.*M_PI; //so that thref is always between -pi and pi
  
          double dtheta = atan2(dy, dx);
  
          double delta_th = dtheta - thref;
  
          if (delta_th > M_PI)
            delta_th = 2.*M_PI - delta_th;
          else if (delta_th < -M_PI)
            delta_th = 2.*M_PI + delta_th;
  
          delta_th = fabs(delta_th);
  
          if ( (delta_th * dR) < All.Highres_delta )
              zoomtarget = All.Highres_targetmass;
        }
    }
  
  if(can_this_cell_be_split(i) && P[i].Mass > 2.0 * zoomtarget)
    return 1;
  
  return 0;                     /* default is not to refine */
}
#endif

#ifdef SPIRAL
int refine_disczoom(int i)
{
  /* This routine refines a segment of a galactic disc that moves with the 
     mean gas. By Rowan Smith 2013. */
  double zoomtarget;
  double wr, thref, dtheta, theta_range, theta_range2,zrange,rinner,Rref,normalise;
  float rfloat, sp_ang;
  int revolve;

  Rref = 7.5 * KILOPARSEC;
  wr = -220.0 * 1.e5 / Rref;    /* vr of 220 kms at 7.5 kpc */
  thref = wr * All.Time * All.UnitTime_in_s;
  revolve = (int) (fabs(thref) / (2.0*M_PI));
  rfloat = (float) revolve;
  thref = thref + rfloat * 2.0 * M_PI;  /*Note this assumes the rotation is clockwise */

  theta_range = M_PI / 4.;
  theta_range2 = theta_range / 2.;
  zrange= 1. * KILOPARSEC/All.UnitLength_in_cm; /*Only refine 1kpc from disc plane*/
  rinner= 4. * KILOPARSEC/All.UnitLength_in_cm; /*Inner radius set to 4 kpc*/
  normalise= 2.3e-13; /*2.3e-13 gives 5 solar mass resolution*/

  double dx = P[i].Pos[0] - boxHalf_X;
  double dy = P[i].Pos[1] - boxHalf_Y;
  double dz = P[i].Pos[2] - boxHalf_Z;
  double rr = sqrt(dx*dx+dy*dy);

  zoomtarget = All.TargetGasMass;

  if (dz < zrange && dz > -1.0*zrange && rr > rinner)
    {
      /*  mpi_printf("Test values boxhalf x %g & z %g zrange %g, rinner %g, rr %g, dx %g, dz %g\n",boxHalf_X,boxHalf_Z,zrange,rinner,rr,dx,dz);*/

      dtheta = atan2(dy, dx) - thref;

      if(dtheta > 2.0 * M_PI)
        dtheta= dtheta - 2.0 *M_PI;

      if(dtheta > M_PI)
        dtheta = 2. * M_PI - dtheta;

      /*      if(dtheta < -1.0 * M_PI)
                {
                  dtheta = 2.0 * M_PI + dtheta;
                }*/
      sp_ang = fabs(dtheta);
#ifdef RAMP_REFINE
      if(sp_ang < theta_range && sp_ang > theta_range2)
        {
          zoomtarget = All.TargetGasMass / 1000.0 * (2490.0 * (sp_ang - M_PI / 8.0) + 5.0);
#ifdef SNE_RAMP_REFINE
          zoomtarget = All.TargetGasMass / 1000.0 * 2.0* (1245.0 * (sp_ang - M_PI / 8.0) + 5.0);  
#endif
        }
      if(zoomtarget < 100.0 * SOLAR_MASS / All.UnitMass_in_g)        /* to avoid sharp contrasts use exponential */
        zoomtarget = normalise * exp(78.2 * sp_ang) * SOLAR_MASS / All.UnitMass_in_g;

      /*Now assign the inner region where resolution is high */
      if(sp_ang <= theta_range2)
        {
          zoomtarget = normalise * exp(78.2 * theta_range2) * SOLAR_MASS / All.UnitMass_in_g;
#ifdef SNE_RAMP_REFINE
          zoomtarget = 10.0 * SOLAR_MASS / All.UnitMass_in_g;
#endif
        }
#else
      if(sp_ang < theta_range)
        zoomtarget = All.TargetGasMass / 10.0;
#endif
    }

  if(can_this_cell_be_split(i) && P[i].Mass > 2.0 * zoomtarget)
    return 1;

  return 0;                     /* default is not to refine */
}

#endif

#ifdef DISC_REFINE_ONLY
int refine_disconly(int i)
{
   double zoomtarget;
   double dx = P[i].Pos[0] - boxHalf_X;
   double dy = P[i].Pos[1] - boxHalf_Y;
   double dz = P[i].Pos[2] - boxHalf_Z; 
   double rr = sqrt(dx*dx+dy*dy);

   double zrange= KILOPARSEC/All.UnitLength_in_cm; /*Only refine 1kpc from plane*/
   double rinner= 3.*KILOPARSEC/All.UnitLength_in_cm;
   double router= 11.*KILOPARSEC/All.UnitLength_in_cm;

   zoomtarget=All.TargetGasMass;
   if (fabs(dz) < zrange && rr > rinner && rr < router)
     {
       zoomtarget=All.TargetGasMass/100.;
     }

      if(can_this_cell_be_split(i) && P[i].Mass > 2.0 * zoomtarget)
	return 1;
  return 0;                     /* default is not to refine */
}
#endif

#ifdef REFINEMENT_AROUND_DM
void dm_particle_create_list()
{
  struct refine_dm_data *DMPartList;
  DMPartList = (struct refine_dm_data *) mymalloc("DMPartList", All.TotPartDM * sizeof(struct refine_dm_data));


  int i, j, nsrc, nimport, ngrp;
  for(i = 0, nsrc = 0; i < NumPart; i++)
    {
      if(P[i].Type == 1)
        {
          DMPartList[nsrc].ID = P[i].ID;

          DMPartList[nsrc].pos[0] = P[i].Pos[0];
          DMPartList[nsrc].pos[1] = P[i].Pos[1];
          DMPartList[nsrc].pos[2] = P[i].Pos[2];

          //DMPartList[nsrc++].softening = All.SofteningTable[P[i].SofteningType];
          DMPartList[nsrc++].softening = All.ForceSoftening[P[i].SofteningType];
        }
    }

  for(j = 0; j < NTask; j++)
    Send_count[j] = nsrc;

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = 0;
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  /* exchange particle data */
  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&DMPartList[Send_offset[recvTask]],
                           Send_count[recvTask] * sizeof(struct refine_dm_data), MPI_BYTE,
                           recvTask, TAG_DENS_A,
                           &DMPartListGlobal[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct refine_dm_data), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(DMPartList);

#ifdef REFINEMENT_VOLUME_LIMIT
  /* now set minimum volume to refine to */
  /* get minimum softening length of all DM particles */
  double minrad = MAX_REAL_NUMBER;
  for(j = 0; j < All.TotPartDM; j++)
    if(minrad > DMPartListGlobal[j].softening / 10)
      minrad = DMPartListGlobal[j].softening / 10;
  /* set volume to 1% of corresponding sphere volume */
  All.MinVolume = 0.01 * 4. * M_PI / 3. * minrad * minrad * minrad;
  mpi_printf("Setting minimum volume to %e\n", All.MinVolume);
#endif
}

void dm_particle_update_list()
{
  struct refine_dm_data *DMPartList;
  DMPartList = (struct refine_dm_data *) mymalloc("DMPartList", All.TotPartDM * sizeof(struct refine_dm_data));

  int i, j, nsrc, nimport, ngrp;
  for(i = 0, nsrc = 0; i < NumPart; i++)
    {
      if(P[i].Type == 1)
        {
          DMPartList[nsrc].ID = P[i].ID;

          DMPartList[nsrc].pos[0] = P[i].Pos[0];
          DMPartList[nsrc].pos[1] = P[i].Pos[1];
          DMPartList[nsrc].pos[2] = P[i].Pos[2];

          //DMPartList[nsrc++].softening = All.SofteningTable[P[i].SofteningType];
          DMPartList[nsrc++].softening = All.ForceSoftening[P[i].SofteningType];
        }
    }

  for(j = 0; j < NTask; j++)
    Send_count[j] = nsrc;

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = 0;
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  /* exchange particle data */
  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&DMPartList[Send_offset[recvTask]],
                           Send_count[recvTask] * sizeof(struct refine_dm_data), MPI_BYTE,
                           recvTask, TAG_DENS_A,
                           &DMPartListGlobal[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct refine_dm_data), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(DMPartList);

}
#endif

#ifdef GMC_REFINEMENT
int gmc_refinement_criteria(int i)
  {
    /* PCC - 08.05.2013
       This file contains the different refinement criteria for
       modelling molecular clouds. At the moment we have a simple Jeans  
       criterion, but we should eventually add a shock criterion. Note
       that the Jeans length definition here is different (smaller) than
       that used in the standard Jeans refinement. This code also ignores
       the cell mass, unlike the other version. 
     */

    double jeans_temperature, jeans_length;
    double cells_per_jeans_length;
    double density;
    double three_over_four_pi, temp_prefactor;
    double cellrad;
    double numdens_ref_min, numdens_ref_max;
    double yn, energy, yntot;

    three_over_four_pi = 0.238732414637843;
    /* 5.*kb/2./gg/2.33/mp * T */
    temp_prefactor = 1.32740615361024e+15;

    /* make sure density is within range */
    density = SphP[i].Density * All.UnitDensity_in_cgs;
    if ( (density < All.GMCRefMinDensity) || (density > All.GMCRefMaxDensity) )
      return 0;

    /* Get the temperature of the gas  */
#ifdef SGCHEM
    yn = density / ((1.0 + 4.0 * ABHE) * PROTONMASS);
    yntot = (1 + ABHE - SphP[i].TracAbund[IH2] + SphP[i].TracAbund[IHP]) * yn;
#else
    yntot = density / 2.33 / PROTONMASS;
#endif

    energy =  (SphP[i].Utherm * All.UnitEnergy_in_cgs/All.UnitMass_in_g) * density;
    jeans_temperature = 2 * energy / (3 * yntot * BOLTZMANN);

    /* Calculate the number of cells per Jeans length */
    jeans_length = 2.0 * sqrt(three_over_four_pi/density) * sqrt(temp_prefactor*jeans_temperature);
    cellrad = 2.0 * get_cell_radius(i) * All.UnitLength_in_cm;
    cells_per_jeans_length = jeans_length / cellrad;

    /* Fix cells per jeans length here... */
    if ( cells_per_jeans_length < All.GMCRefCellsPerJeansLength )
      return(1);

    return(0);
  }
#endif

#endif
