/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/amr/amr_refinement_criterion.c
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

#include "../allvars.h"
#include "../proto.h"

#ifdef AMR
#ifdef REFINEMENT

int amr_should_this_node_be_split(int no)
{
  if(no < Ngb_MaxPart)
    {
#ifdef REFINEMENT_SPLIT_CELLS
#ifdef DG

      if(fabs(SphP[no].Weights[1][W_RHO]) > All.DG_SlopeRangeFactor * All.DG_TargetSlope        //x slope
         || fabs(SphP[no].Weights[2][W_RHO]) > All.DG_SlopeRangeFactor * All.DG_TargetSlope     //y slope
         DG3D(||fabs(SphP[no].Weights[3][W_RHO]) > All.DG_SlopeRangeFactor * All.DG_TargetSlope))       //z slope
        {
          return 1;
        }
#else
      if(P[no].Mass > 2.0 * All.TargetGasMass)
        return 1;
#endif
#else
      return 0;
#endif
    }
  else if(no < Ngb_MaxPart + Ngb_MaxNodes)
    {

#ifdef REFINEMENT_MERGE_CELLS
#ifdef DG

#ifdef DEREFINE_X_LARGER

      if(Ngb_Nodes[no].Center[0] > DEREFINE_X_LARGER)
        {
          return 0;
        }
#endif

#ifdef DEREFINE_X_SMALLER

      if(Ngb_Nodes[no].Center[0] < DEREFINE_X_SMALLER)
        {
          return 0;
        }
#endif

#ifdef DEREFINE_Y_LARGER

      if(Ngb_Nodes[no].Center[1] > DEREFINE_Y_LARGER)
        {
          return 0;
        }
#endif

#ifdef DEREFINE_Y_SMALLER

      if(Ngb_Nodes[no].Center[1] < DEREFINE_Y_SMALLER)
        {
          return 0;
        }
#endif

#ifdef DEREFINE_Z_LARGER

      if(Ngb_Nodes[no].Center[2] > DEREFINE_Z_LARGER)
        {
          return 0;
        }
#endif

#ifdef DEREFINE_Z_SMALLER

      if(Ngb_Nodes[no].Center[2] < DEREFINE_Z_SMALLER)
        {
          return 0;
        }
#endif

#if defined(DEREFINE_X_LARGER) || defined(DEREFINE_X_SMALLER) || defined(DEREFINE_Y_LARGER) || defined(DEREFINE_Y_SMALLER) || defined(DEREFINE_Z_LARGER) || defined(DEREFINE_Z_SMALLER)
      return 1;
#endif

      if(fabs(Ngb_Nodes[no].hydro.Weights[1][W_RHO]) > All.DG_TargetSlope / All.DG_SlopeRangeFactor     //x slope
         || fabs(Ngb_Nodes[no].hydro.Weights[2][W_RHO]) > All.DG_TargetSlope / All.DG_SlopeRangeFactor  //y slope
         DG3D(||fabs(Ngb_Nodes[no].hydro.Weights[3][W_RHO]) > All.DG_TargetSlope / All.DG_SlopeRangeFactor))    //z slope
        {
          return 1;
        }
#else
      if(Ngb_Nodes[no].hydro.mass > 0.5 * All.TargetGasMass)    //if the node is heavy enough: keep it refined, otherwise: merge subnodes.
        return 1;
#endif
#else
      return 1;
#endif
    }

  return 0;
}
#endif
#endif
