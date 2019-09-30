/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dg/dg_refinement.c
 * \date        08/2015
 * \author      Kevin Schaal
 * \brief       (de) refinement routines for the Galerkin (DG) module
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include "../allvars.h"
#include "../proto.h"

#ifdef DG

/*!
 * Project the hydro quantities of a node onto a a subnode
 * \param sun_num The position of the subnode (0: back left, 1: back right, 2: front left, 3: front right, 4: top back left, 5: top back right, 6: top front left, 7: top front right)
 * \param sun_index the particle index of the subnode
 */
void dg_split_hydro(int sun_num, int sun_index)
{
  assert(sun_index < Ngb_MaxPart);

  double (*P_X)[NOF_BASE_FUNCTIONS];

  switch (sun_num)
    {                           //speedup: P_X=sun_num
    case 0:
      P_X = P_A;
      break;
    case 1:
      P_X = P_B;
      break;
    case 2:
      P_X = P_C;
      break;
    case 3:
      P_X = P_D;
      break;
#ifndef TWODIMS
    case 4:
      P_X = P_E;
      break;
    case 5:
      P_X = P_F;
      break;
    case 6:
      P_X = P_G;
      break;
    case 7:
      P_X = P_H;
      break;
#endif
    default:
      terminate("ERROR in amr_accumulate_node_data: invalid subnode number!\n");
    }

  //store the old weights, set new weights to zero
  int l, j, k;

  for(l = 0; l < Nof_base_functions; l++)
    {
      for(k = 0; k < 5; k++)
        {
          Temp_weights[l][k] = SphP[sun_index].Weights[l][k];
          SphP[sun_index].Weights[l][k] = 0;
        }
    }


  //calculate new weights with a L2 projection
  for(l = 0; l < Nof_base_functions; l++)
    {
      for(j = 0; j < Nof_base_functions; j++)
        {
          SphP[sun_index].Weights[l][0] += Temp_weights[j][0] * P_X[j][l];
          SphP[sun_index].Weights[l][1] += Temp_weights[j][1] * P_X[j][l];
          SphP[sun_index].Weights[l][2] += Temp_weights[j][2] * P_X[j][l];
          SphP[sun_index].Weights[l][3] += Temp_weights[j][3] * P_X[j][l];
          SphP[sun_index].Weights[l][4] += Temp_weights[j][4] * P_X[j][l];
        }
    }

  //set the conserved variables

  P[sun_index].Mass = SphP[sun_index].Weights[0][W_RHO] * SphP[sun_index].Volume;
  SphP[sun_index].OldMass = 0.;
  SphP[sun_index].Momentum[0] = SphP[sun_index].Weights[0][W_PX] * SphP[sun_index].Volume;
  SphP[sun_index].Momentum[1] = SphP[sun_index].Weights[0][W_PY] * SphP[sun_index].Volume;
  SphP[sun_index].Momentum[2] = SphP[sun_index].Weights[0][W_PZ] * SphP[sun_index].Volume;
  SphP[sun_index].Energy = SphP[sun_index].Weights[0][W_E] * SphP[sun_index].Volume;
}


/*!
 * Project the hydro quantities of a subnode onto a node
 * \param new_cell The index of the new cell
 * \param sun_num The position of the subnode (0: lower left, 1: lower right, 2: upper left, 3: upper right)
 * \param sun_index the particle index
 */

void dg_sum_hydro(int new_cell, int sun_num, int sun_index)
{
  assert(new_cell < Ngb_MaxPart);
  assert(sun_index < Ngb_MaxPart);

  MyDouble(*P_X)[NOF_BASE_FUNCTIONS];
  int l, j, k;

  switch (sun_num)
    {                           //speedup: P_X=sun_num
    case 0:
      P_X = P_A;
      break;
    case 1:
      P_X = P_B;
      break;
    case 2:
      P_X = P_C;
      break;
    case 3:
      P_X = P_D;
      break;
#ifndef TWODIMS
    case 4:
      P_X = P_E;
      break;
    case 5:
      P_X = P_F;
      break;
    case 6:
      P_X = P_G;
      break;
    case 7:
      P_X = P_H;
      break;
#endif
    default:
      terminate("ERROR in dg_sum_hydro: invalid subnode number!\n");
    }

  if(sun_num == 0)              // => reset weights
    {

      for(l = 0; l < Nof_base_functions; l++)
        {
          for(k = 0; k < 5; k++)
            {
              Temp_weights[l][k] = 0;
            }
        }
    }

  for(l = 0; l < Nof_base_functions; l++)
    {
      for(j = 0; j < Nof_base_functions; j++)
        {

          Temp_weights[l][0] += 1. / DG_PROJ_NORM * SphP[sun_index].Weights[j][0] * P_X[l][j];
          Temp_weights[l][1] += 1. / DG_PROJ_NORM * SphP[sun_index].Weights[j][1] * P_X[l][j];
          Temp_weights[l][2] += 1. / DG_PROJ_NORM * SphP[sun_index].Weights[j][2] * P_X[l][j];
          Temp_weights[l][3] += 1. / DG_PROJ_NORM * SphP[sun_index].Weights[j][3] * P_X[l][j];
          Temp_weights[l][4] += 1. / DG_PROJ_NORM * SphP[sun_index].Weights[j][4] * P_X[l][j];
        }
    }

  if(sun_num == (2 << (NUMDIMS - 1)) - 1)       //3 for 2d, 7 for 3d
    {

      for(l = 0; l < Nof_base_functions; l++)
        {
          for(k = 0; k < 5; k++)
            {
              SphP[new_cell].Weights[l][k] = Temp_weights[l][k];
            }
        }

      P[new_cell].Mass = SphP[new_cell].Weights[0][W_RHO] * SphP[new_cell].Volume;
      SphP[new_cell].Momentum[0] = SphP[new_cell].Weights[0][W_PX] * SphP[new_cell].Volume;
      SphP[new_cell].Momentum[1] = SphP[new_cell].Weights[0][W_PY] * SphP[new_cell].Volume;
      SphP[new_cell].Momentum[2] = SphP[new_cell].Weights[0][W_PZ] * SphP[new_cell].Volume;
      SphP[new_cell].Energy = SphP[new_cell].Weights[0][W_E] * SphP[new_cell].Volume;

      if(P[new_cell].Mass < 0 || SphP[new_cell].Energy < 0)
        {
          terminate("ERROR in dg_sum_hydro: negative values: mass=%f, E=%f\n", P[new_cell].Mass, SphP[new_cell].Energy);
        }
    }
}

#endif
