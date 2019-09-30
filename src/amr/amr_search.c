/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/amr/amr_search.c
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
#include "../domain.h"

#ifdef AMR

/*!
 * \brief Perform a point location via tree walk, returns the task and the cell index on this task containing the point (x,y,z).
 *
 * \param x x-coordinate of the point we are looking for
 * \param y y-coordinate of the point we are looking for
 * \param z z-coordinate of the point we are looking for
 * \param task output parameter, the task which contains the point (x,y,z)
 * \param cell output parameter, the cell index on this task which contains the point (x,y,z)
 */


void amr_point_location(tessellation * T, double x, double y, double z, int *task, int *cell)
{
  if(ThisTask != 0)
    {
      terminate("Parallelize me!\n");
    }

  double cx, cy, cz, len, offset;
  peanokey key, morton;
  int subnode = 0, shift, rep;
  int th, nn, no;
  len = DomainLen;

  //global center
  cx = DomainCenter[0];
  cy = DomainCenter[1];
  cz = DomainCenter[2];

  rep = 0;

  int xb = domain_double_to_int(((x - DomainCorner[0]) * DomainInverseLen) + 1.0);
  int yb = domain_double_to_int(((y - DomainCorner[1]) * DomainInverseLen) + 1.0);
  int zb = domain_double_to_int(((z - DomainCorner[2]) * DomainInverseLen) + 1.0);

  key = peano_and_morton_key(xb, yb, zb, BITS_PER_DIMENSION, &morton);


  shift = 3 * (BITS_PER_DIMENSION - 1);

  no = 0;
  while(TopNodes[no].Daughter >= 0)
    {
      int sub = (key - TopNodes[no].StartKey) / (TopNodes[no].Size / 8);

      no = TopNodes[no].Daughter + sub;

      subnode = ((morton >> shift) & 7);

      shift -= 3;
      rep++;

      len *= 0.5;

      offset = 0.5 * len;

      if(subnode & 1)
        cx = cx + offset;
      else
        cx = cx - offset;

      if(subnode & 2)
        cy = cy + offset;
      else
        cy = cy - offset;

      if(subnode & 4)
        cz = cz + offset;
      else
        cz = cz - offset;
    }

  no = TopNodes[no].Leaf;
  th = Ngb_DomainNodeIndex[no];

  while(1)
    {
      if(th < Ngb_MaxPart)      //local cell
        {
          printf("Cell: (%f,%f,%f)\n", T->DP[th].x, T->DP[th].y, T->DP[th].z);
          return;
        }
      else if(th < Ngb_MaxPart + Ngb_MaxNodes)  //local node
        {
          if(shift >= 0)
            {
              subnode = ((morton >> shift) & 7);
            }
          else
            {
              subnode = 0;
              if(x > cx)
                subnode += 1;
              if(y > cy)
                subnode += 2;
              if(z > cz)
                subnode += 4;
            }

          nn = Ngb_Nodes[th].u.suns[subnode];

          shift -= 3;

          if(nn >= 0)           // subnode exists
            {
              th = nn;
              rep++;

              len *= 0.5;

              offset = 0.5 * len;

              if(subnode & 1)
                cx = cx + offset;
              else
                cx = cx - offset;

              if(subnode & 2)
                cy = cy + offset;
              else
                cy = cy - offset;

              if(subnode & 4)
                cz = cz + offset;
              else
                cz = cz - offset;
            }
          else
            {
              assert(0);
            }
        }
      else                      //ghost cell/node
        {
          terminate("Parallelize me!\n");
        }
    }                           //end while

  return;

}

#endif //AMR
