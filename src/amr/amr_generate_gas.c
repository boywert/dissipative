/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/amr/amr_refinement.c
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
#ifdef GENERATE_GAS_IN_ICS

void amr_generate_gas_in_ics()
{
  int count;
  double fac, d, a, b, rho;

  if(RestartFlag == 0)
    {
      header.flag_entropy_instead_u = 0;

      Ngb_MaxPart = All.MaxPartSph;
      Ngb_MaxNodes = (int) (All.NgbTreeAllocFactor * All.MaxPartSph) + NTopnodes;
      ngb_treeallocate();
      ngb_treebuild(NumGas);

      MyIDType ids_offset = determine_ids_offset();

      fac = All.OmegaBaryon / All.Omega0;
      rho = All.Omega0 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

      int i;

      for(i = 0; i < NumPart; i++)
#ifdef SPLIT_PARTICLE_TYPE
        if((1 << P[i].Type) & (SPLIT_PARTICLE_TYPE))
#else
        if(P[i].Type == 1)
#endif
          {
            double cx, cy, cz, len, offset;
            peanokey key, morton;
            int k, subnode = 0, shift, rep;
            int th, nn, no;
            len = DomainLen;

            cx = DomainCenter[0];
            cy = DomainCenter[1];
            cz = DomainCenter[2];

            rep = 0;

            peano1D xb = domain_double_to_int(((P[i].Pos[0] - DomainCorner[0]) * DomainInverseLen) + 1.0);
            peano1D yb = domain_double_to_int(((P[i].Pos[1] - DomainCorner[1]) * DomainInverseLen) + 1.0);
            peano1D zb = domain_double_to_int(((P[i].Pos[2] - DomainCorner[2]) * DomainInverseLen) + 1.0);

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
                if(th >= Ngb_MaxPart)   /* we are dealing with a ghost node */
                  {
                    if(shift >= 0)
                      {
                        subnode = ((morton >> shift) & 7);
                      }
                    else
                      {
                        subnode = 0;
                        if(P[i].Pos[0] > cx)
                          subnode += 1;
                        if(P[i].Pos[1] > cy)
                          subnode += 2;
                        if(P[i].Pos[2] > cz)
                          subnode += 4;
                      }



                    nn = Ngb_Nodes[th].u.suns[subnode];

                    shift -= 3;


                    if(nn >= 0) /* ok, something is in the daughter slot already, need to continue */
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
                        terminate("");
                      }
                  }
                else
                  {
                    double dmass = P[i].Mass * fac;
                    double oldmass = P[no].Mass;

                    P[no].Mass += dmass;
                    P[i].Mass -= dmass;

                    P[no].Vel[0] = (oldmass * P[no].Vel[0] + dmass * P[i].Vel[0]) / P[no].Mass;
                    P[no].Vel[1] = (oldmass * P[no].Vel[1] + dmass * P[i].Vel[1]) / P[no].Mass;
                    P[no].Vel[2] = (oldmass * P[no].Vel[2] + dmass * P[i].Vel[2]) / P[no].Mass;

                  }
              }
          }

      ngb_treefree();

      All.MassTable[0] = 0;

#ifdef SPLIT_PARTICLE_TYPE
      for(i = 1; i < NTYPES; i++)
        if((1 << i) & (SPLIT_PARTICLE_TYPE))
          All.MassTable[i] *= (1 - fac);
#else
      All.MassTable[1] *= (1 - fac);
#endif
    }
}
#endif
#endif
