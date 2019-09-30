/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/amr/amr_ngb.c
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

void amr_link_ngb_node(int new_node, int node, int subnode)
{
  int r;
  int ngb;
  int dir;

  for(dir = 0; dir < NUMDIMS; dir++)
    {
#if defined(LONG_X) || defined(LONG_Y) ||defined(LONG_Z)
      if(Ngb_Nodes[node].level < amr_last_long_level[dir])
        {
          Ngb_Nodes[new_node].neighbors[2 * dir] = new_node;
          Ngb_Nodes[new_node].neighbors[2 * dir + 1] = new_node;

          continue;
        }
#endif

      int ref_boundary = 0;

      r = (subnode >> dir) % 2; //left or right side?

      //Handle the external side
      ngb = Ngb_Nodes[node].neighbors[2 * dir + r];
#ifdef REFLECTIVE_X
      if(dir == 0 && ngb == node)
        {
          ngb = new_node;
          ref_boundary = 1;
        }
#endif
#ifdef REFLECTIVE_Y
      if(dir == 1 && ngb == node)
        {
          ngb = new_node;
          ref_boundary = 1;
        }
#endif
#ifdef REFLECTIVE_Z
      if(dir == 2 && ngb == node)
        {
          ngb = new_node;
          ref_boundary = 1;
        }
#endif
      if(ref_boundary == 0)
        {
          if(ngb >= Ngb_FirstNonTopLevelNode || (ngb >= Ngb_MaxPart && ngb < Ngb_FirstNonTopLevelNode && (Ngb_DomainTask[ngb] == ThisTask || Ngb_DomainTask[ngb] == -1)))       // a local neighbor node, follow to is daughter
            {
              ngb = Ngb_Nodes[ngb].u.suns[subnode ^ (1 << dir)];
            }
        }

      Ngb_Nodes[new_node].neighbors[2 * dir + r] = ngb;
      if(ngb < 0)
        terminate("e1");

      //handle the internal side
      ngb = Ngb_Nodes[node].u.suns[subnode ^ (1 << dir)];

      Ngb_Nodes[new_node].neighbors[2 * dir + (r ^ 1)] = ngb;
      if(ngb < 0)
        terminate("e2");
    }
}

void amr_link_ngb_particle(int particle, int node, int subnode)
{
  int r;
  int ngb;
  int dir;

  for(dir = 0; dir < NUMDIMS; dir++)
    {
      int ref_boundary = 0;

      r = (subnode >> dir) % 2; //left or right side?

      //Handle the external side
      ngb = Ngb_Nodes[node].neighbors[2 * dir + r];
#ifdef REFLECTIVE_X
      if(dir == 0 && ngb == node)
        {
          ngb = particle;
          ref_boundary = 1;
        }
#endif
#ifdef REFLECTIVE_Y
      if(dir == 1 && ngb == node)
        {
          ngb = particle;
          ref_boundary = 1;
        }
#endif
#ifdef REFLECTIVE_Z
      if(dir == 2 && ngb == node)
        {
          ngb = particle;
          ref_boundary = 1;
        }
#endif
      if(ref_boundary == 0)
        {
          if(ngb >= Ngb_FirstNonTopLevelNode || (ngb >= Ngb_MaxPart && (Ngb_DomainTask[ngb] == ThisTask || Ngb_DomainTask[ngb] == -1))) // a local neighbor node, follow to is daughter, or a leaf top node as neighbor
            {
              ngb = Ngb_Nodes[ngb].u.suns[subnode ^ (1 << dir)];
            }
        }

      Mesh.DP[particle].neighbors[2 * dir + r] = ngb;


      //handle the internal side
      ngb = Ngb_Nodes[node].u.suns[subnode ^ (1 << dir)];
      Mesh.DP[particle].neighbors[2 * dir + (r ^ 1)] = ngb;
    }
}

void amr_link_ngb_particle_ghost(int particle, int node, int subnode)
{
  int r;
  int ngb;
  int dir;

  for(dir = 0; dir < NUMDIMS; dir++)
    {
      r = (subnode >> dir) % 2; //left or right side?
      //Handle the external side
      ngb = Ngb_Nodes[node].neighbors[2 * dir + r];
      if(ngb >= Ngb_MaxPart && ngb < Ngb_MaxPart + Mesh.nodes_total)
        ngb = Ngb_Nodes[ngb].u.suns[subnode ^ (1 << dir)];

      Mesh.DP[particle - Mesh.nodes_total].neighbors[2 * dir + r] = ngb;

      if(ngb >= 0)
        {
          if(ngb < Ngb_MaxPart && Mesh.DP[ngb].level == Mesh.DP[particle - Mesh.nodes_total].level)
            {
              Mesh.DP[ngb].neighbors[2 * dir + (r ^ 1)] = particle;
            }
          else if(ngb >= Ngb_MaxPart && ngb < Ngb_MaxPart + Ngb_MaxNodes && (ngb >= Ngb_FirstNonTopLevelNode || (Ngb_DomainTask[ngb] == ThisTask || Ngb_DomainTask[ngb] == -1)))        // a local node
            {
              Ngb_Nodes[ngb].neighbors[2 * dir + (r ^ 1)] = particle;

              int i = 0;
              int j = 0;

#if NUMDIMS >= 2
              for(i = 0; i < 2; i++)
#endif
#if NUMDIMS == 3
                for(j = 0; j < 2; j++)
#endif
                  {
                    {
                      int subnode = (i << ((dir + 1) % NUMDIMS)) + (j << ((dir + 2) % NUMDIMS)) + ((r ^ 1) << dir);
                      int cell = Ngb_Nodes[ngb].u.suns[subnode];
                      assert(cell < Ngb_MaxPart);       //must be a local particle, otherwise the mesh is broken
                      Mesh.DP[cell].neighbors[2 * dir + (r ^ 1)] = particle;
                    }
                  }

            }
        }

      //handle the internal side
      ngb = Ngb_Nodes[node].u.suns[subnode ^ (1 << dir)];

      Mesh.DP[particle - Mesh.nodes_total].neighbors[2 * dir + (r ^ 1)] = ngb;
    }
}

void amr_link_ngb_node_ghost(int new_node, int node, int subnode)
{
  int r;
  int ngb;
  int dir;

  for(dir = 0; dir < NUMDIMS; dir++)
    {
      r = (subnode >> dir) % 2; //left or right side?
      //Handle the external side
      ngb = Ngb_Nodes[node].neighbors[2 * dir + r];
      if(ngb >= Ngb_MaxPart && ngb < Ngb_MaxPart + Mesh.nodes_total)
        ngb = Ngb_Nodes[ngb].u.suns[subnode ^ (1 << dir)];

      Ngb_Nodes[new_node].neighbors[2 * dir + r] = ngb;

      if(ngb >= 0)
        {
          if(ngb < Ngb_MaxPart)
            {
              Mesh.DP[ngb].neighbors[2 * dir + (r ^ 1)] = new_node;
            }
          else if(ngb < Ngb_MaxPart + Ngb_MaxNodes)
            {
              Ngb_Nodes[ngb].neighbors[2 * dir + (r ^ 1)] = new_node;
            }
        }


      //handle the internal side
      ngb = Ngb_Nodes[node].u.suns[subnode ^ (1 << dir)];

      Ngb_Nodes[new_node].neighbors[2 * dir + (r ^ 1)] = ngb;
    }

}


void amr_link_ngb(int no)
{
  int j, p;

  if(no >= Ngb_MaxPart)         /* internal node */
    {


      for(j = 0; j < (1 << NUMDIMS); j++)
        {
          p = Ngb_Nodes[no].u.suns[j];

          //skip local particles and nodes (already linked)
          if(p >= 0 && p < Ngb_MaxPart)
            {
              continue;
            }

          else if(p >= Ngb_FirstNonTopLevelNode && p < Ngb_MaxPart + Ngb_MaxNodes)
            {
              continue;
            }

          else if(p >= Ngb_MaxPart + Ngb_MaxNodes && p < Ngb_MaxPart + Mesh.nodes_total)
            {
              amr_link_ngb_node_ghost(p, no, j);
            }
          else if(p >= Ngb_MaxPart + Mesh.nodes_total)
            {
              amr_link_ngb_particle_ghost(p, no, j);
              continue;
            }

          amr_link_ngb(p);
        }
    }
}


void amr_link_toptree(int no)
{
  int j, p;

  if(no >= Ngb_MaxPart && no < Ngb_MaxPart + Ngb_MaxNodes)      /* internal node */
    {
      if(!(no >= Ngb_MaxPart && no < Ngb_FirstNonTopLevelNode))
        terminate("can't be");

      for(j = 0; j < (1 << NUMDIMS); j++)
        {
          if((p = Ngb_Nodes[no].u.suns[j]) >= 0)
            {

              if(p >= Ngb_MaxPart)
                {
                  if(p < Ngb_MaxPart + Ngb_MaxNodes)    // not a pseudo particle
                    {

                      if(p < Ngb_FirstNonTopLevelNode)
                        {
                          amr_set_node_data(p, no, j, 0, &Mesh);
                          amr_link_ngb_node(p, no, j);

                          amr_link_toptree(p);
                        }
                    }
                }
            }
          else
            {
              if(no >= Ngb_FirstNonTopLevelNode)
                {
                  terminate("invalid amr mesh: not all subnodes are occupied in node %d\n", no);
                }
            }
        }
    }
  else
    {
      terminate("bad");
    }
}
#endif //AMR
