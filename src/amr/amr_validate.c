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
#include <string.h>

#include "../allvars.h"
#include "../proto.h"

#ifdef AMR


void amr_validate_mesh(int node, int parent, int subnode)
{
  if(node < Ngb_MaxPart)
    {
      if(subnode & 1)
        {
          assert(Ngb_Nodes[parent].Center[0] == P[node].Pos[0] - amr_length[Mesh.DP[node].level + 1]);
        }
      else
        {
          assert(Ngb_Nodes[parent].Center[0] == P[node].Pos[0] + amr_length[Mesh.DP[node].level + 1]);
        }
#if NUMDIMS > 1
      if(subnode & 2)
        {
          assert(Ngb_Nodes[parent].Center[1] == P[node].Pos[1] - amr_length[Mesh.DP[node].level + 1]);
        }
      else
        {
          assert(Ngb_Nodes[parent].Center[1] == P[node].Pos[1] + amr_length[Mesh.DP[node].level + 1]);
        }
#if NUMDIMS>2
      if(subnode & 4)
        {
          assert(Ngb_Nodes[parent].Center[2] == P[node].Pos[2] - amr_length[Mesh.DP[node].level + 1]);
        }
      else
        {
          assert(Ngb_Nodes[parent].Center[2] == P[node].Pos[2] + amr_length[Mesh.DP[node].level + 1]);
        }
#endif
#endif

      assert(P[node].ID == Mesh.DP[node].ID);
    }
  else if(node < Ngb_MaxPart + Ngb_MaxNodes)
    {
      int i;
      for(i = 0; i < 1 << NUMDIMS; i++)
        {
          int child = Ngb_Nodes[node].u.suns[i];



          if(child < 0 && (node >= Ngb_FirstNonTopLevelNode || Ngb_DomainTask[node] == ThisTask))
            printf("missing child %d on task %d, node %d\n", i, ThisTask, node);

          if(child >= 0)
            amr_validate_mesh(child, node, i);
        }
    }
  else
    {
      //terminate("implement this!");
    }



}

void amr_validate_lists()
{

  int global_maxlevel;

  MPI_Allreduce(&Mesh.maxlevel, &global_maxlevel, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  int *visited = mymalloc("visited", sizeof(int) * (Ngb_MaxPart + Ngb_MaxNodes));
  memset(visited, 0, sizeof(int) * (Ngb_MaxPart + Ngb_MaxNodes));

  int level;

  for(level = global_maxlevel; level >= 0; level--)
    {
      int node = Mesh.lastinlevel[level];
      //int prevnode = -1;
      while(node >= 0)
        {
          if(visited[node] != 0)
            terminate("cycle at %d detected", node);
          visited[node] = 1;

          if(node >= Ngb_MaxPart)
            {
              if(Ngb_Nodes[node].level != level)
                terminate("node %d is in wrong level list", node);
              //prevnode = node;
              node = Ngb_Nodes[node].nextinlevel;
            }
          else
            {
              if(Mesh.DP[node].level != level)
                terminate("node %d is in wrong level list", node);
              //prevnode = node;
              node = Mesh.DP[node].nextinlevel;
            }
        }
    }

  myfree(visited);

}

#endif
