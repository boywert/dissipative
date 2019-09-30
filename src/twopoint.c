/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/twopoint.c
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
#include <mpi.h>
#include <gsl/gsl_rng.h>

#include "allvars.h"
#include "proto.h"



/*! \file twopoint.c
 *  \brief computes the two-point mass correlation function on the fly
 */

static MyFloat *RsList;

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];
  MyFloat Rs;

  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

/* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in * in, int i, int firstnode)
{
#ifdef CELL_CENTER_GRAVITY
  if(P[i].Type == 0)
    {
      in->Pos[0] = SphP[i].Center[0];
      in->Pos[1] = SphP[i].Center[1];
      in->Pos[2] = SphP[i].Center[2];
    }
  else
#endif
    {
      in->Pos[0] = P[i].Pos[0];
      in->Pos[1] = P[i].Pos[1];
      in->Pos[2] = P[i].Pos[2];
    }

  in->Rs = RsList[i];

  in->Firstnode = firstnode;
}

/* local data structure that holds results acquired on remote processors */
typedef struct
{
  char unused;

} data_out;

static data_out *DataResult, *DataOut;

/* routine to store or combine result data */
static void out2particle(data_out * out, int i, int mode)
{

}

#include "generic_comm_helpers2.h"




#define FACT1 0.366025403785    /* FACT1 = 0.5 * (sqrt(3)-1) */
#define FACT2 0.86602540        /* FACT2 = 0.5 * sqrt(3) */

/* Note: This routine will only work correctly for particles of equal mass ! */


#define BINS_TP  40             /* number of bins used */
#define ALPHA  -1.0             /* slope used in randomly selecting radii around target particles */

#ifndef FRACTION_TP
#define FRACTION_TP  0.2
#endif /* fraction of particles selected for sphere
          placement. Will be scaled with total
          particle number so that a fixed value
          should give roughly the same noise level
          in the meaurement, indpendent of
          simulation size */



#define SQUARE_IT(x) ((x)*(x))



static double ThreadCount[NUM_THREADS][BINS_TP];
static double ThreadCountSpheres[NUM_THREADS][BINS_TP];
static double Count[BINS_TP];
static double CountSpheres[BINS_TP];
static double Xi[BINS_TP];
static double Rbin[BINS_TP];

static double R0, R1;           /* inner and outer radius for correlation function determination */

static double logR0;
static double binfac;
static double PartMass;
static double scaled_frac;


static int twopoint_count_local(int target, int mode, int threadid);


static void kernel_local(void)
{
  int i, bin;
  double p, rs;
#ifdef GENERIC_ASYNC
  int flag = 0;
#endif

#pragma omp parallel private(i, p, bin, rs)
  {
    int j, threadid = get_thread_num();
#ifdef GENERIC_ASYNC
    int count = 0;
#endif

    for(j = 0; j < NTask; j++)
      Thread[threadid].Exportflag[j] = -1;

    while(1)
      {
        if(Thread[threadid].ExportSpace < MinSpace)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              if(generic_polling_primary(count, NumPart))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        i = NextParticle++;

        if(i >= NumPart)
          break;

        if(gsl_rng_uniform(random_generator) < scaled_frac)
          {
            p = gsl_rng_uniform(random_generator);

            rs = pow(pow(R0, ALPHA) + p * (pow(R1, ALPHA) - pow(R0, ALPHA)), 1 / ALPHA);

            bin = (int) ((log(rs) - logR0) * binfac);

            rs = exp((bin + 1) / binfac + logR0);

            RsList[i] = rs;

            twopoint_count_local(i, MODE_LOCAL_PARTICLES, threadid);

            for(j = 0; j <= bin; j++)
              ThreadCountSpheres[threadid][j]++;
          }
      }
  }
}

static void kernel_imported(void)
{
  /* now do the particles that were sent to us */
  int i, cnt = 0;
#pragma omp parallel private(i)
  {
    int threadid = get_thread_num();
#ifdef GENERIC_ASYNC
    int count = 0;
#endif

    while(1)
      {
#pragma omp atomic capture
        i = cnt++;

        if(i >= Nimport)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              generic_polling_secondary();
          }

        count++;
#endif

        twopoint_count_local(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

/*  This function computes the two-point function.
 */
void twopoint(void)
{
  int i, k;
  double vol;
  double tstart, tend;
  double mass, masstot;
  void *state_buffer;


  if(ThisTask == 0)
    {
      printf("begin two-point correlation function...\n");
      fflush(stdout);
    }

  tstart = second();


  /* set inner and outer radius for the bins that are used for the correlation function estimate */
  R0 = get_default_softening_of_particletype(1);        /* we assume that type=1 is the primary type */
  R1 = All.BoxSize / 2;


  scaled_frac = FRACTION_TP * 1.0e7 / All.TotNumPart;

  logR0 = log(R0);
  binfac = BINS_TP / (log(R1) - log(R0));


  for(i = 0, mass = 0; i < NumPart; i++)
    mass += P[i].Mass;

  MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  PartMass = masstot / All.TotNumPart;

  for(k = 0; k < NUM_THREADS; k++)
    for(i = 0; i < BINS_TP; i++)
      {
        ThreadCount[k][i] = 0;
        ThreadCountSpheres[k][i] = 0;
      }

  /* allocate buffers to arrange communication */

  RsList = (MyFloat *) mymalloc("RsList", NumPart * sizeof(MyFloat));


  generic_set_MaxNexport();


  state_buffer = mymalloc("state_buffer", gsl_rng_size(random_generator));

  memcpy(state_buffer, gsl_rng_state(random_generator), gsl_rng_size(random_generator));

  gsl_rng_set(random_generator, P[0].ID + ThisTask);    /* seed things with first particle ID to make sure we are
                                                           different on each CPU */


  generic_comm_pattern(NumPart, kernel_local, kernel_imported);


  memcpy(gsl_rng_state(random_generator), state_buffer, gsl_rng_size(random_generator));
  myfree(state_buffer);



  myfree(RsList);


  /* Now compute the actual correlation function */

  for(k = 1; k < NUM_THREADS; k++)
    for(i = 0; i < BINS_TP; i++)
      {
        ThreadCount[0][i] += ThreadCount[k][i];
        ThreadCountSpheres[0][i] += ThreadCountSpheres[k][i];
      }

  MPI_Allreduce(ThreadCount[0], Count, BINS_TP, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(ThreadCountSpheres[0], CountSpheres, BINS_TP, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


  for(i = 0; i < BINS_TP; i++)
    {
      vol = 4 * M_PI / 3.0 * (pow(exp((i + 1.0) / binfac + logR0), 3) - pow(exp((i + 0.0) / binfac + logR0), 3));

      if(CountSpheres[i] > 0)
        Xi[i] = -1 + Count[i] / ((double) CountSpheres[i]) / (All.TotNumPart / pow(All.BoxSize, 3)) / vol;
      else
        Xi[i] = 0;

      Rbin[i] = exp((i + 0.5) / binfac + logR0);
    }

  twopoint_save();

  tend = second();

  mpi_printf("TWO-POINT: end two-point: Took=%g seconds.\n", timediff(tstart, tend));
}





/*! This function counts the pairs in a sphere
 */
int twopoint_count_local(int target, int mode, int threadid)
{
  int no, k, p, bin, bin2, numnodes, *firstnode;
  double r2, r, ri, ro;
  MyDouble dx, dy, dz, dist;
#ifdef PERIODIC
  MyDouble xtmp, ytmp, ztmp;
#endif
  MyDouble *pos;
  MyFloat rs;

  data_in local, *target_data;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      target_data = &local;

      numnodes = 1;
      firstnode = NULL;
    }
  else
    {
      target_data = &DataGet[target];

      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  pos = target_data->Pos;
  rs = target_data->Rs;


  /* Now start the actual tree-walk for this particle */

  for(k = 0; k < numnodes; k++)
    {
      if(mode == MODE_LOCAL_PARTICLES)
        {
          no = Tree_MaxPart;    /* root node */
        }
      else
        {
          no = firstnode[k];
          no = Nodes[no].u.d.nextnode;  /* open it */
        }

      while(no >= 0)
        {

          if(no < Tree_MaxPart) /* single particle */
            {
              p = no;
              no = Nextnode[no];

              dx = NGB_PERIODIC_LONG_X(Tree_Pos_list[3 * p + 0] - pos[0]);
              dy = NGB_PERIODIC_LONG_Y(Tree_Pos_list[3 * p + 1] - pos[1]);
              dz = NGB_PERIODIC_LONG_Z(Tree_Pos_list[3 * p + 2] - pos[2]);

            }
          else if(no < Tree_MaxPart + Tree_MaxNodes)    /* internal node */
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                {
                  if(no < Tree_FirstNonTopLevelNode)    /* we reached a top-level node again, which means that we are done with the branch */
                    break;
                }

              int current_no = no;
              struct NODE *current = &Nodes[no];

              no = current->u.d.sibling;        /* in case the node can be discarded */

              dist = rs + 0.5 * current->len;
              dx = NGB_PERIODIC_LONG_X(current->center[0] - pos[0]);
              if(dx > dist)
                continue;
              dy = NGB_PERIODIC_LONG_Y(current->center[1] - pos[1]);
              if(dy > dist)
                continue;
              dz = NGB_PERIODIC_LONG_Z(current->center[2] - pos[2]);
              if(dz > dist)
                continue;
              /* now test against the minimal sphere enclosing everything */
              dist += FACT1 * current->len;
              if((r2 = dx * dx + dy * dy + dz * dz) > dist * dist)
                continue;

              r = sqrt(r2);

              ri = r - FACT2 * current->len;
              ro = r + FACT2 * current->len;

              if(ri >= R0 && ro < R1)
                {
                  if(ro < rs)
                    {
                      bin = (int) ((log(ri) - logR0) * binfac);
                      bin2 = (int) ((log(ro) - logR0) * binfac);
                      if(bin == bin2)
                        {
                          if(mode == MODE_IMPORTED_PARTICLES)
                            {
                              if(current_no < Tree_FirstNonTopLevelNode)
                                continue;
                            }
                          ThreadCount[threadid][bin] += (long long) (current->u.d.mass / PartMass);
                          continue;
                        }
                    }
                }

              no = current->u.d.nextnode;       /* ok, we need to open the node */
              continue;
            }
          else if(no >= Tree_ImportedNodeOffset)        /* point from imported nodelist */
            {
              int n = no - Tree_ImportedNodeOffset;

              dx = NGB_PERIODIC_LONG_X(Tree_Points[n].Pos[0] - pos[0]);
              dy = NGB_PERIODIC_LONG_Y(Tree_Points[n].Pos[1] - pos[1]);
              dz = NGB_PERIODIC_LONG_Z(Tree_Points[n].Pos[2] - pos[2]);

              no = Nextnode[no - Tree_MaxNodes];
            }
          else                  /* pseudo particle */
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                terminate("mode == MODE_IMPORTED_PARTICLES");

              if(target >= 0)
                tree_treefind_export_node_threads(no, target, threadid);


              no = Nextnode[no - Tree_MaxNodes];
              continue;
            }

          r2 = dx * dx + dy * dy + dz * dz;

          if(r2 >= R0 * R0 && r2 < R1 * R1)
            {
              if(r2 < rs * rs)
                {
                  bin = (int) ((log(sqrt(r2)) - logR0) * binfac);
                  if(bin < BINS_TP)
                    ThreadCount[threadid][bin]++;
                }
            }
        }

    }


  return 0;
}




void twopoint_save(void)
{
  FILE *fd;
  char buf[500];
  int i;

  if(ThisTask == 0)
    {
      if(Argc < 6)
        sprintf(buf, "%s/correl_%03d.txt", All.OutputDir, RestartSnapNum);
      else
        {
          printf("Writing twopoint to %s\n", Argv[5]);
          sprintf(buf, "%s/correl_%03d.txt", Argv[5], RestartSnapNum);
        }

      if(!(fd = fopen(buf, "w")))
        {
          printf("can't open file `%s`\n", buf);
          endrun();
        }

      fprintf(fd, "%g\n", All.Time);
      i = BINS_TP;
      fprintf(fd, "%d\n", i);

      for(i = 0; i < BINS_TP; i++)
        fprintf(fd, "%g %g %g %g\n", Rbin[i], Xi[i], (double) Count[i], (double) CountSpheres[i]);

      fclose(fd);
    }
}
