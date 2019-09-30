/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/tracer_trajectory.c
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
#include <sys/stat.h>
#include "allvars.h"
#include "proto.h"

#ifdef TRACER_TRAJECTORY

typedef struct
{
  float Pos[3];
  float Rho;
  float Temp;
  float Utherm;
  MyIDType ID;
#ifdef TRACER_TRAJECTORY_EXTENDED_OUTPUT
  float Composition[EOS_NSPECIES];
  float Dedt;
  float Vel[3];
#endif
} tracer_data;

#ifdef TRACER_TRAJECTORY_GENERATE
static double generate_tracers(void);
#else
static double load_tracers(void);
#endif
static void generate_single_tracer_in_cell(int idx);

void tracer_init(void)
{
#ifdef TRACER_TRAJECTORY_GENERATE
#ifndef TETRA_INDEX_IN_FACE
#error "The option TRACER_TRAJECTORY_GENERATE requires TETRA_INDEX_IN_FACE."
#endif
  double TracerMass = generate_tracers();
#else
  double TracerMass = load_tracers();
#endif

  tracer_init_output(TracerMass);

  All.TracerNextOutput = -1.;
}

#ifdef TRACER_TRAJECTORY_GENERATE
double generate_tracers(void)
{
  mpi_printf("Generating %d tracer particles...\n", All.NumberOfTracersToGenerate);

  int i;
  double mTask = 0;
  for(i = 0; i < NumGas; i++)
    mTask += P[i].Mass;

  int *nTracerTasks = NULL;
  double *mTasks = NULL;
  double tracerMass;

  if(ThisTask == 0)
    {
      mTasks = (double *) mymalloc("mTasks", NTask * sizeof(double));
      nTracerTasks = (int *) mymalloc("nTracerTask", NTask * sizeof(int));
    }

  MPI_Gather(&mTask, 1, MPI_DOUBLE, mTasks, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      double mtot = 0;
      for(i = 0; i < NTask; i++)
        mtot += mTasks[i];

      printf("mtot=%g\n", mtot);
      tracerMass = mtot / All.NumberOfTracersToGenerate;

      int count = 0;
      for(i = 0; i < NTask; i++)
        {
          nTracerTasks[i] = floor(mTasks[i] / tracerMass);
          count += nTracerTasks[i];
        }

      for(i = 1; i < NTask; i++)
        mTasks[i] += mTasks[i - 1];     // now holds cumulative masse

      int left = All.NumberOfTracersToGenerate - count;

      for(i = 0; i < left; i++)
        {
          double mPos = get_random_number() * mtot;

          if(mPos < mTasks[0])
            {
              nTracerTasks[0]++;
              continue;
            }

          int Left = 0;
          int Right = NTask - 1;

          while(!(Left + 1 == Right && mTasks[Left] <= mPos && mTasks[Right] > mPos))
            {
              int idx = (Left + Right) / 2;
              if(mTasks[idx] <= mPos)
                Left = idx;
              else if(mTasks[idx] > mPos)
                Right = idx;
            }

          nTracerTasks[Right]++;
        }
    }

  MPI_Bcast(&tracerMass, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  mpi_printf("TracerMass=%g\n", tracerMass);

  int nTracersLocal;
  MPI_Scatter(nTracerTasks, 1, MPI_INT, &nTracersLocal, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      myfree(nTracerTasks);
      myfree(mTasks);
    }

  double *mCumulative = (double *) mymalloc("mCumulative", NumGas * sizeof(double));
  mCumulative[0] = P[0].Mass;
  for(i = 1; i < NumGas; i++)
    mCumulative[i] = mCumulative[i - 1] + P[i].Mass;

  double mtot = mCumulative[NumGas - 1];

  mpi_printf("Starting to %d generate tracers locally...\n", nTracersLocal);

  for(i = 0; i < nTracersLocal; i++)
    {
      double mPos = get_random_number() * mtot;

      if(mPos >= mCumulative[0])
        {
          int Left = 0;
          int Right = NumGas - 1;

          while(!(Left + 1 == Right && mCumulative[Left] <= mPos && mCumulative[Right] > mPos))
            {
              int idx = (Left + Right) / 2;
              if(mCumulative[idx] <= mPos)
                Left = idx;
              else if(mCumulative[idx] > mPos)
                Right = idx;
            }

          generate_single_tracer_in_cell(Right);
        }
      else
        generate_single_tracer_in_cell(0);
      P[NumPart - 1].ID = 2000000000 + ThisTask * All.NumberOfTracersToGenerate + i;
    }

  myfree(mCumulative);

  All.TotNumPart += All.NumberOfTracersToGenerate;
  All.NTracerTot = All.NumberOfTracersToGenerate;

  mpi_printf("All tracers generated.\n");

  return tracerMass;
}

void generate_single_tracer_in_cell(int idx)
{
  const int edge_start[6] = { 0, 0, 0, 1, 1, 2 };
  const int edge_end[6] = { 1, 2, 3, 2, 3, 3 };
  const int edge_opposite[6] = { 3, 1, 2, 3, 0, 1 };
  const int edge_nexttetra[6] = { 2, 3, 1, 0, 2, 0 };

  /* find all corners of the cell to determine max search radius for random points */
  double min[3], max[3];

  int i;
  for(i = 0; i < 3; i++)
    {
      min[i] = +MAX_DOUBLE_NUMBER;
      max[i] = -MAX_DOUBLE_NUMBER;
    }

  int q = SphP[idx].first_connection;

  while(q >= 0)
    {
      int vf = DC[q].vf_index;
      int p1 = Mesh.VF[vf].p1;
      int p2 = Mesh.VF[vf].p2;

      if(p1 < 0 || p2 < 0)
        {
          q = DC[q].next;
          continue;
        }

      int tt = Mesh.VF[vf].dt_index;
      tetra *t = &Mesh.DT[tt];

      int nr;
      for(nr = 0; nr < 6; nr++)
        {
          int start_index = t->p[edge_start[nr]];
          int end_index = t->p[edge_end[nr]];

          if((start_index == p1 && end_index == p2) || (start_index == p2 && end_index == p1))
            break;
        }

      if(nr == 6)
        terminate("fail");

      tetra *prev, *next;
      tetra_center *nextc;

      int i, j, k, l, m, ii, jj, kk, ll, nn;

      i = edge_start[nr];
      j = edge_end[nr];
      k = edge_opposite[nr];
      l = edge_nexttetra[nr];

      prev = t;

      do
        {
          nn = prev->t[l];
          next = &Mesh.DT[nn];
          nextc = &Mesh.DTC[nn];

          if(nextc->cx < min[0])
            min[0] = nextc->cx;
          if(nextc->cy < min[1])
            min[1] = nextc->cy;
          if(nextc->cz < min[2])
            min[2] = nextc->cz;

          if(nextc->cx > max[0])
            max[0] = nextc->cx;
          if(nextc->cy > max[1])
            max[1] = nextc->cy;
          if(nextc->cz > max[2])
            max[2] = nextc->cz;

          for(m = 0, ll = ii = jj = -1; m < 4; m++)
            {
              if(next->p[m] == prev->p[k])
                ll = m;
              if(next->p[m] == prev->p[i])
                ii = m;
              if(next->p[m] == prev->p[j])
                jj = m;
            }

          if(ll < 0 || ii < 0 || jj < 0)
            terminate("inconsistency");

          kk = 6 - (ll + ii + jj);

          prev = next;

          i = ii;
          l = ll;
          j = jj;
          k = kk;
        }
      while(next != t);

      if(q == SphP[idx].last_connection)
        break;

      q = DC[q].next;
    }

  /* find random point in cell */
  double dx = max[0] - min[0];
  double dy = max[1] - min[1];
  double dz = max[2] - min[2];

  double px, py, pz, dist2;
  int found = 0, iter = 0;
  while(!found)
    {
      px = min[0] + get_random_number() * dx;
      py = min[1] + get_random_number() * dy;
      pz = min[2] + get_random_number() * dz;

      found = 1;
      dist2 = (px - P[idx].Pos[0]) * (px - P[idx].Pos[0]) + (py - P[idx].Pos[1]) * (py - P[idx].Pos[1]) + (pz - P[idx].Pos[2]) * (pz - P[idx].Pos[2]);

      int q = SphP[idx].first_connection;
      while(q >= 0)
        {
          int dp = DC[q].dp_index;
          int particle = Mesh.DP[dp].index;

          if(particle < 0)
            {
              q = DC[q].next;
              continue;
            }

          double d2 = (px - Mesh.DP[dp].x) * (px - Mesh.DP[dp].x) + (py - Mesh.DP[dp].y) * (py - Mesh.DP[dp].y) + (pz - Mesh.DP[dp].z) * (pz - Mesh.DP[dp].z);

          if(d2 < dist2)
            {
              /* found a neighbor that is closer to this point -> it is not in our cell */
              found = 0;
              if(iter >= 80)
                printf("fail: d=%g, dist=%g, px=%g, py=%g, pz=%g\n", sqrt(d2), sqrt(dist2), px, py, pz);
              break;
            }

          if(q == SphP[idx].last_connection)
            break;

          q = DC[q].next;
        }

      iter++;
      if(iter > 100)
        {
          printf("There is a problem with cell %d (ID=%d): pos=%g,%g,%g min=%g,%g,%g max=%g,%g,%g, cellvol=%g, searchvol=%g\n", idx, P[idx].ID,
                 P[idx].Pos[0], P[idx].Pos[1], P[idx].Pos[2], min[0], min[1], min[2], max[0], max[1], max[2], SphP[idx].Volume, dx * dy * dz);
          terminate("fail");
        }
    }

  if(NumPart == All.MaxPart)
    terminate("No free particle slots left.");

  int p = NumPart;

  P[p].Pos[0] = px;
  P[p].Pos[1] = py;
  P[p].Pos[2] = pz;
  P[p].Vel[0] = 0;
  P[p].Vel[1] = 0;
  P[p].Vel[2] = 0;
  P[p].Mass = 0;
  P[p].Type = TRACER_PARTICLE;
  P[p].TracerHsml = 1.1 * sqrt(dist2);

  P[p].TimeBinHydro = 0;
  P[p].TimeBinGrav = 0;

  NumPart++;
}

#else

double load_tracers()
{
  int NTracerTot, NTracer, Left, NTracerPart;
  int part, p;
  double Mass, MassAll, TracerMass;
  double *pos;
  size_t start;
  MyIDType IDOffset;

  FILE *fp = fopen(All.TracerInitFile, "r");
  fread(&NTracerTot, sizeof(int), 1, fp);

  mpi_printf("Found %d tracers to load.\n", NTracerTot);

  NTracerPart = NTracerTot / NTask;
  Left = NTracerTot - NTracerPart * NTask;

  NTracer = NTracerPart;
  if(ThisTask == NTask - 1)
    NTracer += Left;

  IDOffset = 2000000000 + ThisTask * NTracerTot;

  Mass = 0;
  for(part = 0; part < NumGas; part++)
    Mass += P[part].Mass;

  MPI_Allreduce(&Mass, &MassAll, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  TracerMass = MassAll / NTracerTot;
  mpi_printf("Tracer Mass: %g\n", TracerMass);

  start = NTracerPart * ThisTask * 3 * sizeof(double) + 4;
  pos = (double *) mymalloc("tracer_pos", NTracer * 3 * sizeof(double));

  fseek(fp, start, 0);
  fread(pos, sizeof(double), NTracer * 3, fp);
  fclose(fp);

  for(part = 0; part < NTracer; part++)
    {
      if(NumPart == All.MaxPart)
        terminate("No free particle slots left.");

      p = NumPart;

      P[p].Pos[0] = pos[part * 3];
      P[p].Pos[1] = pos[part * 3 + 1];
      P[p].Pos[2] = pos[part * 3 + 2];
      P[p].Vel[0] = 0;
      P[p].Vel[1] = 0;
      P[p].Vel[2] = 0;
      P[p].Mass = 0;
      P[p].Type = TRACER_PARTICLE;
      P[p].ID = IDOffset + part;
      P[p].TracerHsml = 1e7;

      P[p].TimeBinHydro = 0;
      P[p].TimeBinGrav = 0;

      NumPart++;
    }

  All.TotNumPart += NTracerTot;
  All.MassTable[TRACER_PARTICLE] = 0;

  myfree(pos);

  All.NTracerTot = NTracerTot;
  mpi_printf("%d Tracer initialized.\n", All.NTracerTot);

  return TracerMass;
}

#endif /* TRACER_TRAJECTORY_GENERATE */

void tracer_init_output_configuration(void)
{
  if(ThisTask == 0)
    {
      FILE *fp;
      if(!(fp = fopen(All.TracerOutputConfFile, "r")))
        terminate("Tracer Output Configuration File %s missing.\n", All.TracerOutputConfFile);

      fscanf(fp, "%d", &All.NTracerOutputSteps);

      All.NTracerOutputTimes = (double *) malloc(sizeof(double) * All.NTracerOutputSteps);
      All.NTracerOutputTimePeriod = (double *) malloc(sizeof(double) * All.NTracerOutputSteps);

      int i;
      for(i = 0; i < All.NTracerOutputSteps; i++)
        {
          fscanf(fp, "%lf %lf", &All.NTracerOutputTimes[i], &All.NTracerOutputTimePeriod[i]);
        }

      fclose(fp);
    }

  MPI_Bcast(&All.NTracerOutputSteps, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if(ThisTask != 0)
    {
      All.NTracerOutputTimes = (double *) malloc(sizeof(double) * All.NTracerOutputSteps);
      All.NTracerOutputTimePeriod = (double *) malloc(sizeof(double) * All.NTracerOutputSteps);
    }

  MPI_Bcast(All.NTracerOutputTimes, All.NTracerOutputSteps, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(All.NTracerOutputTimePeriod, All.NTracerOutputSteps, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      printf("TRACER_TRAJECTORY: %d steps.\n", All.NTracerOutputSteps);
      int i;
      for(i = 0; i < All.NTracerOutputSteps; i++)
        {
          printf("TRACER_TRAJECTORY: Until t=%g, we write a tracer output every dt=%g.\n", All.NTracerOutputTimes[i], All.NTracerOutputTimePeriod[i]);
        }
    }

  if(All.TracerNextOutput < 0)
    {
      All.TracerNextOutput = All.Time;
      tracer_write_output_if_needed();
    }
}

void tracer_init_output(double tmass)
{
  if(ThisTask == 0)
    {
      char buf[1000];
      sprintf(buf, "%s/%s", All.OutputDir, All.TracerOutputFile);

      int NTracer = All.NTracerTot;
      double *tmasses = (double *) mymalloc("tmasses", NTracer * sizeof(double));

      int tr;
      for(tr = 0; tr < NTracer; tr++)
        tmasses[tr] = tmass;

      mkdir(All.OutputDir, 02755);

      FILE *fd = fopen(buf, "w");

      int dummy = 4;
      fwrite(&dummy, sizeof(int), 1, fd);
      fwrite(&NTracer, sizeof(int), 1, fd);
      fwrite(&dummy, sizeof(int), 1, fd);

      dummy = NTracer * sizeof(double);
      fwrite(&dummy, sizeof(int), 1, fd);
      fwrite(tmasses, sizeof(double), NTracer, fd);
      fwrite(&dummy, sizeof(int), 1, fd);

      fclose(fd);

      myfree(tmasses);

      printf("TRACER_TRAJECTORY: Tracer output initialized.\n");
    }
}

int compare_tracers(const void *a, const void *b)
{
  if(((tracer_data *) a)->ID < ((tracer_data *) b)->ID)
    return -1;

  if(((tracer_data *) a)->ID > ((tracer_data *) b)->ID)
    return +1;

  return 0;
}

void tracer_write_output_if_needed(void)
{
  if(All.Time > All.TracerNextOutput)
    {
      tracer_write_output();

      int idx = 0;
      while(idx < All.NTracerOutputSteps - 1 && All.Time >= All.NTracerOutputTimes[idx])
        idx++;

      while(All.TracerNextOutput < All.Time)
        All.TracerNextOutput += All.NTracerOutputTimePeriod[idx];
      mpi_printf("TRACER_TRAJECTORY: Wrote tracer output, will write next at t=%g.\n", All.TracerNextOutput);
    }
}

void tracer_write_output(void)
{
  tracer_data *tdata;
  int NTracer = All.NTracerTot;
  int nt, p;

  tdata = (tracer_data *) mymalloc("tdata", NTracer * sizeof(tracer_data));

  nt = 0;
  for(p = 0; p < NumPart; p++)
    {
      if(P[p].Type == TRACER_PARTICLE)
        {
          tdata[nt].Pos[0] = P[p].Pos[0];
          tdata[nt].Pos[1] = P[p].Pos[1];
          tdata[nt].Pos[2] = P[p].Pos[2];
          tdata[nt].Rho = P[p].tRho;
          tdata[nt].Temp = P[p].tTemp;
          tdata[nt].Utherm = P[p].tUtherm;
          tdata[nt].ID = P[p].ID;
#ifdef TRACER_TRAJECTORY_EXTENDED_OUTPUT
          int k;
          for(k = 0; k < EOS_NSPECIES; k++)
            tdata[nt].Composition[k] = P[p].tComposition[k];
          tdata[nt].Dedt = P[p].tDedt;
          tdata[nt].Vel[0] = P[p].Vel[0];
          tdata[nt].Vel[1] = P[p].Vel[1];
          tdata[nt].Vel[2] = P[p].Vel[2];
#endif
          nt++;
        }
    }

  if(ThisTask == 0)
    {
      MPI_Status status;
      int Task, count;
      for(Task = 1; Task < NTask; Task++)
        {
          MPI_Recv(&count, 1, MPI_INT, Task, 666, MPI_COMM_WORLD, &status);
          MPI_Recv(&tdata[nt], count * sizeof(tracer_data), MPI_BYTE, Task, 667, MPI_COMM_WORLD, &status);
          nt += count;
        }
    }
  else
    {
      MPI_Send(&nt, 1, MPI_INT, 0, 666, MPI_COMM_WORLD);
      MPI_Send(tdata, nt * sizeof(tracer_data), MPI_BYTE, 0, 667, MPI_COMM_WORLD);
    }

  if(ThisTask == 0)
    {
      printf("TRACER_TRAJECTORY: Collected %d tracer particles.\n", nt);

      mysort(tdata, NTracer, sizeof(tracer_data), compare_tracers);

      char buf[500];
      sprintf(buf, "%s/tracer.dat", All.OutputDir);

      FILE *fd = fopen(buf, "a");
      printf("TRACER_TRAJECTORY: Opened file %s. Start writing at position %ld.\n", buf, ftell(fd));

#ifndef TRACER_TRAJECTORY_EXTENDED_OUTPUT
      int dummy = NTracer * 6 * sizeof(float) + 8;
#else
      int dummy = NTracer * (10 + EOS_NSPECIES) * sizeof(float) + 8;
#endif
      fwrite(&dummy, sizeof(int), 1, fd);
      fwrite(&All.Time, sizeof(double), 1, fd);

      float *data = mymalloc("tdata_single", NTracer * sizeof(float));

      for(nt = 0; nt < NTracer; nt++)
        data[nt] = (float) tdata[nt].Pos[0];
      fwrite(data, sizeof(float), NTracer, fd);

      for(nt = 0; nt < NTracer; nt++)
        data[nt] = (float) tdata[nt].Pos[1];
      fwrite(data, sizeof(float), NTracer, fd);

      for(nt = 0; nt < NTracer; nt++)
        data[nt] = (float) tdata[nt].Pos[2];
      fwrite(data, sizeof(float), NTracer, fd);

      for(nt = 0; nt < NTracer; nt++)
        data[nt] = (float) tdata[nt].Rho;
      fwrite(data, sizeof(float), NTracer, fd);

      for(nt = 0; nt < NTracer; nt++)
        data[nt] = (float) tdata[nt].Temp;
      fwrite(data, sizeof(float), NTracer, fd);

      for(nt = 0; nt < NTracer; nt++)
        data[nt] = (float) tdata[nt].Utherm;
      fwrite(data, sizeof(float), NTracer, fd);

#ifdef TRACER_TRAJECTORY_EXTENDED_OUTPUT
      for(nt = 0; nt < NTracer; nt++)
        data[nt] = (float) tdata[nt].Vel[0];
      fwrite(data, sizeof(float), NTracer, fd);

      for(nt = 0; nt < NTracer; nt++)
        data[nt] = (float) tdata[nt].Vel[1];
      fwrite(data, sizeof(float), NTracer, fd);

      for(nt = 0; nt < NTracer; nt++)
        data[nt] = (float) tdata[nt].Vel[2];
      fwrite(data, sizeof(float), NTracer, fd);

      for(nt = 0; nt < NTracer; nt++)
        data[nt] = (float) tdata[nt].Dedt;
      fwrite(data, sizeof(float), NTracer, fd);

      int k;
      for(k = 0; k < EOS_NSPECIES; k++)
        {
          for(nt = 0; nt < NTracer; nt++)
            data[nt] = (float) tdata[nt].Composition[k];
          fwrite(data, sizeof(float), NTracer, fd);
        }
#endif

      myfree(data);

      fwrite(&dummy, sizeof(int), 1, fd);

      fclose(fd);
    }

  myfree(tdata);
}

#endif /* TRACER_TRAJECTORY */
