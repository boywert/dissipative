/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/OTVET/otvet_CGmethod.c
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
#include <gsl/gsl_math.h>

#include "../allvars.h"
#include "../proto.h"
#include "otvet_proto.h"

#ifdef OTVET

#define MAX_ITER 10000
#define ACCURACY 1.0e-2
#define EPSILON 1.0e-5
#define tiny 1e-8

/*structures for radtransfer*/
struct radtransferdata_in
{
  int NodeList[NODELISTLENGTH];
  MyDouble Pos[3];
  MyFloat Hsml;
  MyFloat ET[6];
  double Kappa, Lambda;
  MyFloat Mass, Density;
}
 *OtvetDataIn, *OtvetDataGet;

struct radtransferdata_out
{
  double Out, Sum;
}
 *OtvetDataResult, *OtvetDataOut;

static double *XVec;
static double *QVec, *DVec, *Residue, *Zvec;
static double *Kappa, *Lambda, *Diag, *Diag2;
static double c_light, dt, a3inv, hubble_a;

void otvet_radtransfer(void)
{
  int i, j, iter;
  double alpha_cg, beta, delta_old, delta_new, sum, min_diag, glob_min_diag, max_diag, glob_max_diag;
  double nH;
  double rel, res, maxrel, glob_maxrel;
  double DQ;
#if defined(OTVET_CHECK_PHOTONCOUNT) && defined(OTVET_CHEMISTRY_PS2011)
  double local_ngamma_in, global_ngamma_in, local_ngamma_src, global_ngamma_src, local_ngamma_out, global_ngamma_out;

  local_ngamma_in = 0.;
  global_ngamma_in = 0.;
  local_ngamma_src = 0.;
  global_ngamma_src = 0.;
  local_ngamma_out = 0.;
  global_ngamma_out = 0.;
#endif


  c_light = CLIGHT / All.UnitVelocity_in_cm_per_s;

  /*  the actual time-step we need to do */
  dt = (All.otvet_Radiation_Ti_endstep - All.otvet_Radiation_Ti_begstep) * All.Timebase_interval;

  if(All.ComovingIntegrationOn)
    {
      a3inv = 1 / (All.Time * All.Time * All.Time);
      hubble_a = hubble_function(All.Time);
      /* in comoving case, timestep is dloga at this point. Convert to dt */
      dt /= hubble_a;
    }
  else
    {
      a3inv = hubble_a = 1.0;
    }

  XVec = (double *) mymalloc("XVec", NumGas * sizeof(double));
  QVec = (double *) mymalloc("QVec", NumGas * sizeof(double));
  DVec = (double *) mymalloc("DVec", NumGas * sizeof(double));
  Residue = (double *) mymalloc("Residue", NumGas * sizeof(double));
  Kappa = (double *) mymalloc("Kappa", NumGas * sizeof(double));
  Lambda = (double *) mymalloc("Lambda", NumGas * sizeof(double));
  Diag = (double *) mymalloc("Diag", NumGas * sizeof(double));
  Zvec = (double *) mymalloc("Zvec", NumGas * sizeof(double));
  Diag2 = (double *) mymalloc("Diag2", NumGas * sizeof(double));

  for(i = 0; i < OT_N_BINS; i++)
    {
      /* initialization for the CG method */
      for(j = 0; j < NumGas; j++)
        if(P[j].Type == 0)
          {
            XVec[j] = 0;
            QVec[j] = 0;
            DVec[j] = 0;
            Residue[j] = 0;
            Kappa[j] = 0;
            Lambda[j] = 0;
            Diag[j] = 0;
            Zvec[j] = 0;
            Diag2[j] = 0;

            XVec[j] = SphP[j].n_gamma[i];

#if defined(OTVET_CHECK_PHOTONCOUNT) && defined(OTVET_CHEMISTRY_PS2011)
            local_ngamma_in += (SphP[j].n_gamma[i] / SphP[j].Density);
#endif

            nH = HYDROGEN_MASSFRAC * SphP[j].Density / PROTONMASS * All.UnitMass_in_g / All.HubbleParam;
            Kappa[j] = a3inv * (SphP[j].HI + tiny) * nH * otvet_sigma_HI[i];

#if defined(OTVET_INCLUDE_HE) && defined(OTVET_MULTI_FREQUENCY)
            Kappa[j] += a3inv * ((SphP[j].HeI + tiny) * nH * otvet_sigma_HeI[i] + (SphP[j].HeII + tiny) * nH * otvet_sigma_HeII[i]);
#endif

            if(All.ComovingIntegrationOn)
              Kappa[j] *= All.Time;


#ifdef OTVET_FLUXLIMITER
            /* now calculate flux limiter */

            if(SphP[j].n_gamma[i] > 0)
              {
                double R = sqrt(SphP[j].Grad_ngamma[0][i] * SphP[j].Grad_ngamma[0][i] +
                                SphP[j].Grad_ngamma[1][i] * SphP[j].Grad_ngamma[1][i] + SphP[j].Grad_ngamma[2][i] * SphP[j].Grad_ngamma[2][i]) / (SphP[j].n_gamma[i] * Kappa[j]);

                if(All.ComovingIntegrationOn)
                  R /= All.Time;

#ifndef OTVET_CHANGEFLUXLIMITER
                R *= 0.1;

                Lambda[j] = (1 + R) / (1 + R + R * R);
#else

                Lambda[j] = (2 + R) / (6 + 3 * R + R * R);
                //      Lambda[j] = (2 + R) / (6 + 3 * R + 0.1 * R * R);
                //      Lambda[j] = (2 + R) / (6 + 3 * R + 0.01 * R * R);
#endif
                if(Lambda[j] < 1e-100)
                  Lambda[j] = 0;
              }
            else
              Lambda[j] = 1.0;
#endif

            /* add the source term */
            SphP[j].n_gamma[i] += dt * SphP[j].Je[i] * P[j].Mass;

#if defined(OTVET_CHECK_PHOTONCOUNT) && defined(OTVET_CHEMISTRY_PS2011)
            local_ngamma_src += (dt * SphP[j].Je[i] * P[j].Mass / SphP[j].Density);
#endif

          }

      otvet_matrix_multiply(XVec, Residue, Diag);

      /* Let's take the diagonal matrix elements as Jacobi preconditioner */

      for(j = 0, min_diag = MAX_REAL_NUMBER, max_diag = -MAX_REAL_NUMBER; j < NumGas; j++)
        if(P[j].Type == 0)
          {
            Residue[j] = SphP[j].n_gamma[i] - Residue[j];

            /* note: in principle we would have to substract the w_ii term, but this is always zero */
            if(Diag[j] < min_diag)
              min_diag = Diag[j];
            if(Diag[j] > max_diag)
              max_diag = Diag[j];

            Zvec[j] = Residue[j] / Diag[j];
            DVec[j] = Zvec[j];
          }

      MPI_Allreduce(&min_diag, &glob_min_diag, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&max_diag, &glob_max_diag, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

      delta_new = otvet_vector_multiply(Zvec, Residue);
      delta_old = delta_new;

#ifndef OTVET_SILENT
      if(ThisTask == 0)
        {
          printf("Begin N_BIN %d/%d\n", i + 1, OT_N_BINS);
          printf("\nBegin CG iteration\nmin-diagonal=%g, max-diagonal=%g\n", glob_min_diag, glob_max_diag);
        }
#endif


      /* begin the CG method iteration */
      iter = 0;

      do
        {
          otvet_matrix_multiply(DVec, QVec, Diag2);

          DQ = otvet_vector_multiply(DVec, QVec);
          if(DQ == 0)
            alpha_cg = 0;
          else
            alpha_cg = delta_new / DQ;


          for(j = 0, maxrel = 0; j < NumGas; j++)
            {
              XVec[j] += alpha_cg * DVec[j];
              Residue[j] -= alpha_cg * QVec[j];

              Zvec[j] = Residue[j] / Diag[j];

              rel = fabs(alpha_cg * DVec[j]) / (XVec[j] + 1.0e-10);
              if(rel > maxrel)
                maxrel = rel;
            }

          delta_old = delta_new;
          delta_new = otvet_vector_multiply(Zvec, Residue);

          sum = otvet_vector_sum(XVec);
          res = otvet_vector_sum(Residue);

          if(delta_old)
            beta = delta_new / delta_old;
          else
            beta = 0;

          for(j = 0; j < NumGas; j++)
            DVec[j] = Zvec[j] + beta * DVec[j];

          MPI_Allreduce(&maxrel, &glob_maxrel, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

#ifndef OTVET_SILENT
          if(ThisTask == 0)
            {
              printf("OTVET radtransfer: iter=%3d  |res|/|x|=%12.6g  maxrel=%12.6g  |x|=%12.6g | res|=%12.6g\n", iter, res / sum, glob_maxrel, sum, res);
              fflush(stdout);
            }
#endif
          iter++;

          if(iter >= MAX_ITER)
            terminate("failed to converge radtransfer\n");
        }
      while((res > ACCURACY * sum && iter < MAX_ITER) || iter < 2);

#ifdef OTVET_SILENT
      if(ThisTask == 0)
        {
          printf("%d iterations performed\n", iter);
          fflush(stdout);
        }
#endif

      /* update the intensity */
      for(j = 0; j < NumGas; j++)
        if(P[j].Type == 0)
          {
            if(XVec[j] < 0 || isnan(XVec[j]))
              XVec[j] = 0;

            SphP[j].n_gamma[i] = XVec[j];
#if defined(OTVET_CHECK_PHOTONCOUNT) && defined(OTVET_CHEMISTRY_PS2011)
            local_ngamma_out += (SphP[j].n_gamma[i] / SphP[j].Density);
#endif
          }
    }

  myfree(Diag2);
  myfree(Zvec);
  myfree(Diag);
  myfree(Lambda);
  myfree(Kappa);
  myfree(Residue);
  myfree(DVec);
  myfree(QVec);
  myfree(XVec);

/* check photon conserving -------- */
#if defined(OTVET_CHECK_PHOTONCOUNT) && defined(OTVET_CHEMISTRY_PS2011)

  MPI_Reduce(&local_ngamma_in, &global_ngamma_in, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_ngamma_src, &global_ngamma_src, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_ngamma_out, &global_ngamma_out, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  double expected = global_ngamma_in + global_ngamma_src;
  mpi_printf("OTVET PHOTONCOUNT: in=%g src=%g expected=%g out=%g\n", global_ngamma_in, global_ngamma_src, expected, global_ngamma_out);
#endif
}

/* internal product of two vectors */
double otvet_vector_multiply(double *a, double *b)
{
  int i;
  double sum, sumall;

  for(i = 0, sum = 0; i < NumGas; i++)
    if(P[i].Type == 0)
      sum += a[i] * b[i];

  MPI_Allreduce(&sum, &sumall, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return sumall;
}


double otvet_vector_sum(double *a)
{
  int i;
  double sum, sumall;

  for(i = 0, sum = 0; i < NumGas; i++)
    if(P[i].Type == 0)
      sum += fabs(a[i]);

  MPI_Allreduce(&sum, &sumall, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return sumall;
}


/* this function computes the vector b(out) given the vector x(in) such as Ax = b, where A is a matrix */
void otvet_matrix_multiply(double *in, double *out, double *sum)
{
  int i, j, k, ngrp, dummy, ndone, ndone_flag;
  int recvTask, nexport, nimport, place;
  double ainv, dt;

  /* allocate buffers to arrange communication */

  Ngblist = (int *) mymalloc("Ngblist", NumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(data_index) + sizeof(struct data_nodelist) +
                                             sizeof(struct radtransferdata_in) + sizeof(struct radtransferdata_out) + sizemax(sizeof(struct radtransferdata_in), sizeof(struct radtransferdata_out))));
  DataIndexTable = (data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(data_index));
  DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

  dt = (All.otvet_Radiation_Ti_endstep - All.otvet_Radiation_Ti_begstep) * All.Timebase_interval;

  if(All.ComovingIntegrationOn)
    {
      ainv = 1.0 / All.Time;
      /* in comoving case, timestep is dloga at this point. Convert to dt */
      dt /= hubble_function(All.Time);
    }
  else
    {
      ainv = 1.0;
    }

  i = 0;

  do                            /* communication loop */
    {

      for(j = 0; j < NTask; j++)
        {
          Send_count[j] = 0;
          Exportflag[j] = -1;
        }

      /* do local particles and prepare export list */
      for(nexport = 0; i < NumGas; i++)
        {
          if(P[i].Type == 0)
            if(radtransfer_evaluate(i, 0, in, out, sum, &nexport, Send_count) < 0)
              break;
        }

      mysort_dataindex(DataIndexTable, nexport, sizeof(data_index), data_index_compare);

      MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

      for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
        {
          nimport += Recv_count[j];

          if(j > 0)
            {
              Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
              Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
            }
        }

      OtvetDataGet = (struct radtransferdata_in *) mymalloc("OtvetDataGet", nimport * sizeof(struct radtransferdata_in));
      OtvetDataIn = (struct radtransferdata_in *) mymalloc("OtvetDataIn", nexport * sizeof(struct radtransferdata_in));

      /* prepare particle data for export */

      for(j = 0; j < nexport; j++)
        {
          place = DataIndexTable[j].Index;
          for(k = 0; k < 3; k++)
            {
              OtvetDataIn[j].Pos[k] = P[place].Pos[k];
              OtvetDataIn[j].ET[k] = SphP[place].ET[k];
              OtvetDataIn[j].ET[k + 3] = SphP[place].ET[k + 3];
            }
          OtvetDataIn[j].Hsml = SphP[place].Hsml;
          OtvetDataIn[j].Kappa = Kappa[place];
          OtvetDataIn[j].Lambda = Lambda[place];
          OtvetDataIn[j].Mass = P[place].Mass;
          OtvetDataIn[j].Density = SphP[place].Density;

          memcpy(OtvetDataIn[j].NodeList, DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
        }

      /* exchange particle data */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
          recvTask = ThisTask ^ ngrp;

          if(recvTask < NTask)
            {
              if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                {
                  /* get the particles */
                  MPI_Sendrecv(&OtvetDataIn[Send_offset[recvTask]],
                               Send_count[recvTask] * sizeof(struct radtransferdata_in), MPI_BYTE,
                               recvTask, TAG_OTVET_A,
                               &OtvetDataGet[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct radtransferdata_in), MPI_BYTE, recvTask, TAG_OTVET_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }

      myfree(OtvetDataIn);
      OtvetDataResult = (struct radtransferdata_out *) mymalloc("OtvetDataResult", nimport * sizeof(struct radtransferdata_out));
      OtvetDataOut = (struct radtransferdata_out *) mymalloc("OtvetDataOut", nexport * sizeof(struct radtransferdata_out));

      /* now do the particles that were sent to us */
      for(j = 0; j < nimport; j++)
        radtransfer_evaluate(j, 1, in, out, sum, &dummy, &dummy);

      if(i < NumGas)
        ndone_flag = 0;
      else
        ndone_flag = 1;

      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      /* get the result */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
          recvTask = ThisTask ^ ngrp;
          if(recvTask < NTask)
            {
              if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                {
                  /* send the results */
                  MPI_Sendrecv(&OtvetDataResult[Recv_offset[recvTask]],
                               Recv_count[recvTask] * sizeof(struct radtransferdata_out),
                               MPI_BYTE, recvTask, TAG_OTVET_B,
                               &OtvetDataOut[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct radtransferdata_out), MPI_BYTE, recvTask, TAG_OTVET_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }

      /* add the result to the local particles */
      for(j = 0; j < nexport; j++)
        {
          place = DataIndexTable[j].Index;
          out[place] += OtvetDataOut[j].Out;
          sum[place] += OtvetDataOut[j].Sum;
        }

      myfree(OtvetDataOut);
      myfree(OtvetDataResult);
      myfree(OtvetDataGet);

    }
  while(ndone < NTask);

  /* do final operations on results */
  for(i = 0; i < NumGas; i++)
    if(P[i].Type == 0)
      {
        /* divide by c_light to get comoving speed of light (because kappa is comoving) */
#ifdef OTVET_CHEMISTRY_PS2009
        if((1 + dt * c_light * ainv * Kappa[i] + sum[i]) < 0)
#else
        if((1 + sum[i]) < 0)    /* solve only transport, not absorption */
#endif
          {
            printf("1 + sum + rate= %g   sum=%g rate=%g i =%d\n", 1 + dt * c_light * ainv * Kappa[i] + sum[i], sum[i], dt * c_light * ainv * Kappa[i], i);
            terminate("11111111");
          }

#ifdef OTVET_CHEMISTRY_PS2009
        sum[i] += 1.0 + dt * c_light * ainv * Kappa[i]; /* LVS: Eq. 67 */
#else
        sum[i] += 1.0;          /* solves only transport */
#endif

        out[i] += in[i] * sum[i];
      }

  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);
}

/* this function evaluates parts of the matrix A */
int radtransfer_evaluate(int target, int mode, double *in, double *out, double *sum, int *nexport, int *nsend_local)
{
  int startnode, numngb, listindex = 0;
  int j, n, k;
  MyFloat *ET_aux, ET_j[6], ET_i[6], ET_ij[6];
  MyFloat kappa_i, kappa_j, kappa_ij;

#ifdef OTVET_FLUXLIMITER
  MyFloat lambda_i, lambda_j;
#endif
  MyDouble *pos;
  MyFloat mass, mass_i, rho, rho_i;
  double sum_out = 0, sum_w = 0, fac = 0;

  double dx, dy, dz;
  double h_j, hinv, hinv4, h_i;
  double dwk_i, dwk_j, dwk;
  double r, r2, r3inv, u, ainv, dt;
#if defined(PERIODIC) && !defined(GRAVITY_NOT_PERIODIC)
  double xtmp, ytmp, ztmp;
#endif

  dt = (All.otvet_Radiation_Ti_endstep - All.otvet_Radiation_Ti_begstep) * All.Timebase_interval;

  if(All.ComovingIntegrationOn)
    {
      ainv = 1.0 / All.Time;
      /* in comoving case, timestep is dloga at this point. Convert to dt */
      dt /= hubble_function(All.Time);
    }
  else
    {
      ainv = 1.0;
    }

  if(mode == 0)
    {
      ET_aux = SphP[target].ET;
      pos = P[target].Pos;
      h_i = SphP[target].Hsml;
      kappa_i = Kappa[target];
#ifdef OTVET_FLUXLIMITER
      lambda_i = Lambda[target];
#endif
      mass_i = P[target].Mass;
      rho_i = SphP[target].Density;
    }
  else
    {
      ET_aux = OtvetDataGet[target].ET;
      pos = OtvetDataGet[target].Pos;
      h_i = OtvetDataGet[target].Hsml;
      kappa_i = OtvetDataGet[target].Kappa;
#ifdef OTVET_FLUXLIMITER
      lambda_i = OtvetDataGet[target].Lambda;
#endif
      mass_i = OtvetDataGet[target].Mass;
      rho_i = OtvetDataGet[target].Density;
    }

#ifdef OTVET_MODIFY_EDDINGTON_TENSOR
  /*modify Eddington tensor */
  ET_i[0] = 2 * ET_aux[0] - 0.5 * ET_aux[1] - 0.5 * ET_aux[2];
  ET_i[1] = 2 * ET_aux[1] - 0.5 * ET_aux[2] - 0.5 * ET_aux[0];
  ET_i[2] = 2 * ET_aux[2] - 0.5 * ET_aux[0] - 0.5 * ET_aux[1];

  for(k = 3; k < 6; k++)
    ET_i[k] = 2.5 * ET_aux[k];
#else
  for(k = 0; k < 6; k++)
    ET_i[k] = ET_aux[k];
#endif

  if(mode == 0)
    {
      startnode = Ngb_MaxPart;
    }
  else
    {
      startnode = OtvetDataGet[target].NodeList[0];
      startnode = Ngb_Nodes[startnode].u.d.nextnode;
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
        {
          numngb = ngb_treefind_variable(pos, h_i, target, &startnode, mode, nexport, nsend_local);

          if(numngb < 0)
            return -1;

          for(n = 0; n < numngb; n++)
            {
              j = Ngblist[n];
              if(P[j].Mass > 0. && P[j].ID != 0)        /* skip cells that have been swallowed or dissolved */
                {
                  dx = GRAVITY_NEAREST_X(pos[0] - P[j].Pos[0]);
                  dy = GRAVITY_NEAREST_Y(pos[1] - P[j].Pos[1]);
                  dz = GRAVITY_NEAREST_Z(pos[2] - P[j].Pos[2]);

                  r2 = dx * dx + dy * dy + dz * dz;
                  r = sqrt(r2);
                  r3inv = 1.0 / (r2 * r);
                  h_j = SphP[j].Hsml;

                  if(r > 0 && (r < h_i || r < h_j))
                    {
                      mass = P[j].Mass;
                      rho = SphP[j].Density;
                      kappa_j = Kappa[j];
#ifdef OTVET_FLUXLIMITER
                      lambda_j = Lambda[j];
#endif

#ifdef OTVET_MODIFY_EDDINGTON_TENSOR
                      ET_aux = SphP[j].ET;

                      /*modify Eddington tensor */
                      ET_j[0] = 2 * ET_aux[0] - 0.5 * ET_aux[1] - 0.5 * ET_aux[2];
                      ET_j[1] = 2 * ET_aux[1] - 0.5 * ET_aux[2] - 0.5 * ET_aux[0];
                      ET_j[2] = 2 * ET_aux[2] - 0.5 * ET_aux[0] - 0.5 * ET_aux[1];

                      for(k = 3; k < 6; k++)
                        ET_j[k] = 2.5 * ET_aux[k];
#else
                      for(k = 0; k < 6; k++)
                        ET_j[k] = SphP[j].ET[k];
#endif

                      for(k = 0; k < 6; k++)
                        ET_ij[k] = 0.5 * (ET_i[k] + ET_j[k]);

                      if(r < h_i)
                        {
                          hinv = 1.0 / h_i;
                          hinv4 = hinv * hinv * hinv * hinv;
                          u = r * hinv;

                          if(u < 0.5)
                            dwk_i = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
                          else
                            dwk_i = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
                        }
                      else
                        dwk_i = 0;

                      if(r < h_j)
                        {
                          hinv = 1.0 / h_j;
                          hinv4 = hinv * hinv * hinv * hinv;
                          u = r * hinv;

                          if(u < 0.5)
                            dwk_j = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
                          else
                            dwk_j = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
                        }
                      else
                        dwk_j = 0;

                      kappa_ij = 0.5 * (1 / kappa_i + 1 / kappa_j);
                      dwk = 0.5 * (dwk_i + dwk_j);
                      mass = 0.5 * (mass + mass_i);
                      rho = 0.5 * (rho + rho_i);

                      double tensor = (ET_ij[0] * dx * dx + ET_ij[1] * dy * dy + ET_ij[2] * dz * dz + 2.0 * ET_ij[3] * dx * dy + 2.0 * ET_ij[4] * dy * dz + 2.0 * ET_ij[5] * dz * dx);

                      if(tensor > 0)
                        {
                          fac = -2.0 * dt * c_light * ainv * (mass / rho) * kappa_ij * dwk * r3inv * tensor;

#ifdef OTVET_FLUXLIMITER
                          fac *= 0.5 * (lambda_i + lambda_j);
#endif

                          sum_out -= fac * in[j];

                          sum_w += fac;
                        }
                    }
                }
            }
        }

      if(mode == 1)
        {
          listindex++;
          if(listindex < NODELISTLENGTH)
            {
              startnode = OtvetDataGet[target].NodeList[listindex];
              if(startnode >= 0)
                startnode = Ngb_Nodes[startnode].u.d.nextnode;
            }
        }
    }

  if(mode == 0)
    {
      out[target] = sum_out;
      sum[target] = sum_w;
    }
  else
    {
      OtvetDataResult[target].Out = sum_out;
      OtvetDataResult[target].Sum = sum_w;
    }

  return 0;
}

/* this function sets up simple initial conditions for a single source in a uniform field of gas with constant density*/
void otvet_set_simple_inits(void)
{
  int i, j;

  for(i = 0; i < NumGas; i++)
    if(P[i].Type == 0)
      {
        for(j = 0; j < OT_N_BINS; j++)
          SphP[i].n_gamma[j] = 0.0;

        /* in code units */
        SphP[i].HII = tiny;
        SphP[i].HI = 1.0 - SphP[i].HII;
        SphP[i].Ne = SphP[i].HII;

#ifdef OTVET_INCLUDE_HE
        double fac = (1 - HYDROGEN_MASSFRAC) / 4.0 / HYDROGEN_MASSFRAC;

        SphP[i].HeIII = tiny * fac;
        SphP[i].HeII = tiny * fac;
        SphP[i].HeI = (1.0 - SphP[i].HeII - SphP[i].HeIII) * fac;

        SphP[i].Ne += SphP[i].HeII + 2.0 * SphP[i].HeIII;
#endif
      }
}

void ot_get_sigma(void)
{
  double fac = 1.0 / All.UnitLength_in_cm / All.UnitLength_in_cm * All.HubbleParam * All.HubbleParam;

#ifndef OTVET_MULTI_FREQUENCY
  otvet_sigma_HI[0] = 6.3e-18 * fac;

#else
  int i, j, integral;
  double e, d_nu, e_start, e_end;
  double sum_HI_sigma, sum_HI_G;
  double hc, T_eff, I_nu;
  double sig, f, fac_two;
#ifdef OTVET_INCLUDE_HE
  double sum_HeI_sigma, sum_HeII_sigma;
  double sum_HeI_G, sum_HeII_G;
#endif

  T_eff = All.star_Teff;
  hc = CLIGHT * PLANCK;

  integral = 10000;

  fac_two = ELECTRONVOLT_IN_ERGS / All.UnitEnergy_in_cgs * All.HubbleParam;

  nu[0] = 13.6;
  nu[1] = 24.6;
  nu[2] = 54.4;
  nu[3] = 70.0;

  sum_HI_sigma = 0.0;
  sum_HI_G = 0.0;
#ifdef OTVET_INCLUDE_HE
  sum_HeI_G = sum_HeII_G = 0.0;
  sum_HeI_sigma = 0.0;
  sum_HeII_sigma = 0.0;
#endif

  for(i = 0; i < OT_N_BINS; i++)
    {
      e_start = nu[i];

      if(i == OT_N_BINS - 1)
        e_end = 500.0;
      else
        e_end = nu[i + 1];

      d_nu = (e_end - e_start) / (float) (integral - 1);

      otvet_sigma_HI[i] = 0.0;
      G_HI[i] = 0.0;

#ifdef OTVET_INCLUDE_HE
      otvet_sigma_HeI[i] = 0.0;
      otvet_sigma_HeII[i] = 0.0;
      G_HeI[i] = G_HeII[i] = 0.0;
#endif

      for(j = 0; j < integral; j++)
        {
          e = e_start + j * d_nu;

          I_nu = 2.0 * pow(e * ELECTRONVOLT_IN_ERGS, 3) / (hc * hc) / (exp(e * ELECTRONVOLT_IN_ERGS / (BOLTZMANN * T_eff)) - 1.0);

          if(nu[i] >= 13.6)
            {
              f = sqrt((e / 13.6) - 1.0);

              if(j == 0)
                sig = 6.3e-18;
              else
                sig = 6.3e-18 * pow(13.6 / e, 4) * exp(4 - (4 * atan(f) / f)) / (1.0 - exp(-2 * M_PI / f));

              otvet_sigma_HI[i] += d_nu * sig * I_nu / e;

              sum_HI_sigma += d_nu * I_nu / e;

              G_HI[i] += d_nu * sig * (e - 13.6) * I_nu / e;

              sum_HI_G += d_nu * sig * I_nu / e;
            }

#ifdef OTVET_INCLUDE_HE
          if(nu[i] >= 24.6)
            {
              f = sqrt((e / 24.6) - 1.0);

              if(j == 0)
                sig = 7.83e-18;
              else
                sig = 7.83e-18 * pow(24.6 / e, 4) * exp(4 - (4 * atan(f) / f)) / (1.0 - exp(-2 * M_PI / f));

              otvet_sigma_HeI[i] += d_nu * sig * I_nu / e;

              sum_HeI_sigma += d_nu * I_nu / e;

              G_HeI[i] += d_nu * sig * (e - 24.6) * I_nu / e;

              sum_HeI_G += d_nu * sig * I_nu / e;
            }

          if(nu[i] >= 54.4)
            {
              f = sqrt((e / 54.4) - 1.0);

              if(j == 0)
                sig = 1.58e-18;
              else
                sig = 1.58e-18 * pow(54.4 / e, 4) * exp(4 - (4 * atan(f) / f)) / (1.0 - exp(-2 * M_PI / f));

              otvet_sigma_HeII[i] += d_nu * sig * I_nu / e;

              sum_HeII_sigma += d_nu * I_nu / e;

              G_HeII[i] += d_nu * sig * (e - 54.4) * I_nu / e;

              sum_HeII_G += d_nu * sig * I_nu / e;
            }
#endif
        }
    }
  for(i = 0; i < OT_N_BINS; i++)
    {
      if(nu[i] >= 13.6)
        {
          otvet_sigma_HI[i] *= fac / sum_HI_sigma;
          G_HI[i] *= fac_two / sum_HI_G;
        }

#ifdef OTVET_INCLUDE_HE
      if(nu[i] >= 24.6)
        {
          otvet_sigma_HeI[i] *= fac / sum_HeI_sigma;
          G_HeI[i] *= fac_two / sum_HeI_G;
        }

      if(nu[i] >= 54.4)
        {
          otvet_sigma_HeII[i] *= fac / sum_HeII_sigma;
          G_HeII[i] *= fac_two / sum_HeII_G;
        }
#endif
    }

  if(ThisTask == 0)
    for(i = 0; i < OT_N_BINS; i++)
      printf("OTVET SIGMA: %g %g | %g %g | %g %g\n", otvet_sigma_HI[i] / fac, G_HI[i] / fac_two, otvet_sigma_HeI[i] / fac, G_HeI[i] / fac_two, otvet_sigma_HeII[i] / fac, G_HeII[i] / fac_two);

#endif
}

#ifdef OTVET_MULTI_FREQUENCY
#if defined(EDDINGTON_TENSOR_STARS) || defined(EDDINGTON_TENSOR_SFR)
void ot_get_lum_stars(void)
{
  int i;
  double T_eff, R_eff, I_nu;
  double hc, sum, d_nu;
  int j, integral;
  double e, e_start, e_end;
  integral = 10000;

  T_eff = All.star_Teff;
  R_eff = 7e11;
  hc = CLIGHT * PLANCK;

  for(i = 0, sum = 0; i < OT_N_BINS; i++)
    {
      e_start = nu[i];

      if(i == OT_N_BINS - 1)
        e_end = 500.0;
      else
        e_end = nu[i + 1];

      d_nu = (e_end - e_start) / (float) (integral - 1);

      lum[i] = 0.0;

      for(j = 0; j < integral; j++)
        {
          e = e_start + j * d_nu;

          I_nu = 2.0 * pow(e * ELECTRONVOLT_IN_ERGS, 3) / (hc * hc) / (exp(e * ELECTRONVOLT_IN_ERGS / (BOLTZMANN * T_eff)) - 1.0);

          lum[i] += 4.0 * M_PI * R_eff * R_eff * M_PI * I_nu / e * d_nu / PLANCK;       // number/s
        }

      sum += lum[i];

      lum[i] *= All.UnitTime_in_s / All.HubbleParam;    //number/time
    }

#ifdef OTVET_TEST_SST
#ifdef EDDINGTON_TENSOR_STARS
  for(i = 0; i < OT_N_BINS; i++)
    lum[i] *= All.IonizingLumPerSolarMass / sum;
#endif
#ifdef EDDINGTON_TENSOR_SFR
  for(i = 0; i < OT_N_BINS; i++)
    lum[i] *= All.IonizingLumPerSFR / sum;
#endif
#endif


  if(ThisTask == 0)
    {
      fprintf(FdOTVETStar, "OTVET T_eff %g\n", All.star_Teff);
      fflush(FdOTVETStar);
      for(i = 0; i < OT_N_BINS; i++)
        {
          fprintf(FdOTVETStar, "%d %g %g\n", i, nu[i], lum[i] / All.UnitTime_in_s * All.HubbleParam);
          fflush(FdOTVETStar);
        }
    }

}
#endif
#endif

#endif
