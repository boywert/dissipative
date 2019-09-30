/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/riemann.c
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

#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>

#include "allvars.h"
#include "proto.h"
#include "voronoi.h"

#if !(defined(RIEMANN_HLLC) || defined(RIEMANN_ROSUNOV) || defined(RIEMANN_HLL) || defined(RIEMANN_HLLD) || defined(RIEMANN_GAMMA))

#define GAMMA_G1 ((GAMMA-1.0)/(2.0*GAMMA))
#define GAMMA_G2 ((GAMMA+1.0)/(2.0*GAMMA))
#define GAMMA_G3 ((2.0*GAMMA/(GAMMA-1.0)))
#define GAMMA_G4 (2.0/(GAMMA-1.0))
#define GAMMA_G5 (2.0/(GAMMA+1.0))
#define GAMMA_G6 ((GAMMA-1.0)/(GAMMA+1.0))
#define GAMMA_G7 (0.5*(GAMMA-1.0))
#define GAMMA_G8 (1.0/GAMMA)
#define GAMMA_G9 (GAMMA-1.0)

#define TOL 1.0e-8


/*

int analytic_riemann_solution(int argc, char *argv[])
{
  int i, N;
  double T, *X;
  double Rho_L, Vel_L, Press_L, Rho_R, Vel_R, Press_R, *Rho, *Vel, *Press;
  double Press_M, Vel_M;

  if(argc != 12)
    {
      fprintf(stderr, "\n\nwrong number of arguments ! (found %d) \n\n", argc);
      return (-1);
    }

  N = *(int *) argv[0];
  X = (double *) argv[1];
  T = *(double *) argv[2];

  Rho_L = *(double *) argv[3];
  Vel_L = *(double *) argv[4];
  Press_L = *(double *) argv[5];

  Rho_R = *(double *) argv[6];
  Vel_R = *(double *) argv[7];
  Press_R = *(double *) argv[8];

  Rho = (double *) argv[9];
  Vel = (double *) argv[10];
  Press = (double *) argv[11];


  double Csnd_L = sqrt(GAMMA * Press_L / Rho_L);
  double Csnd_R = sqrt(GAMMA * Press_R / Rho_R);

  riemann(Rho_L, Vel_L, Press_L, Csnd_L, Rho_R, Vel_R, Press_R, Csnd_R, &Press_M, &Vel_M);
  
  for(i = 0; i < N; i++)
    sample_solution(X[i] / T,
		    Rho_L, Vel_L, Press_L, Csnd_L,
		    Rho_R, Vel_R, Press_R, Csnd_R, Press_M, Vel_M, &Rho[i], &Vel[i], &Press[i]);

  return 0;
}
*/


double godunov_flux_3d(struct state *st_L, struct state *st_R, struct state_face *st_face)
{
  double Vel;

#ifndef ISOTHERM_EQS
#ifdef LOCALLY_ISOTHERM_DISK
  if((st_L->inDiskFlag != 1) || (st_R->inDiskFlag != 1))
#endif
    {
      if(st_L->press == 0 && st_R->press == 0)
        {
          /* vacuum state */
          st_face->velx = 0;
          st_face->rho = 0;
          st_face->press = 0;
          st_face->vely = 0;
          st_face->velz = 0;
#ifdef MAXSCALARS
          st_face->scalars = NULL;
#endif
#ifdef TRACER
          st_face->tracer = 0;
#endif
          return 0;
        }

      if(st_L->rho > 0 && st_R->rho > 0)
        {
          st_L->csnd = sqrt(GAMMA * st_L->press / st_L->rho);
          st_R->csnd = sqrt(GAMMA * st_R->press / st_R->rho);

          double Press;

          if(riemann(st_L, st_R, &Press, &Vel))
            {
              sample_solution_3d(0.0,   /* S=x/t */
                                 st_L, st_R, Press, Vel, st_face);
              return Press;
            }
          else
            {
              /* ICs lead to vacuum, need to sample vacuum solution */

              sample_solution_vaccum_generate_3d(0.0,   /* S=x/t */
                                                 st_L, st_R, st_face);
              return 0;
            }
        }
      else
        {
          if(st_L->rho == 0 && st_R->rho > 0)
            {
              sample_solution_vaccum_left_3d(0.0,       /* S=x/t */
                                             st_R, st_face);
              return 0;
            }
          else if(st_R->rho == 0 && st_L->rho > 0)
            {
              sample_solution_vaccum_right_3d(0.0,      /* S=x/t */
                                              st_L, st_face);
              return 0;
            }
          else if(st_R->rho == 0 && st_L->rho == 0)
            {
              /* vacuum state */
              st_face->velx = 0;
              st_face->rho = 0;
              st_face->press = 0;
              st_face->vely = 0;
              st_face->velz = 0;
#ifdef MAXSCALARS
              st_face->scalars = NULL;
#endif
#ifdef TRACER
              st_face->tracer = 0;
#endif
            }
          else
            {
              terminate("one of the densities is negative\n");
            }
          return 0;
        }
    }
#ifdef LOCALLY_ISOTHERM_DISK
  else
    {
      double Rho;
      double csnd;
      csnd = 0.5 * (st_L->localSoundSpeed + st_R->localSoundSpeed);

      riemann_isotherm(st_L, st_R, &Rho, &Vel, csnd);

      sample_solution_isothermal3d(0.0, /* S=x/t */
                                   st_L, st_R, Rho, Vel, st_face, csnd);

      st_face->press = st_face->rho * csnd * csnd;

      return 0;
    }
#endif

#else
  double Rho;
  double csnd;

  csnd = All.IsoSoundSpeed;
  riemann_isotherm(st_L, st_R, &Rho, &Vel, csnd);

  sample_solution_isothermal3d(0.0,     /* S=x/t */
                               st_L, st_R, Rho, Vel, st_face, csnd);

  st_face->press = st_face->rho * csnd * csnd;

  return 0;
#endif
}



void sample_solution_vaccum_left_3d(double S, struct state *st_R, struct state_face *st_face)
{
  double Csnd;

  double Sr = st_R->velx - 2 * st_R->csnd / GAMMA_MINUS1;

  st_face->vely = st_R->vely;
  st_face->velz = st_R->velz;
#ifdef MAXSCALARS
  st_face->scalars = st_R->scalars;
#endif
#ifdef TRACER
  st_face->tracer = st_R->tracer;
#endif

  if(S >= Sr)
    {
      /* right fan */

      double shr = st_R->velx + st_R->csnd;

      if(S >= shr)              /* right data state */
        {
          st_face->rho = st_R->rho;
          st_face->velx = st_R->velx;
          st_face->press = st_R->press;
        }
      else
        {
          /* rarefaction fan right state */
          st_face->velx = GAMMA_G5 * (-st_R->csnd + GAMMA_G7 * st_R->velx + S);
          Csnd = GAMMA_G5 * (st_R->csnd - GAMMA_G7 * (st_R->velx - S));
          st_face->rho = st_R->rho * pow(Csnd / st_R->csnd, GAMMA_G4);
          st_face->press = st_R->press * pow(Csnd / st_R->csnd, GAMMA_G3);
        }
    }
  else
    {
      /* vacuum state */
      st_face->velx = Sr;
      st_face->rho = 0;
      st_face->press = 0;
    }
}

void sample_solution_vaccum_right_3d(double S, struct state *st_L, struct state_face *st_face)
{
  double Csnd;

  double Sl = st_L->velx + 2 * st_L->csnd / GAMMA_MINUS1;

  st_face->vely = st_L->vely;
  st_face->velz = st_L->velz;
#ifdef MAXSCALARS
  st_face->scalars = st_L->scalars;
#endif
#ifdef TRACER
  st_face->tracer = st_L->tracer;
#endif

  if(S <= Sl)
    {
      /* left fan */

      double shl = st_L->velx - st_L->csnd;

      if(S <= shl)
        {
          /* left data state */
          st_face->rho = st_L->rho;
          st_face->velx = st_L->velx;
          st_face->press = st_L->press;
        }
      else
        {
          /* rarefaction fan left state */
          st_face->velx = GAMMA_G5 * (st_L->csnd + GAMMA_G7 * st_L->velx + S);
          Csnd = GAMMA_G5 * (st_L->csnd + GAMMA_G7 * (st_L->velx - S));
          st_face->rho = st_L->rho * pow(Csnd / st_L->csnd, GAMMA_G4);
          st_face->press = st_L->press * pow(Csnd / st_L->csnd, GAMMA_G3);
        }
    }
  else
    {
      /* vacuum in between */
      st_face->velx = Sl;
      st_face->rho = 0;
      st_face->press = 0;
    }
}




void sample_solution_vaccum_generate_3d(double S, struct state *st_L, struct state *st_R, struct state_face *st_face)
{
  double Csnd;

  double Sl = st_L->velx + 2 * st_L->csnd / GAMMA_MINUS1;
  double Sr = st_R->velx - 2 * st_R->csnd / GAMMA_MINUS1;

  if(S <= Sl)
    {
      /* left fan */

      st_face->vely = st_L->vely;
      st_face->velz = st_L->velz;
#ifdef MAXSCALARS
      st_face->scalars = st_L->scalars;
#endif
#ifdef TRACER
      st_face->tracer = st_L->tracer;
#endif

      double shl = st_L->velx - st_L->csnd;

      if(S <= shl)
        {
          /* left data state */
          st_face->rho = st_L->rho;
          st_face->velx = st_L->velx;
          st_face->press = st_L->press;
        }
      else
        {
          /* rarefaction fan left state */
          st_face->velx = GAMMA_G5 * (st_L->csnd + GAMMA_G7 * st_L->velx + S);
          Csnd = GAMMA_G5 * (st_L->csnd + GAMMA_G7 * (st_L->velx - S));
          st_face->rho = st_L->rho * pow(Csnd / st_L->csnd, GAMMA_G4);
          st_face->press = st_L->press * pow(Csnd / st_L->csnd, GAMMA_G3);
        }
    }
  else if(S >= Sr)
    {
      /* right fan */

      double shr = st_R->velx + st_R->csnd;

      st_face->vely = st_R->vely;
      st_face->velz = st_R->velz;
#ifdef MAXSCALARS
      st_face->scalars = st_R->scalars;
#endif
#ifdef TRACER
      st_face->tracer = st_R->tracer;
#endif

      if(S >= shr)              /* right data state */
        {
          st_face->rho = st_R->rho;
          st_face->velx = st_R->velx;
          st_face->press = st_R->press;
        }
      else
        {
          /* rarefaction fan right state */
          st_face->velx = GAMMA_G5 * (-st_R->csnd + GAMMA_G7 * st_R->velx + S);
          Csnd = GAMMA_G5 * (st_R->csnd - GAMMA_G7 * (st_R->velx - S));
          st_face->rho = st_R->rho * pow(Csnd / st_R->csnd, GAMMA_G4);
          st_face->press = st_R->press * pow(Csnd / st_R->csnd, GAMMA_G3);
        }
    }
  else
    {
      /* vacuum in between */
      st_face->velx = S;
      st_face->rho = 0;
      st_face->press = 0;

      st_face->vely = st_L->vely + (st_R->vely - st_L->vely) * (S - Sl) / (Sr - Sl);
      st_face->velz = st_L->velz + (st_R->velz - st_L->velz) * (S - Sl) / (Sr - Sl);

#ifdef MAXSCALARS
      st_face->scalars = NULL;
#endif
#ifdef TRACER
      st_face->tracer = 0;
#endif
    }
}


void get_mach_numbers(struct state *st_L, struct state *st_R, double Press)
{
#ifdef TGCHEM
  double gamma_g1, gamma_g2;
#endif

#if (defined ENTROPY_MACH_THRESHOLD) || (defined GODUNOV_STATS) || (TRACER_MC_MACHMAX) || (TRACER_PART_MACHMAX)
  if(Press <= st_L->press)      /* left fan */
    {
      st_L->mach = 0;
    }
  else                          /* left shock */
    {
      double pml = Press / st_L->press;
#ifdef TGCHEM
      gamma_g1 = (st_L->gamma - 1) / 2 / st_L->gamma;
      gamma_g2 = (st_L->gamma + 1) / 2 / st_L->gamma;

      st_L->mach = sqrt(gamma_g2 * pml + gamma_g1);
#else
      st_L->mach = sqrt(GAMMA_G2 * pml + GAMMA_G1);
#endif
    }

  if(Press > st_R->press)       /* right shock */
    {
      double pmr = Press / st_R->press;
#ifdef TGCHEM
      gamma_g1 = (st_R->gamma - 1) / 2 / st_R->gamma;
      gamma_g2 = (st_R->gamma + 1) / 2 / st_R->gamma;

      st_R->mach = sqrt(gamma_g2 * pmr + gamma_g1);
#else
      st_R->mach = sqrt(GAMMA_G2 * pmr + GAMMA_G1);
#endif
    }
  else
    {
      st_R->mach = 0;
    }
#endif
}




void sample_solution_3d(double S, struct state *st_L, struct state *st_R, double Press, double Vel, struct state_face *st_face)
{
  double Csnd;

  if(S <= Vel)                  /* sample point is left of contact */
    {
      st_face->vely = st_L->vely;
      st_face->velz = st_L->velz;
#ifdef MAXSCALARS
      st_face->scalars = st_L->scalars;
#endif
#ifdef TRACER
      st_face->tracer = st_L->tracer;
#endif

      if(Press <= st_L->press)  /* left fan */
        {
          double shl = st_L->velx - st_L->csnd;

          if(S <= shl)          /* left data state */
            {
              st_face->rho = st_L->rho;
              st_face->velx = st_L->velx;
              st_face->press = st_L->press;
            }
          else
            {
              double cml = st_L->csnd * pow(Press / st_L->press, GAMMA_G1);
              double stl = Vel - cml;

              if(S > stl)       /* middle left state */
                {
                  st_face->rho = st_L->rho * pow(Press / st_L->press, GAMMA_G8);
                  st_face->velx = Vel;
                  st_face->press = Press;
                }
              else              /* left state inside fan */
                {
                  st_face->velx = GAMMA_G5 * (st_L->csnd + GAMMA_G7 * st_L->velx + S);
                  Csnd = GAMMA_G5 * (st_L->csnd + GAMMA_G7 * (st_L->velx - S));
                  st_face->rho = st_L->rho * pow(Csnd / st_L->csnd, GAMMA_G4);
                  st_face->press = st_L->press * pow(Csnd / st_L->csnd, GAMMA_G3);
                }
            }
        }
      else                      /* left shock */
        {
          if(st_L->press > 0)
            {

              double pml = Press / st_L->press;
              double sl = st_L->velx - st_L->csnd * sqrt(GAMMA_G2 * pml + GAMMA_G1);

              if(S <= sl)       /* left data state */
                {
                  st_face->rho = st_L->rho;
                  st_face->velx = st_L->velx;
                  st_face->press = st_L->press;
                }
              else              /* middle left state behind shock */
                {
                  st_face->rho = st_L->rho * (pml + GAMMA_G6) / (pml * GAMMA_G6 + 1.0);
                  st_face->velx = Vel;
                  st_face->press = Press;
                }
            }
          else
            {
              st_face->rho = st_L->rho / GAMMA_G6;
              st_face->velx = Vel;
              st_face->press = Press;
            }
        }
    }
  else                          /* right of contact */
    {
      st_face->vely = st_R->vely;
      st_face->velz = st_R->velz;
#ifdef MAXSCALARS
      st_face->scalars = st_R->scalars;
#endif
#ifdef TRACER
      st_face->tracer = st_R->tracer;
#endif

      if(Press > st_R->press)   /* right shock */
        {
          if(st_R->press > 0)
            {

              double pmr = Press / st_R->press;
              double sr = st_R->velx + st_R->csnd * sqrt(GAMMA_G2 * pmr + GAMMA_G1);

              if(S >= sr)       /* right data state */
                {
                  st_face->rho = st_R->rho;
                  st_face->velx = st_R->velx;
                  st_face->press = st_R->press;
                }
              else              /* middle right state behind shock */
                {
                  st_face->rho = st_R->rho * (pmr + GAMMA_G6) / (pmr * GAMMA_G6 + 1.0);
                  st_face->velx = Vel;
                  st_face->press = Press;
                }
            }
          else
            {
              st_face->rho = st_R->rho / GAMMA_G6;
              st_face->velx = Vel;
              st_face->press = Press;
            }
        }
      else                      /* right fan */
        {
          double shr = st_R->velx + st_R->csnd;

          if(S >= shr)          /* right data state */
            {
              st_face->rho = st_R->rho;
              st_face->velx = st_R->velx;
              st_face->press = st_R->press;
            }
          else
            {
              double cmr = st_R->csnd * pow(Press / st_R->press, GAMMA_G1);
              double str = Vel + cmr;

              if(S <= str)      /* middle right state */
                {
                  st_face->rho = st_R->rho * pow(Press / st_R->press, GAMMA_G8);
                  st_face->velx = Vel;
                  st_face->press = Press;
                }
              else              /* fan right state */
                {
                  st_face->velx = GAMMA_G5 * (-st_R->csnd + GAMMA_G7 * st_R->velx + S);
                  Csnd = GAMMA_G5 * (st_R->csnd - GAMMA_G7 * (st_R->velx - S));
                  st_face->rho = st_R->rho * pow(Csnd / st_R->csnd, GAMMA_G4);
                  st_face->press = st_R->press * pow(Csnd / st_R->csnd, GAMMA_G3);
                }
            }
        }
    }
}

int riemann(struct state *st_L, struct state *st_R, double *Press, double *Vel)
{
  double F_L, FD_L, F_R, FD_R, pold;

  double dVel = st_R->velx - st_L->velx;

  double critVel = GAMMA_G4 * (st_L->csnd + st_R->csnd) - dVel;

  if(critVel < 0)
    {
      /*
         printf("ICs lead to vacuum. stopping. Csnd_L=%g Csnd_R=%g dVel=%g\n", Csnd_L, Csnd_R, dVel);
       */
      return 0;
    }

  double p = guess_for_pressure(st_L, st_R);

  int iter = 0;

  do                            /* newton-raphson scheme */
    {
      pold = p;

      pressure_function(p, st_L, &F_L, &FD_L);
      pressure_function(p, st_R, &F_R, &FD_R);

      if(iter < MAXITER / 2)
        p -= (F_L + F_R + dVel) / (FD_L + FD_R);
      else
        p -= 0.5 * (F_L + F_R + dVel) / (FD_L + FD_R);

      if(p < 0.1 * pold)
        p = 0.1 * pold;

      pressure_function(p, st_L, &F_L, &FD_L);
      pressure_function(p, st_R, &F_R, &FD_R);

      if(iter < MAXITER / 2)
        p -= (F_L + F_R + dVel) / (FD_L + FD_R);
      else
        p -= 0.5 * (F_L + F_R + dVel) / (FD_L + FD_R);

      if(p < 0.1 * pold)
        p = 0.1 * pold;

      iter++;
    }
  while(2 * fabs((p - pold) / (p + pold)) > TOL && iter < MAXITER);

  if(iter >= MAXITER)
    {
      printf("Task=%d: Warning: ICs for riemann solver lead to divergence.\n", ThisTask);
      printf("Rho_L=%g  Vel_L=%g  Press_L=%g  Csnd_L=%g\n", st_L->rho, st_L->velx, st_L->press, st_L->csnd);
      printf("Rho_R=%g  Vel_R=%g  Press_R=%g  Csnd_R=%g\n", st_R->rho, st_R->velx, st_R->press, st_R->csnd);
      printf("Adopted solution: Press=%g Vel=%g\n", p, 0.5 * (st_L->velx + st_R->velx + F_R - F_L));

      FILE *fd;

      if((fd = fopen("riemann.dat", "w")))
        {
          fwrite(&st_L->rho, sizeof(double), 1, fd);
          fwrite(&st_L->velx, sizeof(double), 1, fd);
          fwrite(&st_L->press, sizeof(double), 1, fd);
          fwrite(&st_L->csnd, sizeof(double), 1, fd);
          fwrite(&st_R->rho, sizeof(double), 1, fd);
          fwrite(&st_R->velx, sizeof(double), 1, fd);
          fwrite(&st_R->press, sizeof(double), 1, fd);
          fwrite(&st_R->csnd, sizeof(double), 1, fd);
          fclose(fd);
        }

      // terminate("failed convergence in Riemann solver");
    }

  /* prepare output values */
  *Press = p;
  *Vel = 0.5 * (st_L->velx + st_R->velx + F_R - F_L);

  return 1;
}

void pressure_function(double P, struct state *st, double *F, double *FD)
{
  if(P <= st->press)            /* rarefaction wave */
    {
      double prat = P / st->press;

      *F = GAMMA_G4 * st->csnd * (pow(prat, GAMMA_G1) - 1.0);
      *FD = (1.0 / (st->rho * st->csnd)) * pow(prat, -GAMMA_G2);
    }
  else                          /* shock wave */
    {
      double ak = GAMMA_G5 / st->rho;
      double bk = GAMMA_G6 * st->press;
      double qrt = sqrt(ak / (bk + P));

      *F = (P - st->press) * qrt;
      *FD = (1.0 - 0.5 * (P - st->press) / (bk + P)) * qrt;
    }
}

double guess_for_pressure(struct state *st_L, struct state *st_R)
{
#define QMAX 2.0

  double pmin, pmax;

  double pv = 0.5 * (st_L->press + st_R->press) - 0.125 * (st_R->velx - st_L->velx) * (st_L->rho + st_R->rho) * (st_L->csnd + st_R->csnd);

  if(st_L->press < st_R->press)
    {
      pmin = st_L->press;
      pmax = st_R->press;
    }
  else
    {
      pmin = st_R->press;
      pmax = st_L->press;
    }

  if(pmin > 0)
    {
      double qrat = pmax / pmin;

      if(qrat <= QMAX && (pmin <= pv && pv <= pmax))
        {
          if(pv < 0)
            {
              printf("pv=%g\n", pv);
              terminate("negative pv");
            }

          return pv;
        }
      else
        {
          if(pv < pmin)         /* use two-rarefaction solution */
            {
              double pnu = (st_L->csnd + st_R->csnd) - GAMMA_G7 * (st_R->velx - st_L->velx);
              double pde = st_L->csnd / pow(st_L->press, GAMMA_G1) + st_R->csnd / pow(st_R->press, GAMMA_G1);

              return pow(pnu / pde, GAMMA_G3);
            }
          else                  /* two-shock approximation  */
            {
              double gel = sqrt((GAMMA_G5 / st_L->rho) / (GAMMA_G6 * st_L->press + pv));
              double ger = sqrt((GAMMA_G5 / st_R->rho) / (GAMMA_G6 * st_R->press + pv));
              double x = (gel * st_L->press + ger * st_R->press - (st_R->velx - st_L->velx)) / (gel + ger);

              if(x < pmin || x > pmax)
                {
                  x = pmin;
                }

              return x;
            }
        }
    }
  else
    {

      return (pmin + pmax) / 2;

    }

}



void riemann_isotherm(struct state *st_L, struct state *st_R, double *Rho, double *Vel, double csnd)
{
  double F_L, FD_L, F_R, FD_R, rhoold, drho;
  double rho;

  double dVel = (st_R->velx - st_L->velx) / csnd;

  if(dVel > 0)
    rho = sqrt(st_L->rho * st_R->rho * exp(-dVel));
  else
    rho = 0.5 * (st_L->rho + st_R->rho);

  int iter = 0;

  if(st_L->rho <= 0 || st_R->rho <= 0)
    terminate("isothermal Riemann solver was called with zero or negative density\n");

  do                            /* newton-raphson scheme */
    {
      isothermal_function(rho, st_L->rho, &F_L, &FD_L);
      isothermal_function(rho, st_R->rho, &F_R, &FD_R);

      rhoold = rho;
      drho = -0.5 * (F_L + F_R + dVel) / (FD_L + FD_R);

      if(fabs(drho) > 0.25 * rho)
        drho = 0.25 * rho * fabs(drho) / drho;

      rho += drho;

      iter++;
    }
  while(2 * fabs(rho - rhoold) / (st_L->rho + st_R->rho) > TOL && iter < MAXITER);

  if(iter >= MAXITER)
    {
#ifndef LONGIDS
      printf("ID_L=%u ID_R=%u  Rho_L=%g Rho_R=%g  Vel_L=%g Vel_R=%g\n", st_L->ID, st_R->ID, st_L->rho, st_R->rho, st_L->velx, st_R->velx);
#else
      printf("ID_L=%llu ID_R=%llu  Rho_L=%g Rho_R=%g  Vel_L=%g Vel_R=%g\n", st_L->ID, st_R->ID, st_L->rho, st_R->rho, st_L->velx, st_R->velx);
#endif
      terminate("ICs for isothermal riemann solver lead to divergence. stopping.");
      /*
       *Rho = 0.5 * (Rho_L + Rho_R);
       *Vel = 0.5 * (Vel_L + Vel_R);
       return;
       */
    }

  /* prepare output values */
  *Rho = rho;
  *Vel = 0.5 * (st_L->velx + st_R->velx + csnd * (F_R - F_L));
}

void isothermal_function(double rhostar, double rho, double *F, double *FD)
{
  if(rhostar <= rho)            /* rarefaction wave */
    {
      *F = log(rhostar / rho);
      *FD = 1.0 / rho;
    }
  else                          /* shock wave */
    {
      *F = (rhostar - rho) / sqrt(rhostar * rho);
      *FD = 0.5 / rhostar * (sqrt(rhostar / rho) + sqrt(rho / rhostar));
    }
}



void sample_solution_isothermal3d(double S, struct state *st_L, struct state *st_R, double Rho, double Vel, struct state_face *st_face, double csnd)
{
  if(S <= Vel)                  /* sample point is left of contact */
    {
      st_face->vely = st_L->vely;
      st_face->velz = st_L->velz;
#ifdef MAXSCALARS
      st_face->scalars = st_L->scalars;
#endif
#ifdef TRACER
      st_face->tracer = st_L->tracer;
#endif

      if(Rho <= st_L->rho)      /* left fan */
        {
          double shl = st_L->velx - csnd;

          if(S <= shl)          /* left data state */
            {
              st_face->rho = st_L->rho;
              st_face->velx = st_L->velx;
            }
          else
            {
              double stl = Vel - csnd;

              if(S > stl)       /* middle left state */
                {
                  st_face->rho = Rho;
                  st_face->velx = Vel;
                }
              else              /* left state inside fan */
                {
                  st_face->velx = S + csnd;
                  st_face->rho = st_L->rho * exp(-((S + csnd) - st_L->velx) / csnd);
                }
            }
        }
      else                      /* left shock */
        {
          double sl = (st_L->rho * st_L->velx - Rho * Vel) / (st_L->rho - Rho);

          if(S <= sl)           /* left data state */
            {
              st_face->rho = st_L->rho;
              st_face->velx = st_L->velx;
            }
          else                  /* left state behind shock */
            {
              st_face->rho = Rho;
              st_face->velx = Vel;
            }
        }
    }
  else                          /* right of contact */
    {
      st_face->vely = st_R->vely;
      st_face->velz = st_R->velz;
#ifdef MAXSCALARS
      st_face->scalars = st_R->scalars;
#endif
#ifdef TRACER
      st_face->tracer = st_R->tracer;
#endif

      if(Rho > st_R->rho)       /* right shock */
        {
          double sr = (st_R->rho * st_R->velx - Rho * Vel) / (st_R->rho - Rho);

          if(S >= sr)           /* right data state */
            {
              st_face->rho = st_R->rho;
              st_face->velx = st_R->velx;
            }
          else                  /* right state behind shock */
            {
              st_face->rho = Rho;
              st_face->velx = Vel;
            }
        }
      else                      /* right fan */
        {
          double shr = st_R->velx + csnd;

          if(S >= shr)          /* right data state */
            {
              st_face->rho = st_R->rho;
              st_face->velx = st_R->velx;
            }
          else
            {
              double str = Vel + csnd;

              if(S <= str)      /* middle right state */
                {
                  st_face->rho = Rho;
                  st_face->velx = Vel;
                }
              else              /* fan right state */
                {
                  st_face->velx = S - csnd;
                  st_face->rho = st_R->rho * exp(((S - csnd) - st_R->velx) / csnd);
                }
            }
        }
    }
}


#endif
