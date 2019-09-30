/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/sidm/sidm_cross.c
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

#ifdef SIDM

static int CrossVbins;
static double CrossUnitFac;
static double Dvlog;
static double **CrossTable;
static double *VelTable;

void sidm_Init_CrossSection(void)
{
  mpi_printf("SIDM: Init (reading cross section files, scatter matrix file, scatter state file...)\n");

  CrossUnitFac = All.UnitMass_in_g / (All.UnitLength_in_cm * All.UnitLength_in_cm);

  mpi_printf("SIDM: CrossUnitFac = %g\n", CrossUnitFac);

  load_scatter_matrix_and_state();

#if !defined(SIDM_CONST_CROSS) && !defined(SIDM_MAXWELLIAN)
  init_cross_table();
#endif

  mpi_printf("SIDM: done.\n");
}


/* load scatter matrix */
void load_scatter_matrix_and_state(void)
{
  FILE *fpin;
  unsigned char state, in_state1, in_state2, out_state1, out_state2;
  unsigned char reaction;
  float delta_mass, initial_fraction;
  float sum_initial_fraction;
  char buf[255], ch;
  int lines;

  //scatter matrix format (each line is a rection with individual cross section)
  //<state of inparticle 1><state of inparticle 2><state of outparticle 1> <state of outparticle 2> 
  //
  //the cross section/mass for each reaction has to be tabulated in
  //a file called sidm_cross_reaction_<num reaction=line number in scatter file>.txt

  sprintf(buf, "scatter_matrix.txt");
  fpin = fopen(buf, "r");

  if(fpin == NULL)
    terminate("SIDM: scatter matrix file 'scatter_matrix.txt' not found.");

  lines = 0;
  while(!feof(fpin))
    {
      ch = fgetc(fpin);
      if(ch == '\n')
        lines++;
    }

  if(SIDM_REACTIONS != lines)
    terminate("SIDM: scatter matrix file 'scatter_matrix.txt' wrong: incorrect number of reactions  lines=%d  SIDM_REACTIONS=%d\n", lines, SIDM_REACTIONS);

  mpi_printf("SIDM: number of reactions in scatter matrix = %d\n", SIDM_REACTIONS);
  fclose(fpin);

  fpin = fopen(buf, "r");
  for(reaction = 0; reaction < SIDM_REACTIONS; reaction++)
    {
      fscanf(fpin, "%hhu %hhu %hhu %hhu", &in_state1, &in_state2, &out_state1, &out_state2);

      SMSIDM[reaction].Out1 = out_state1;
      SMSIDM[reaction].Out2 = out_state2;
      SMSIDM[reaction].In1 = in_state1;
      SMSIDM[reaction].In2 = in_state2;

      if((in_state1 >= SIDM_STATES) || (in_state2 >= SIDM_STATES) || (out_state1 >= SIDM_STATES) || (out_state2 >= SIDM_STATES))
        terminate("SIDM: scatter matrix file 'scatter_matrix.txt' wrong: incorrect states  reaction=%d  in_state1=%d  in_state2=%d  out_state1=%d  out_state2=%d\n", reaction, in_state1, in_state2, out_state1, out_state2);

      mpi_printf("SIDM: Scatter Matrix: reaction= %02d:   %02d  %02d -->  %02d  %02d\n", reaction, in_state1, in_state2, out_state1, out_state2);
    }
  fclose(fpin);


  //scatter states


  sprintf(buf, "scatter_states.txt");
  fpin = fopen(buf, "r");

  if(fpin == NULL)
    terminate("SIDM: scatter states file 'scatter_states.txt' not found.");

  lines = 0;
  while(!feof(fpin))
    {
      ch = fgetc(fpin);
      if(ch == '\n')
        lines++;
    }

  if(SIDM_STATES != lines)
    terminate("SIDM: scatter states file 'scatter_states.txt' wrong: incorrect number of states\n");

  mpi_printf("SIDM: number of states in scatter states = %d\n", SIDM_STATES);
  fclose(fpin);

  fpin = fopen(buf, "r");
  for(sum_initial_fraction = 0.0, state = 0; state < SIDM_STATES; state++)
    {
      fscanf(fpin, "%g %g", &delta_mass, &initial_fraction);

      STSIDM[state].DeltaMass = delta_mass;
      STSIDM[state].InitialFraction = initial_fraction;

      sum_initial_fraction += initial_fraction;

      if ((state == 0) && (delta_mass != 0.0))
       terminate("SIDM: scatter states file 'scatter_states.txt' wrong: ground state delta_mass = %g (expected 0.0 since the first state should be the ground state)\n", delta_mass);

      mpi_printf("SIDM: Scatter State: state=%d:  delta_mass=%g\n", state, delta_mass);
    }
  if (fabs(sum_initial_fraction-1.0) > 1e-5)
   terminate("SIDM: sum_initial_fraction-1.0=%g (should be 0.0)\n", sum_initial_fraction-1.0);
  fclose(fpin);


}




/* evaluate cross section */
#if defined(SIDM_CONST_CROSS) 
MyDouble sidm_cross_sigma(MyDouble rel_vel, unsigned char reaction)
{
  return CrossUnitFac * All.CrossSectionPerMass_in_cgs;
}
#elif defined(SIDM_MAXWELLIAN)
MyDouble sidm_cross_sigma(MyDouble rel_vel, unsigned char reaction)
{
  return CrossUnitFac * All.CrossSectionPerMass_in_cgs / ((rel_vel > 0.0)? (rel_vel):(0.0));
}
#else
MyDouble sidm_cross_sigma(MyDouble rel_vel, unsigned char reaction)
{
  int bin = (int) ( (log(rel_vel) - log(VelTable[0])) / Dvlog );

  /* out of velocity range -> set cross section to boundary values */
  /*
  if (rel_vel < VelTable[0]) 
    bin = 0;
  if (rel_vel > VelTable[CrossVbins - 1]) 
    bin = CrossVbins - 1;
  */

  /* out of velocity range -> set cross section to zero */
  if (rel_vel < VelTable[0])
    return 0.0;
  if (rel_vel > VelTable[CrossVbins - 1])
    return 0.0;

  return CrossTable[reaction][bin];
}
#endif



/* load cross section */
void init_cross_table(void)
{
  int i;
  double tmp_cross, tmp_rel_vel;
  char buf[255];
  int lines = 0;
  char ch;
  FILE *fp = NULL, *fpin = NULL;
  unsigned char reaction;

  //the cross section/mass for each reaction has to be tabulated in
  //  //a file called sidm_cross_reaction_<num reaction=line number in scatter file>.txt
  mpi_printf("SIDM: reading cross sections\n");

  sprintf(buf, "sidm_cross_reaction_0.txt");
  fpin = fopen(buf, "r");

  if(fpin == NULL)
    terminate("SIDM: cross section file 'sidm_cross_reaction_0.txt' not found");

  while(!feof(fpin))
    {
      ch = fgetc(fpin);
      if(ch == '\n')
        lines++;
    }

   fclose(fpin);

   CrossVbins = lines;

   mpi_printf("SIDM: read %d lines\n", lines);
   fflush(stdout);



  /* allocate on all tasks */
  VelTable = (double *) mymalloc("VelTable", CrossVbins * sizeof(double *));
  CrossTable = (double **) mymalloc("CrossTable", SIDM_REACTIONS * sizeof(double *));
  for(reaction = 0; reaction < SIDM_REACTIONS; reaction++)
    CrossTable[reaction] = (double *) mymalloc("CrossTable", CrossVbins * sizeof(double));

  {

    for(reaction = 0; reaction < SIDM_REACTIONS; reaction++)
      {
        sprintf(buf, "sidm_cross_reaction_%d.txt", reaction);
        fpin = fopen(buf, "r");

        if(fpin == NULL)
          terminate("SIDM: cross section file '%s' not found", buf);


        fpin = fopen(buf, "r");

        for(i = 0; i < CrossVbins; i++)
          {
            fscanf(fpin, "%lf %lf", &tmp_rel_vel, &tmp_cross);

            VelTable[i] = tmp_rel_vel;  //km/s
            CrossTable[reaction][i] = tmp_cross * CrossUnitFac; //cm^2/g

            if(i > 0)
              if(VelTable[i - 1] >= VelTable[i])
                terminate("SIDM: cross section file '%s' wrong:  velocity bins are not monotonic", buf);

          }

        fclose(fpin);
      }

    Dvlog = log(VelTable[CrossVbins - 1] / VelTable[0]) / CrossVbins;

    if(ThisTask == 0)
      {
        for(reaction = 0; reaction < SIDM_REACTIONS; reaction++)
          {
            sprintf(buf, "%s/sidm_cross_reaction_%d.txt", All.OutputDir, reaction);
            fp = fopen(buf, "w");
            for(i = 0; i < CrossVbins; i++)
              fprintf(fp, "%d %g %g %g %g %g\n", i, VelTable[i], CrossTable[reaction][i], CrossTable[reaction][i] / CrossUnitFac, sidm_cross_sigma(VelTable[i], reaction),
                      sidm_cross_sigma(VelTable[i], reaction) / CrossUnitFac);
            fclose(fp);
          }
      }
    mpi_printf("SIDM: read cross section  -->  CrossVbins=%d   Dvlog=%g   minvel=%g   maxvel=%g\n", CrossVbins, Dvlog, VelTable[0], VelTable[CrossVbins - 1]);
  }
}


double sidm_scatter_P(double phys_rho, double phys_rel_vel, double Ekin, unsigned char reaction, int *retval)
{
  double delta_E_1, delta_E_2, delta_E;

  if (sidm_cross_sigma(phys_rel_vel, reaction) == 0.0)
   {
     *retval = -2;     
     return 0.0; 
   }   

  sidm_get_delta_energy(reaction, &delta_E_1, &delta_E_2);

  delta_E = delta_E_1 + delta_E_2;

  if ( (Ekin + delta_E) >= 0)
   {
     *retval = 0;
     return phys_rho * sidm_cross_sigma(phys_rel_vel, reaction) * phys_rel_vel / 2.0  *  1.0 / (1.0 + STSIDM[SMSIDM[reaction].In2].DeltaMass); //FIXME: last factor takes into account that sigma/m has to be rescaled for different particle masses
   }
  else
   {
     *retval = -1;
     return 0.0;
   }
}


/* note: we do not divide masses by HubbleParam; those cancel out during the calculation */
void sidm_get_delta_energy(unsigned char reaction, double *delta_E_1, double *delta_E_2)
{
  /* convert mass split to energy */
  /* note negative sign: endothermic = upscattering = negative sign   exothermic = downscattering = positive sign */
  *delta_E_1 = -(STSIDM[SMSIDM[reaction].Out1].DeltaMass - STSIDM[SMSIDM[reaction].In1].DeltaMass)  *  All.SIDM_clight * All.SIDM_clight * All.SIDM_GroundStateMass;
  *delta_E_2 = -(STSIDM[SMSIDM[reaction].Out2].DeltaMass - STSIDM[SMSIDM[reaction].In2].DeltaMass)  *  All.SIDM_clight * All.SIDM_clight * All.SIDM_GroundStateMass;
}

/* note: we do not divide masses by HubbleParam; those cancel out during the calculation */
void sidm_get_delta_mass(unsigned char reaction, double in_mass1, double in_mass2, double *out_mass1, double *out_mass2)
{
#ifndef SIDM_NO_MASSCHANGE
  *out_mass1 = (1.0 + STSIDM[SMSIDM[reaction].Out1].DeltaMass) * All.SIDM_GroundStateMass;
  *out_mass2 = (1.0 + STSIDM[SMSIDM[reaction].Out2].DeltaMass) * All.SIDM_GroundStateMass;
#else
  *out_mass1 = in_mass1; 
  *out_mass2 = in_mass2;
#endif
  if ((*out_mass1 == 0.0) || (*out_mass2 == 0.0))
    terminate("SIDM: wrong masses *out_mass1=%g *out_mass2=%g  reaction=%d   All.SIDM_GroundStateMass=%g\n", *out_mass1, *out_mass2, reaction, All.SIDM_GroundStateMass);
}


#endif // SIDM or not
