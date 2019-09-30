#include "allvars.h"
#include "proto.h"

#ifdef MHD

static void do_mhd_source_terms(void);
#ifdef MHD_POWELL_SPLIT
static void do_mhd_powell_source_terms(void);
#endif

void do_mhd_source_terms_first_half(void)
{
#ifdef MHD_POWELL_SPLIT
  do_mhd_powell_source_terms();
#endif
  do_mhd_source_terms();
  update_primitive_variables();
}

void do_mhd_source_terms_second_half(void)
{
  do_mhd_source_terms();
#ifdef MHD_POWELL_SPLIT
  do_mhd_powell_source_terms();
#endif
  update_primitive_variables();
}

void do_mhd_source_terms(void)
{
  TIMER_START(CPU_MHD);

  if(All.ComovingIntegrationOn)
    {
      double atime = All.Time;
      double hubble_a = hubble_function(atime);

      int idx, i;
      for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
        {
          i = TimeBinsHydro.ActiveParticleList[idx];
          if(i < 0)
            continue;

          double dt_cell = 0.5 * (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval / hubble_a; /* half the timestep of the cell */
          SphP[i].Energy += dt_cell * 0.5 * (SphP[i].B[0] * SphP[i].B[0] + SphP[i].B[1] * SphP[i].B[1] + SphP[i].B[2] * SphP[i].B[2]) * SphP[i].Volume * atime * hubble_a;
        }
    }
  
  TIMER_STOP(CPU_MHD);
}

#ifdef MHD_POWELL_SPLIT
void do_mhd_powell_source_terms(void)
{
  TIMER_START(CPU_MHD);

  if(All.ComovingIntegrationOn)
    terminate( "do_mhd_powell_source_terms still lacks the cosmological factors" );

  set_cosmo_factors_for_current_time();
  calculate_gradients(); // this updates divB
  
  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;
      
      double dt_cell = 0.5 * (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval / All.cf_hubble_a; /* half the timestep of the cell */
      dt_cell *= SphP[i].DivB;
      
      for(int k = 0; k < 3; k++)
        {
          double dMomentum = - SphP[i].B[k] * dt_cell;
          SphP[i].Momentum[k] += dMomentum;
          All.Powell_Momentum[k] += dMomentum;
          
          double dEnergy = - P[i].Vel[k] * SphP[i].B[k] * dt_cell;
          SphP[i].Energy += dEnergy;
          All.Powell_Energy += dEnergy;
          
          double dB = - P[i].Vel[k] * dt_cell;
          SphP[i].BConserved[k] += dB;
        }
    }
  
  TIMER_STOP(CPU_MHD);
}
#endif

#endif
