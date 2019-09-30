/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/GFM/stellar_evolution_proto.h
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

#ifndef STELLAR_EVOLUTION_H
#define STELLAR_EVOLUTION_H

/* helper functions */
int element_index(char *element_name);
int get_element_index(char **table, int size, char *element_name);
void get_z_indicies(MyFloat log_metallicity, MyFloat * metal_values, int N_Z, int *iz_low, int *iz_high, MyFloat * dz);
void convert_metal_abundances_to_masses(void);
MyFloat ejecta_mass_for_bin(MyFloat ** ejecta_spline, MyFloat *** yield_spline, int elem_index, MyFloat * initial_metals, int iz_low, int iz_high, int imass, MyFloat dz);

/* IMF and lifetime */
void init_imf(void);
double define_imf(double imf_param);
double calc_FactorSN(void);
double calc_WindEgySpecSN(void);
void init_SNIa_rates(void);
void get_imf_bins(MyDouble log_min_mass, MyDouble log_max_mass, int *ilow, int *ihigh);
MyDouble integrate_imf(MyDouble log_min_mass, MyDouble log_max_mass, int mode, double *imf_integrand_mass);
MyFloat get_dying_mass_in_Msun(MyFloat time_in_Gyr, MyFloat metallicity);
MyFloat get_lifetime_in_Gyr(MyDouble mass, MyFloat metallicity);
void get_initial_mass_fractions(MyDouble * mass_fractions);

/* stellar population */
void evolve_SNIa(MyDouble log_min_mass, MyDouble log_max_mass, MyFloat age_of_star_in_Gyr, MyFloat dt_in_Gyr, MyFloat metallicity, stellar_evolution_data * sed);
void evolve_NSNS(MyDouble log_min_mass, MyDouble log_max_mass, MyFloat age_of_star_in_Gyr, MyFloat dt_in_Gyr, MyFloat metallicity, MyFloat initial_mass, stellar_evolution_data * sed);
void evolve_SNII(MyDouble log_min_mass, MyDouble log_max_mass, MyFloat log_metallicity, MyFloat * initial_metals, stellar_evolution_data * sed);
void evolve_AGB(MyDouble log_min_mass, MyDouble log_max_mass, MyFloat log_metallicity, MyFloat * initial_metals, stellar_evolution_data * sed);

void do_stellar_evolution(MyFloat age_of_star_in_Gyr, MyFloat dtime_in_Gyr, int iPart, stellar_evolution_data * sed);
void test_stellar_evolution(void);
void evolve_active_stars(void);

/* log files */
void output_stellar_evolution_statistics(void);
#if defined(FM_STAR_FEEDBACK) && defined(OUTPUT_STELLAR_FEEDBACK)
void update_cumulative_feedback_energy_of_converted_cells(int index);
void output_stellar_feedback_statistics(void);
#endif

/* yields I/O and init */
int read_yield_tables(void);
void init_yields(void);

/* cell finding and enrichment */
void find_cells_dump(int npart);
void find_cells_to_enrich(void);

void do_chemical_enrichment(void);

void find_star_cells(void);
int find_star_cells_evaluate(int target, int mode, int *nexport, int *nsend_local);

#ifdef GFM_DUST
double get_cell_dtime_Gyr(int i);
void dust_growth_and_destruction(void);
#endif

#if defined(FM_STAR_FEEDBACK) && defined(DELAYED_COOLING)
void find_feedback_cells(void);
int find_feedback_cells_evaluate(int target, int mode, int thread_id);
#endif

/* stellar feedback */
void do_stellar_feedback(void);
int stellar_feedback_evaluate(int target, int mode, int thread_id);

/* allocation and setup */
void start_enrichment(void);
void end_enrichment(void);

#ifdef GFM_PREENRICH
void gfm_preenrich_gas(void);
void gfm_read_preenrich_table(char *fname);
#endif

#ifdef GFM_WINDS_STRIPPING
void find_cells_to_strip(void);
void start_stripping(void);
void end_stripping(void);
void do_chemical_stripping(void);
void strip_active_winds(void);
#endif

/* interpolation */
static inline MyFloat interpol_1d(MyFloat * table, int i, float dx);
static inline MyFloat interpol_2d(MyFloat ** table, int i, int j, double dx, double dy);
static inline MyFloat interpol_3d(MyFloat *** table, int i, int j, int k, double dx, double dy, double dz);
static inline MyFloat interpol_4d(MyFloat **** table, int i, int j, int k, int l, double dx, double dy, double dz, double dw);
#ifndef RADCOOL
static inline MyFloat interpol_5d(MyFloat ***** table, int i, int j, int k, int l, int m, double dx, double dy, double dz, double dw, double dv);
#endif

#ifdef RADCOOL
static inline MyFloat interpol_5d(float *****table, int i, int j, int k, int l, int m, double dx, double dy, double dz, double dw, double dv);
#ifdef RADCOOL_HOTHALO
static inline MyFloat interpol_8d(float ********table, int ir, int io, int in, int i6, int i7, int i8, int id, int it, double dir, double dio, double din, double di6, double di7, double di8,
                                  double did, double dit);
#endif
#endif

static inline MyFloat interpol_1d(MyFloat * table, int i, float dx)
{
  return (1 - dx) * table[i] + dx * table[i + 1];
}

static inline MyFloat interpol_2d(MyFloat ** table, int i, int j, double dx, double dy)
{
  return (1 - dx) * (1 - dy) * table[i][j] + (1 - dx) * dy * table[i][j + 1] + dx * (1 - dy) * table[i + 1][j] + dx * dy * table[i + 1][j + 1];
}

static inline MyFloat interpol_3d(MyFloat *** table, int i, int j, int k, double dx, double dy, double dz)
{
  int il = i, jl = j, kl = k;
  int ir = i + 1, jr = j + 1, kr = k + 1;

  double dxl = 1 - dx, dyl = 1 - dy, dzl = 1 - dz;
  double dxr = dx, dyr = dy, dzr = dz;

  if(dxr == 0)
    ir = i;

  if(dyr == 0)
    jr = j;

  if(dzr == 0)
    kr = k;

  return dxl * dyl * dzl * table[il][jl][kl] +
    dxl * dyl * dzr * table[il][jl][kr] +
    dxl * dyr * dzl * table[il][jr][kl] +
    dxl * dyr * dzr * table[il][jr][kr] + dxr * dyl * dzl * table[ir][jl][kl] + dxr * dyl * dzr * table[ir][jl][kr] + dxr * dyr * dzl * table[ir][jr][kl] + dxr * dyr * dzr * table[ir][jr][kr];
}

static inline MyFloat interpol_4d(MyFloat **** table, int i, int j, int k, int l, double dx, double dy, double dz, double dw)
{
  int il = i, jl = j, kl = k, ll = l;
  int ir = i + 1, jr = j + 1, kr = k + 1, lr = l + 1;

  double dxl = 1 - dx, dyl = 1 - dy, dzl = 1 - dz, dwl = 1 - dw;
  double dxr = dx, dyr = dy, dzr = dz, dwr = dw;
  if(dxr == 0)
    ir = i;

  if(dyr == 0)
    jr = j;

  if(dzr == 0)
    kr = k;

  if(dwr == 0)
    lr = l;

  return
    dxl * dyl * dzl * dwl * table[il][jl][kl][ll] +
    dxl * dyl * dzl * dwr * table[il][jl][kl][lr] +
    dxl * dyl * dzr * dwl * table[il][jl][kr][ll] +
    dxl * dyl * dzr * dwr * table[il][jl][kr][lr] +
    dxl * dyr * dzl * dwl * table[il][jr][kl][ll] +
    dxl * dyr * dzl * dwr * table[il][jr][kl][lr] +
    dxl * dyr * dzr * dwl * table[il][jr][kr][ll] +
    dxl * dyr * dzr * dwr * table[il][jr][kr][lr] +
    dxr * dyl * dzl * dwl * table[ir][jl][kl][ll] +
    dxr * dyl * dzl * dwr * table[ir][jl][kl][lr] +
    dxr * dyl * dzr * dwl * table[ir][jl][kr][ll] +
    dxr * dyl * dzr * dwr * table[ir][jl][kr][lr] +
    dxr * dyr * dzl * dwl * table[ir][jr][kl][ll] + dxr * dyr * dzl * dwr * table[ir][jr][kl][lr] + dxr * dyr * dzr * dwl * table[ir][jr][kr][ll] + dxr * dyr * dzr * dwr * table[ir][jr][kr][lr];
}

#ifndef RADCOOL
inline MyFloat interpol_5d(MyFloat ***** table, int i, int j, int k, int l, int m, double dx, double dy, double dz, double dw, double dv)
#else
inline MyFloat interpol_5d(float *****table, int i, int j, int k, int l, int m, double dx, double dy, double dz, double dw, double dv)
#endif
{
  int il = i, jl = j, kl = k, ll = l, ml = m;
  int ir = i + 1, jr = j + 1, kr = k + 1, lr = l + 1, mr = m + 1;


  double dxl = 1 - dx, dyl = 1 - dy, dzl = 1 - dz, dwl = 1 - dw, dvl = 1 - dv;
  double dxr = dx, dyr = dy, dzr = dz, dwr = dw, dvr = dv;

  if(dxr == 0)
    ir = i;

  if(dyr == 0)
    jr = j;

  if(dzr == 0)
    kr = k;

  if(dwr == 0)
    lr = l;

  if(dvr == 0)
    mr = m;

  return
    dxl * dyl * dzl * dwl * dvl * table[il][jl][kl][ll][ml] +
    dxl * dyl * dzl * dwl * dvr * table[il][jl][kl][ll][mr] +
    dxl * dyl * dzl * dwr * dvl * table[il][jl][kl][lr][ml] +
    dxl * dyl * dzl * dwr * dvr * table[il][jl][kl][lr][mr] +
    dxl * dyl * dzr * dwl * dvl * table[il][jl][kr][ll][ml] +
    dxl * dyl * dzr * dwl * dvr * table[il][jl][kr][ll][mr] +
    dxl * dyl * dzr * dwr * dvl * table[il][jl][kr][lr][ml] +
    dxl * dyl * dzr * dwr * dvr * table[il][jl][kr][lr][mr] +
    dxl * dyr * dzl * dwl * dvl * table[il][jr][kl][ll][ml] +
    dxl * dyr * dzl * dwl * dvr * table[il][jr][kl][ll][mr] +
    dxl * dyr * dzl * dwr * dvl * table[il][jr][kl][lr][ml] +
    dxl * dyr * dzl * dwr * dvr * table[il][jr][kl][lr][mr] +
    dxl * dyr * dzr * dwl * dvl * table[il][jr][kr][ll][ml] +
    dxl * dyr * dzr * dwl * dvr * table[il][jr][kr][ll][mr] +
    dxl * dyr * dzr * dwr * dvl * table[il][jr][kr][lr][ml] +
    dxl * dyr * dzr * dwr * dvr * table[il][jr][kr][lr][mr] +
    dxr * dyl * dzl * dwl * dvl * table[ir][jl][kl][ll][ml] +
    dxr * dyl * dzl * dwl * dvr * table[ir][jl][kl][ll][mr] +
    dxr * dyl * dzl * dwr * dvl * table[ir][jl][kl][lr][ml] +
    dxr * dyl * dzl * dwr * dvr * table[ir][jl][kl][lr][mr] +
    dxr * dyl * dzr * dwl * dvl * table[ir][jl][kr][ll][ml] +
    dxr * dyl * dzr * dwl * dvr * table[ir][jl][kr][ll][mr] +
    dxr * dyl * dzr * dwr * dvl * table[ir][jl][kr][lr][ml] +
    dxr * dyl * dzr * dwr * dvr * table[ir][jl][kr][lr][mr] +
    dxr * dyr * dzl * dwl * dvl * table[ir][jr][kl][ll][ml] +
    dxr * dyr * dzl * dwl * dvr * table[ir][jr][kl][ll][mr] +
    dxr * dyr * dzl * dwr * dvl * table[ir][jr][kl][lr][ml] +
    dxr * dyr * dzl * dwr * dvr * table[ir][jr][kl][lr][mr] +
    dxr * dyr * dzr * dwl * dvl * table[ir][jr][kr][ll][ml] +
    dxr * dyr * dzr * dwl * dvr * table[ir][jr][kr][ll][mr] + dxr * dyr * dzr * dwr * dvl * table[ir][jr][kr][lr][ml] + dxr * dyr * dzr * dwr * dvr * table[ir][jr][kr][lr][mr];
}


#ifdef RADCOOL_HOTHALO
inline MyFloat interpol_8d(float ********table, int ir, int io, int in, int i6, int i7, int i8, int id, int it, double dir, double dio, double din, double di6, double di7, double di8, double did,
                           double dit)
{
  int idx[8][2];
  double dval[8][2];

  int cn0, cn1, cn2, cn3, cn4, cn5, cn6, cn7;


  int i1;


  idx[0][0] = ir;
  idx[1][0] = io;
  idx[2][0] = in;
  idx[3][0] = i6;
  idx[4][0] = i7;
  idx[5][0] = i8;
  idx[6][0] = id;
  idx[7][0] = it;

  dval[0][1] = dir;
  dval[1][1] = dio;
  dval[2][1] = din;
  dval[3][1] = di6;
  dval[4][1] = di7;
  dval[5][1] = di8;
  dval[6][1] = did;
  dval[7][1] = dit;

  for(i1 = 0; i1 < 8; i1++)
    {
      dval[i1][0] = 1 - dval[i1][1];
      if(dval[i1][1] == 0)
        idx[i1][1] = idx[i1][0];
      else
        idx[i1][1] = idx[i1][0] + 1;
    }

  double sum = 0.0;
  for(cn0 = 0; cn0 < 2; cn0++)
    for(cn1 = 0; cn1 < 2; cn1++)
      for(cn2 = 0; cn2 < 2; cn2++)
        for(cn3 = 0; cn3 < 2; cn3++)
          for(cn4 = 0; cn4 < 2; cn4++)
            for(cn5 = 0; cn5 < 2; cn5++)
              for(cn6 = 0; cn6 < 2; cn6++)
                for(cn7 = 0; cn7 < 2; cn7++)
                  sum +=
                    dval[0][cn0] * dval[1][cn1] * dval[2][cn2] * dval[3][cn3] * dval[4][cn4] * dval[5][cn5] * dval[6][cn6] * dval[7][cn7] *
                    table[idx[0][cn0]][idx[1][cn1]][idx[2][cn2]][idx[3][cn3]][idx[4][cn4]][idx[5][cn5]][idx[6][cn6]][idx[7][cn7]];

  return sum;
}
#endif

#endif
