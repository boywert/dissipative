#ifndef SGCHEM_PROTO_H
#define SGCHEM_PROTO_H
#include "f2c.h"

void CALC_DUST_TEMP(double* yn, double* chi_mean, double* temp, double* abHI, double* RH2, double* Rdust, int* nocalc);
void CALC_SHIELDING(double* yn, double* dl, double* temp, double* abh2, double* abco, double column_density_projection[NPIX],
		    double column_density_projection_H2[NPIX], double column_density_projection_co[NPIX], double* fshield_H2,
		    double* fshield_CO, double* AV_mean, double* chi_mean);
void CHEMINMO(void);
void COOLINMO(void);
void do_chemical_abund_checks(struct particle_data *P, struct sph_particle_data *SphP, int i, struct pv_update_data *pvd);
void EVOLVE_ABUNDANCES(double* time, double* dl, double* yn, double* divv, double* energy, double* redshift,
		       double non_eq_abundances[SGCHEM_NUM_SPECIES], double breakdown_of_thermal_rates[SGCHEM_NUM_THERMAL_RATES], 
		       double column_density_projection[NPIX], double column_density_projection_H2[NPIX], 
		       double column_density_projection_co[NPIX], double* dust_temp, int* id);
void evolve_chemistry(void);
void evolve_chemistry_for_single_cell(double dt, int index);
void init_chemistry(void);
void INIT_CHEMISTRY_PARAMETERS(double* deuterium,
#ifndef SGCHEM_VARIABLE_Z
                               double* carbon, double* oxygen, double* misc_metal, double* Z_atom,
#endif
                               double* tdust,
#ifndef SGCHEM_VARIABLE_ISRF
                               double* UVfield,
#endif
                               long* lwbgtype, double* lwbgstart,
#ifndef SGCHEM_VARIABLE_Z 
                               double* dust_to_gas,
#endif
#ifndef SGCHEM_VARIABLE_CRION
                               double* cr_ion, 
#endif
                               double* redshift, double* AV_ext, double* h2_form_ex,
                               double* h2_form_kin, long* iphoto, long* isrf, long* iflag_atom,
                               long* iflag_H2_opacity);
void lwbg(void);
void INIT_TOLERANCES(void);
void LOAD_H2_TABLE(void);
void SET_INDEX_ID_FOR_CHEMISTRY(int* index, int* id);
void SET_LOCAL_DUST_ABUNDANCE(double* dust_to_gas);
void SET_LOCAL_ELEMENT_ABUNDANCES(double* carb_abund, double* oxy_abund, double* m_abund, double* Z_atom);
void SET_LOCAL_ISRF(double* local_ISRF);
void SET_LOCAL_CR_ION_RATE(double* CR_energy_density, double* local_CR_ion_rate);
int tgset_jeans_ref(int i);
#ifdef SGCHEM_VARIABLE_ISRF
void find_and_set_local_ISRF(int index);
#endif
#ifdef SGCHEM_VARIABLE_CRION
void find_and_set_local_CRION(int index);
#endif

#ifdef MCMA
void CMA_INIT(void);
void CMA_CORRECT(double species_flux[NSPEC_CMA], double element_flux[NELEM_CMA], long* sgflux_count);
#endif

#endif
