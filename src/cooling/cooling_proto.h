/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/cooling/cooling_proto.h
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

#ifndef INLINE_FUNC
#ifdef INLINE
#define INLINE_FUNC inline
#else
#define INLINE_FUNC
#endif
#endif

void SetOutputGasState(int i, double *ne_guess, double *nH0, double *coolrate);

double convert_u_to_temp(double u, double rho, double *ne_guess);
double CoolingRate(double logT, double rho, double *nelec);
double CoolingRateFromU(double u, double rho, double *ne_guess);
double DoCooling(double u_old, double rho, double dt, double *ne_guess);
double GetCoolingTime(double u_old, double rho, double *ne_guess);
double DoInstabilityCooling(double m_old, double u, double rho, double dt, double fac, double *ne_guess);
double DoSimpleCoolingHeating(int i, double dtcool);

#ifdef RADCOOL
void init_radcool_units(void);
void set_radcool_units_for_current_time(void);
#endif
void find_abundances_and_rates(double logT, double rho, double *ne_guess);
void InitCool(void);
void InitCoolMemory();
void IonizeParamsUVB(void);
void IonizeParams(void);
void IonizeParamsFunction(void);
void IonizeParamsTable(void);
double INLINE_FUNC LogTemp(double u, double ne);
void MakeCoolingTable(void);
void ReadIonizeParams(char *fname, int which);
void SetZeroIonization(void);
void TestCool(void);

#ifdef GFM_AGN_RADIATION
void IonizeParamsAGN(void);
#endif
#ifdef RADCOOL
void IonizeParamsRADCOOL(void);
#endif
#if !defined(RADCOOL) && ( defined(GFM_AGN_RADIATION) || defined(GFM_UVB_CORRECTIONS) || defined(UVB_SELF_SHIELDING))
void update_radiation_state(MyFloat rho, MyFloat xh, MyFloat bolIntensity);
void reset_radiation_state(void);
#endif
#ifdef RADCOOL
void update_radiation_state(MyFloat rho, MyFloat xh, MyFloat Phios, MyFloat Phins
#ifdef RADCOOL_HOTHALO
                            , MyFloat PhiT6, MyFloat PhiT7, MyFloat PhiT8
#endif
  );
void reset_radiation_state(void);
#endif
#ifdef GFM_COOLING_METAL
void update_gas_state(MyFloat rho, MyFloat xh, MyFloat log_metallicity);
#endif

#ifdef UVB_SELF_SHIELDING
double calc_self_shielding_factor(double nH);
void ReadSelfShieldingParams(char *fname);
#endif

#ifdef GFM_LAMBDA
void cooling_set_partindex(int i);
#endif
