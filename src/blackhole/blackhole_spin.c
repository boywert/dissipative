/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/blackhole/blackhole_spin.c
 * \date        09/2016
 * \author		Sebastian Bustamante
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - 16.09.2016 Creation of new module
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "../allvars.h"
#include "../proto.h"

/*! \file blackhole_spin.c
 *  \brief evolves BHs spin
 */


#if defined(BLACK_HOLES) && defined(BH_SPIN_EVOLUTION)

static int blackhole_do_prolonged_spin_evolution(int i, MyFloat CurrentTime);
static int blackhole_do_chaotic_spin_evolution(int i, MyFloat CurrentTime);
static MyFloat blackhole_do_spin_evolution(MyFloat a0, MyFloat Mass0, MyFloat Massf);

/* calculate spin properties once a new accretion episode takes place */
void blackhole_new_accretion_episode(int i, MyFloat CurrentTime, MyFloat BH_Mass)
{
  int k;
  MyFloat LambEdd, Xi;
  MyFloat a = BPP(i).BH_SpinParameter;
  MyFloat Z1, Z2, Rlso;
  MyFloat Mass = BH_Mass*All.UnitMass_in_g/1.989e41;
  MyFloat Rwarp, Rself;
  MyFloat v2v1, eps01, eps;
  MyFloat Alpha = All.ShakuraSunyaevParameter;
  MyFloat BH_TotalAngMomentum, BH_DirAngMomentum[3], NormalAccDisk[3];
  MyFloat Jbh, Jd;
  MyFloat DotProd;
  MyFloat BH_Theta;

  //Calculating viscosity ratio
  v2v1 = 2*( 1 + 7*Alpha)/( Alpha*Alpha*(4 + Alpha*Alpha) );

  //Updating radiative efficiency
  Z1 = 1 + pow(1 - a*a,1/3.)*(pow(1+a,1/3.0) + pow(1-a,1/3.0));
  Z2 = sqrt(3*a*a + Z1*Z1);
  Rlso = 3 + Z2 - fabs(a)/a*sqrt((3 - Z1)*(3 + Z1 + 2*Z2));
  BPP(i).BlackHoleRadiativeEfficiency = 1 - sqrt( 1 - 2/(3*Rlso) );
  eps = BPP(i).BlackHoleRadiativeEfficiency;
  LambEdd = eps*BPP(i).BH_Mdot/blackhole_luminosity_eddington(BH_Mass);
  Xi = BPP(i).BH_Mdot/blackhole_mdot_eddington(BH_Mass);
  eps01 = BPP(i).BlackHoleRadiativeEfficiency/0.1;

  //BH angular momentum
  Jbh = fabs(a)*Mass*Mass;

  //Disk angular momentum (Dubois, et al. 2014)
  //Rwarp = 4e2*pow(fabs(a),5/8.)*pow(Mass,1/8.)*pow(eps01/Xi,1/4.)*pow(v2v1/85.,-5/8.)*pow(Alpha/0.1,-1/2.);
  //Rself = 5e2*pow(Mass,-52/45.)*pow(eps01/Xi,22/45.)*pow(Alpha/0.1,28/45.);

  //Disk angular momentum (Fanidakis, et al. 2010)
  Rwarp = 3.6e3*pow(fabs(a),5/8.)*pow(Mass,1/8.)*pow(LambEdd,-1/4.)*pow(v2v1/85.,-5/8.)*pow(Alpha,-1/2.);
  Rself = 1.5e3*pow(eps,8/27.)*pow(Mass,-26/27.)*pow(LambEdd,-8/27.)*pow(Alpha,14/27.);

  //Switching spin modes
#if (BH_SPIN_MODEL == 3)
  if( Rself<Rwarp )
    BPP(i).BH_SpinModel = 1;
  else
	BPP(i).BH_SpinModel = 0;
#endif

  printf("HERE!!    %e  %e\t %e   %e   %e   %e   %e   %e \n", Rwarp, Rself, fabs(a), Mass, eps01, LambEdd, v2v1, Alpha );
  //Evolving spin
  if( BPP(i).BH_FlagOngAccEpis == 1 )
    {
	  if( BPP(i).BH_SpinModel == 0 )
	   {
		 MyFloat AngMomGasCells_Norm = 0;
		 for(k=0; k<3; k++)
		  AngMomGasCells_Norm += BPP(i).BH_AngMomGasCells[k]*BPP(i).BH_AngMomGasCells[k];
		 AngMomGasCells_Norm = sqrt(AngMomGasCells_Norm);
		 for(k=0; k<3; k++)
		  NormalAccDisk[k] = BPP(i).BH_AngMomGasCells[k]/AngMomGasCells_Norm;
	   }
	  else
	   {
		 double Phi, CosTheta, SinTheta;
		 Phi = 2*M_PI*get_random_number();
		 CosTheta = 1-2*get_random_number();
		 SinTheta = sqrt(1-CosTheta*CosTheta);

		 NormalAccDisk[0] = SinTheta*cos(Phi);
		 NormalAccDisk[1] = SinTheta*sin(Phi);
		 NormalAccDisk[2] = CosTheta;
	   }

	  if(BPP(i).BH_SpinModel == 1)
		  Jd = Mass*(Mass - BPP(i).BH_Mass_Previous*All.UnitMass_in_g/1.989e41)*pow(Rself,1/2.);
	  else
		  Jd = Mass*(Mass - BPP(i).BH_Mass_Previous*All.UnitMass_in_g/1.989e41)*pow(Rwarp,1/2.);

	  DotProd = 0;
	  for(k=0; k<3; k++)
	   DotProd += BPP(i).BH_SpinOrientation[k]*NormalAccDisk[k];
	  BH_Theta = DotProd;

	  if(BH_Theta<-0.5*Jd/Jbh)
		BPP(i).BH_SpinParameter = -blackhole_do_spin_evolution(-a, BPP(i).BH_Mass_Previous, BH_Mass);
	  else
		BPP(i).BH_SpinParameter = blackhole_do_spin_evolution(a, BPP(i).BH_Mass_Previous, BH_Mass);

	  BH_TotalAngMomentum = 0;
	  for(k=0; k<3; k++)
	   {
		 BH_DirAngMomentum[k] = Jbh*BPP(i).BH_SpinOrientation[k] + Jd*NormalAccDisk[k];
		 BH_TotalAngMomentum += BH_DirAngMomentum[k]*BH_DirAngMomentum[k];
	   }
	  BH_TotalAngMomentum = sqrt(BH_TotalAngMomentum);
	  for(k=0; k<3; k++)
	   BPP(i).BH_SpinOrientation[k] = BH_DirAngMomentum[k]/BH_TotalAngMomentum;
    }

  BPP(i).BH_TimeAccretion_Previous = CurrentTime;
  BPP(i).BH_Mass_Previous = BH_Mass;

  //Updating accretion time in prolonged mode and accreted mass in chaotic mode (Dubois, et al. 2014)
  //BPP(i).BH_DTimeAccretion_Current = 5.3e5*pow(fabs(a),7/8.)*pow(Mass,11/8.)*pow(eps01/Xi,3/4.)*pow(v2v1,-7/8.)*pow(Alpha,-3/2.)*SEC_PER_YEAR/SEC_PER_GIGAYEAR; //in yrs
  //BPP(i).BH_DMass_Current = 6.0e5*pow(Alpha/0.1,-1/45.)*pow(Mass,34/45.)*pow(Xi/eps01,4/45.)*SOLAR_MASS/All.UnitMass_in_g;

  //Updating accretion time in prolonged mode and accreted mass in chaotic mode (Fanidakis, et al. 2010)
  BPP(i).BH_DTimeAccretion_Current = 3.0e6*pow(fabs(a),7/8.)*pow(Mass,11/8.)*pow(LambEdd,-3/4.)*pow(v2v1,-7/8.)*pow(Alpha,-3/2.)*SEC_PER_YEAR/SEC_PER_GIGAYEAR; //in yrs
  BPP(i).BH_DMass_Current = 2.13e5*pow(eps,-5/27.)*pow(Mass,23/27.)*pow(LambEdd,5/27.)*pow(Alpha,-2/17.)*SOLAR_MASS/All.UnitMass_in_g;

}


/* calculate spin parameter */
MyFloat blackhole_do_spin_evolution(MyFloat a0, MyFloat Mass0, MyFloat Massf)
{
  MyFloat Z1, Z2, Rlso;
  MyFloat a;

  //calculating Radius of last stable orbit
  Z1 = 1 + pow(1 - a0*a0,1/3.)*(pow(1+a0,1/3.0) + pow(1-a0,1/3.0));
  Z2 = sqrt(3*a0*a0 + Z1*Z1);
  Rlso = 3 + Z2 - fabs(a0)/a0*sqrt((3 - Z1)*(3 + Z1 + 2*Z2));

  //new spin parameter
  if(Massf/Mass0 <= sqrt(Rlso) && a0 < 0.998)
	a = (1/3.)*sqrt(Rlso)*Mass0/Massf*(4 - sqrt(3*Rlso*pow(Mass0/Massf,2)-2));
  else
	a = 0.998;

  return a;
}


/* calculate prolonged spin model */
int blackhole_do_prolonged_spin_evolution(int i, MyFloat CurrentTime)
{
  MyFloat BH_Mass;
  MyFloat DTimeCurrentEpisode = CurrentTime - BPP(i).BH_TimeAccretion_Previous;

  if(DTimeCurrentEpisode >= BPP(i).BH_DTimeAccretion_Current)
   {
     BH_Mass = BPP(i).BH_Mass_Previous + (BPP(i).BH_Mass-BPP(i).BH_Mass_Previous)*BPP(i).BH_DTimeAccretion_Current/DTimeCurrentEpisode;
     blackhole_new_accretion_episode(i, BPP(i).BH_TimeAccretion_Previous + BPP(i).BH_DTimeAccretion_Current, BH_Mass);
     return 1;
   }
  else
    return 0;
}


/* calculate chaotic spin model */
int blackhole_do_chaotic_spin_evolution(int i, MyFloat CurrentTime)
{
  MyFloat BH_Mass;
  MyFloat DMassCurrentEpisode = BPP(i).BH_Mass - BPP(i).BH_Mass_Previous;

  if(DMassCurrentEpisode >= BPP(i).BH_DMass_Current)
   {
	  BH_Mass = BPP(i).BH_Mass_Previous + BPP(i).BH_DMass_Current;
	  blackhole_new_accretion_episode(i, CurrentTime, BH_Mass);
	  return 1;
   }
  else
    return 0;
}


/* spin evolution for all active BHs */
void blackhole_spin_evolution(void)
{
  MyFloat MdotNorm;
  MyFloat CurrentTime;
  int TimeBreak = 1;

  for(int idx = 0; idx < TimeBinsBHAccretion.NActiveParticles; idx++)
   {
      int i = TimeBinsBHAccretion.ActiveParticleList[idx];
      if(i < 0)
        continue;
      if(BPP(i).BH_Mass < All.SeedBlackHoleMass)
      	continue;

      if(All.ComovingIntegrationOn)
        CurrentTime = All.TimeBegin * exp(P[i].Ti_Current * All.Timebase_interval);
      else
        CurrentTime = All.TimeBegin + P[i].Ti_Current * All.Timebase_interval;
      CurrentTime = get_time_difference_in_Gyr(All.TimeBegin, CurrentTime);

      MdotNorm = BPP(i).BH_Mdot/blackhole_mdot_eddington(BPP(i).BH_Mass);

      if(BPP(i).BH_FlagOngAccEpis == 0 && MdotNorm>=1e-5)
       {
         blackhole_new_accretion_episode(i, CurrentTime, BPP(i).BH_Mass);
         BPP(i).BH_FlagOngAccEpis = 1;
       }
      if(MdotNorm<1e-5)
       {
  	     //ADAF MODEL OR ALTERNATIVE MODEL FOR LOW ACCRETION REGIMES
         BPP(i).BH_FlagOngAccEpis = 0;
         continue;
       }

      //Evolving spin for current timestep
      if(BPP(i).BH_FlagOngAccEpis == 1)
        while(TimeBreak == 1)
         {
           fprintf(FdBlackHolesSpin, "BH=%llu %g %g %g %g %g %g %g %g %g %g %d\n",
        		   (long long) P[i].ID, All.Time,
        		   BPP(i).BH_SpinParameter, BPP(i).BH_Mass, MdotNorm,
				   BPP(i).BH_SpinOrientation[0], BPP(i).BH_SpinOrientation[1], BPP(i).BH_SpinOrientation[2],
				   BPP(i).BH_AngMomGasCells[0], BPP(i).BH_AngMomGasCells[1], BPP(i).BH_AngMomGasCells[2],
				   BPP(i).BH_SpinModel);
#if (BH_SPIN_MODEL == 0)
    	   BPP(i).BH_SpinModel = 0;
    	   TimeBreak = blackhole_do_prolonged_spin_evolution(i, CurrentTime);
#endif
#if (BH_SPIN_MODEL == 1)
           BPP(i).BH_SpinModel = 1;
           TimeBreak = blackhole_do_chaotic_spin_evolution(i, CurrentTime);
#endif
#if (BH_SPIN_MODEL == 2 || BH_SPIN_MODEL == 3)
           if(BPP(i).BH_SpinModel == 0)
            {
              TimeBreak = blackhole_do_prolonged_spin_evolution(i, CurrentTime);
            }
           else
             TimeBreak = blackhole_do_chaotic_spin_evolution(i, CurrentTime);
#endif
         }
   }
}

#endif
