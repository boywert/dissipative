/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/barpot/barpot.cc
 * \date        09/2016
 * \author      Mattia Sormani, Robin Tress 
 * \brief       This code calculates the contribution of different galactic components to the potential.
 *
 * \details     
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */

#include "common.h"
#include "potential.h"
#include "../allvars.h"

ExpBarPotential* Barpot;
ExpBarPotential* Barpot_axi;
ExpDiskPotential* Diskpot1;
ExpDiskPotential* Diskpot2;
DPLPotential* Halopot;
ModifiedMcMillanBulgePotential* Bulgepot;

#ifdef __cplusplus
extern "C" 
{
#endif
  void galpot_init(void)
  {
    double G = All.G;   

  /*Initializing bar potential*/
    double rminbar = 0., rmaxbar = 0.;
    int nrbar = 500, npolybar = 6, ngaussbar = 8;
    double rho0bar = 5.e6, a0bar = 7.5, qbar = 0.5;

    ExpBarDensity   rhobar(rho0bar, a0bar, qbar, rminbar, rmaxbar);
    Barpot = new ExpBarPotential(G, rhobar, nrbar, npolybar, ngaussbar);
    
    ExpBarDensity   rhobar_axi(rho0bar*qbar*qbar, a0bar, 1., rminbar, rmaxbar);
    Barpot_axi = new ExpBarPotential(G, rhobar_axi, nrbar, npolybar, ngaussbar);

  /*Initializing Bulge potential*/
    double rminbulge = 0., rmaxbulge = 0.;
    int nrbulge = 500, npolybulge = 6, ngaussbulge = 8;
    double rho0bulge = 8.0e5, a0bulge=10.0, acutbulge = 10.0;
    double alfbulge = 1.7, qbulge = 0.5;

    ModifiedMcMillanBulgeDensity   rhobulge(rho0bulge, alfbulge, a0bulge, acutbulge, qbulge, rminbulge, rmaxbulge);
    Bulgepot = new ModifiedMcMillanBulgePotential(G, rhobulge, nrbulge, npolybulge, ngaussbulge);

  /*Initializing exponential disk potential 
   *TODO: multipole expansion not ideal for disks. Consider switching to GalPot by Paul McMillan and Walter Dehnen */
    /*Thick stellar disk*/
    double rmind1 = 0., rmaxd1 = 0.;
    int nrd1 = 500, npolyd1 = 20, ngaussd1 = 8;
    double Sigma0d1 = 1.74e6, Rd1 = 30.2, zd1 = 9.;

    ExpDiskDensity   rhodisk1(Sigma0d1, zd1, Rd1, rmind1, rmaxd1);
    Diskpot1 = new ExpDiskPotential(G, rhodisk1, nrd1, npolyd1, ngaussd1);

    /*Thin stellar disk*/
    double rmind2 = 0., rmaxd2 = 0.;
    int nrd2 = 500, npolyd2 = 20, ngaussd2 = 8;
    double Sigma0d2 = 8.5e6, Rd2 = 25.0, zd2 = 3.0;

    ExpDiskDensity   rhodisk2(Sigma0d2, zd2, Rd2, rmind2, rmaxd2);
    Diskpot2 = new ExpDiskPotential(G, rhodisk2, nrd2, npolyd2, ngaussd2);

  /*Initializing Halo potential*/
    double rminh = 0., rmaxh = 0.;
    int nrh = 500, npolyh = 6, ngaussh = 8;
    double rho0h = 8.11e3, ah = 196., alphah = 1.0, betah = 3.0;

    DPLDensity   rhohalo(rho0h, ah, alphah, betah, rminh, rmaxh);
    Halopot = new DPLPotential(G, rhohalo, nrh, npolyh, ngaussh);

    return;
  }

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" 
{
#endif
  void galpot_dPhidx(double x, double y, double z, double time, double *dPhidx)
  {
    Vec3 X(x,y,z);

    dPhidx[0] = 0.; 
    dPhidx[1] = 0.; 
    dPhidx[2] = 0.;

  /*Bar potential*/
    double omega_b = -4.; //internal units
    double t_grow = 1.5; //internal units

    double costheta = cos(omega_b * time);
    double sintheta = sin(omega_b * time);
    Vec3 Xbar((x*costheta + y*sintheta), (-x*sintheta + y*costheta), z);

    Vec3 dBarPhidx = Barpot->dPhidx(Xbar);  
    Vec3 dBarPhidx_axi = Barpot_axi->dPhidx(Xbar);
    Vec3 bar_dPhidx;
    
    if (time > t_grow)
      {
        bar_dPhidx[0] = dBarPhidx[0];
        bar_dPhidx[1] = dBarPhidx[1];  
        bar_dPhidx[2] = dBarPhidx[2];
      }
    else
      {
        bar_dPhidx[0] = (time / t_grow) * dBarPhidx[0] + (1. - (time / t_grow)) * dBarPhidx_axi[0];
        bar_dPhidx[1] = (time / t_grow) * dBarPhidx[1] + (1. - (time / t_grow)) * dBarPhidx_axi[1];
        bar_dPhidx[2] = (time / t_grow) * dBarPhidx[2] + (1. - (time / t_grow)) * dBarPhidx_axi[2];
      }

    dPhidx[0] += bar_dPhidx[0] * costheta - bar_dPhidx[1] * sintheta;
    dPhidx[1] += bar_dPhidx[0] * sintheta + bar_dPhidx[1] * costheta;
    dPhidx[2] += bar_dPhidx[2];

  /*Bulge potential*/
    Vec3 dBulgePhidx = Bulgepot->dPhidx(X);
    
    dPhidx[0] += dBulgePhidx[0];
    dPhidx[1] += dBulgePhidx[1];
    dPhidx[2] += dBulgePhidx[2];

  /*Disk potential*/ 
    /*Thick disk*/
    Vec3 dDisk1Phidx = Diskpot1->dPhidx(X);

    dPhidx[0] += dDisk1Phidx[0];
    dPhidx[1] += dDisk1Phidx[1];
    dPhidx[2] += dDisk1Phidx[2];

    /*Thin disk*/
    Vec3 dDisk2Phidx = Diskpot2->dPhidx(X);

    dPhidx[0] += dDisk2Phidx[0];
    dPhidx[1] += dDisk2Phidx[1];
    dPhidx[2] += dDisk2Phidx[2];

  /*Halo potential*/
    Vec3 dDPLPhidx = Halopot->dPhidx(X);

    dPhidx[0] += dDPLPhidx[0];
    dPhidx[1] += dDPLPhidx[1];
    dPhidx[2] += dDPLPhidx[2];

    return;
  }
#ifdef __cplusplus
}
#endif

