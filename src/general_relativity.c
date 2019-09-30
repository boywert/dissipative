/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/general_relativity.c
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

/* TODO or TOCHECK

 -> limits that specific energy remains positive
 -> higher order extrapolation of cell motion for better angular momentum conservation
 -> EoS call in HLLC
 -> is it problematic to unsplit update when mesh's moving?
 -> if we unsplit the time evolution, is this compatible with moving mesh???
 -> ACHTUNG bei der Zeitintegration haben wir womoeglich ueberall einen Fehler 2 drin, so dass das nicht auffaefllt aber sich einfach die Zeit anders entwickelt. Das koennte dann aber beim moving mesh zu problem fuehren, je nachdem wie es da entwickelt wird.

 -> microphys EoS
 -> couple metric solver
 -> GW extraction

 */


#include "allvars.h"
#include "proto.h"
#include <math.h>

#ifdef GENERAL_RELATIVITY

// todo if adaiabatic, we need to think how the pressure is computed on the faces, consistent with EoS?

//static double BH_x = 27.392*0.5-0.105353846154; //15.;   // this will be used as center if we use a fixed numerical metric that is read in
//static double BH_y = 27.392*0.5-0.105353846154; //15.;   // -> so it will be overwritten
//static double BH_z = 27.392*0.5-0.105353846154; //15.;
static double BH_x = 27.392*0.5;
static double BH_y = 27.392*0.5;
static double BH_z = 27.392*0.5;

static double BH_Mass = 1.;
static double BH_Spin = 0.0;

static int GR_counter=0;

static double kpoly_gr = 100.;  // for atmosphere treatment or for adiabatic evolution

static double minlapse =1000.;

static double psimax = 0.0;

#ifdef ATMOSPHERE_GENERAL_RELATIVITY

static double rhomax_gr;
static double atm_fct = 1.e-7;
static double thr_fct = 1.e2;//org: 1.5; // for stronger damping 1.e2; org 1.5

#endif

#if METRIC_TYPE == 3
struct oned_metric   // this is used for a non-equidistant grid in radial direction
{
    int nsize;      // size of one-dimensional arrays
    double *r;    // radial coordiante measured from center, can be non-equidsitant
    double *alp;
    double *psi;
    double *bx;
    double *by;
    double *bz;
    //double mgrav;   // gravitational mass of the system to continue beyond the metric grid with Schwarzschild metric
    //-> we use the available BH variables for the gravitational mass and center of the metric
    
    // the numerical shift vector is actually not used because menaingless in 1d

    
};

// global pointer to 1d metric table
static struct oned_metric *oned_metric;

#endif

#if METRIC_TYPE == 4
struct oned_metric    // this is used for a equidistant grid in radial direction
{
    int nsize;      // size of one-dimensional arrays
    double dr;      // stepsize of radial grid
    double r0;      // smallest grid value, usually 0.
    double rmax;    // largest grid value
    double *alp;
    double *psi;

    //double mgrav;   // gravitational mass of the system to continue beyond the metric grid with Schwarzschild metric
    //-> we use the available BH variables for the gravitational mass and center of the metric

};

// global pointer to 1d metric table
static struct oned_metric *oned_metric;

#endif
void get_metric_general_relativity(double x, double y, double z, double *alp, double *psi, double *bx, double *by, double *bz)
{

    
#if METRIC_TYPE == 1

    // this routine returns the metric of a schwarzschild black hole in isotropic coordinates, i.e. cfc
    
    // possibly one could directly hand over the particle positions; use a black hole mass from the prepo to be sure t
    // that it's the same everywhere, and the same for BH position!!!!
    
    double mbh = BH_Mass;  // maybe use black hole mass from prepro
    
    // we assume that the BH is located 
    
    double xbh = BH_x;  // use preproc or read in BH position
    double ybh = BH_y;
    double zbh = BH_z;
    
    double r = sqrt( (x-xbh)*(x-xbh) + (y-ybh)*(y-ybh) + (z-zbh)*(z-zbh) );
    
    *psi  =  1. + 0.5 * mbh / r;
    
    *alp  =  ( 1. - 0.5 * mbh / r ) / *psi;
    
    *bx   = 0.0;
    *by   = 0.0;
    *bz   = 0.0;

    return;

#endif        
    
#if METRIC_TYPE == 2

    // this routine returns the metric of a Kerr black hole in isotropic coordinates, i.e. cfc
    // this implies that it's not a perfect Kerr BH, but very, very close
    // possibly one could directly hand over the particle positions; use a black hole mass from the prepo to be sure t
    // that it's the same everywhere, and the same for BH position!!!!
    
    // metric taken from Grandclement et al. 2002

    double mbh = BH_Mass;  // maybe use black hole mass from prepro
    
    double a = BH_Spin;

    // we assume that the BH is located 
    
    double xbh = BH_x;  // use preproc or read in BH position
    double ybh = BH_y;
    double zbh = BH_z;
    
    double r = sqrt( (x-xbh)*(x-xbh) + (y-ybh)*(y-ybh) + (z-zbh)*(z-zbh) );
    
    double phi = atan2(y-ybh,x-xbh);

    double theta = acos( (z-zbh) / r );

    double mu = (z-zbh) / r;

    double capr = r + 0.25 * (mbh*mbh-a*a) / r + mbh;

    double sigma = capr * capr + a*a*mu*mu;

    *alp = 1. - 2.*mbh*capr/sigma 
          + (4.*a*a*mbh*mbh*capr*capr*sin(theta)*sin(theta)) 
             / ( sigma*sigma*(capr*capr+a*a) + 2.*a*a*sigma*mbh*capr*sin(theta)*sin(theta) );

    *alp = sqrt(*alp);

    *psi = 1. + 2.*mbh/r + (3.*mbh*mbh+a*a*mu*mu) / (2.*r*r) 
              + (mbh*mbh-a*a)*mbh / (2.*r*r*r) 
              + (mbh*mbh-a*a)*(mbh*mbh-a*a) / (16.*r*r*r*r); // that's psi^4

    *psi = sqrt(*psi);
    *psi = sqrt(*psi);   // since the above formula gives psi^4

    double bphi = (2.*a*mbh*capr) 
                  / ( sigma*(capr*capr+a*a) + 2.*a*a*mbh*capr*sin(theta)*sin(theta) );
    
    *bx   = - bphi * sin(phi);
    *by   = + bphi * cos(phi);
    *bz   = 0.0;
     
    return;
#endif

#if METRIC_TYPE == 3

    double mbh = BH_Mass;  // this is read in !
    
    // we assume that the center is located 
    
    double xbh = BH_x;  // this is read in
    double ybh = BH_y;
    double zbh = BH_z;
    
    double r = sqrt( (x-xbh)*(x-xbh) + (y-ybh)*(y-ybh) + (z-zbh)*(z-zbh) );
    
    int nmax = oned_metric->nsize - 1;
    double rmin = oned_metric->r[0];
    double rmax = oned_metric->r[nmax];

    
    if ( r <= rmax && r >= rmin)  // we are inside the radial range provided by table
    {
        int n1 = 0;
        int n2 = nmax;
        
        int k;
        for(k = 0; k < oned_metric->nsize; k++)   // we don't need so many iterations but to be sure
        {
            int i = ( n2 - n1 ) / 2;
            
            if ( (n2 - n1) == 1)
               break;
            else if ( oned_metric->r[i] > r )
                n2 = i;
            else if ( oned_metric->r[i] <= r )
                n1 = i;
        }
        
        // now we found the right indices and can interpolate linearly
        
        double r1 = oned_metric->r[n1];
        double r2 = oned_metric->r[n2];       
        
        *alp  =  (oned_metric->alp[n2]-oned_metric->alp[n1]) / (r2-r1) *  (r-r1) + oned_metric->alp[n1];
        *psi  =  (oned_metric->psi[n2]-oned_metric->psi[n1]) / (r2-r1) *  (r-r1) + oned_metric->psi[n1];

          
    }
    if ( r < rmin )  // let's extrapolate to lower values
    {
        double r1 = oned_metric->r[0];
        double r2 = oned_metric->r[1];
        
        *alp =  (oned_metric->alp[1]-oned_metric->alp[0]) / (r2-r1) *  (r-r1) + oned_metric->alp[0];
        *psi =  (oned_metric->psi[1]-oned_metric->psi[0]) / (r2-r1) *  (r-r1) + oned_metric->psi[0];       
    }
    if ( r > rmax ) // we use the Schwarzschild metric in isotropic coordinates
    {               // but we should be careful that the transitoin is smooth
                    // one day we may introduce a check (todo)
    
    *psi  =  1. + 0.5 * mbh / r;
    
    *alp  =  ( 1. - 0.5 * mbh / r ) / *psi;
        
    }
    
    
    *bx   =  0.0;
    *by   =  0.0;   // todo actually it doens't make much sense to use numerical shift in 1d
    *bz   =  0.0;
    
     /*   if ( !(*alp>0 && *alp <1) )
        {
            printf("cenetr %g %g %g %g\n",xbh,ybh,zbh,mbh); 
            printf("alp %g %g %g %g %g\n",*alp,*psi,r,rmax,rmin);
            printf("\n");            
            int k;
            for(k = 0; k < oned_metric->nsize; k++)
            {
                if (k % 53){
                    printf("%d %g %g %g\n",k,&oned_metric->r[k],&oned_metric->alp[k],&oned_metric->psi[k]);
                }
            }
            printf("\n");
            terminate("gude");
        }*/

    
#endif

#if METRIC_TYPE == 4

    double mbh = BH_Mass;  // this is read in !
    
    // we assume that the center is located 
    
    double xbh = BH_x;  // this is read in
    double ybh = BH_y;
    double zbh = BH_z;
    
    double r = sqrt( (x-xbh)*(x-xbh) + (y-ybh)*(y-ybh) + (z-zbh)*(z-zbh) );
    
    int nmax    = oned_metric->nsize - 1;
    double rmin = oned_metric->r0;
    double dr   = oned_metric->dr;
    double rmax = oned_metric->rmax;

    // determine cell where our radial coordinate lies
            
    int k = floor( r / dr );
            
    if ( r < rmax && r >= rmin)  // we are inside the radial range provided by table
    {
        
        double r1 = (double) k  * dr ;
        
        *alp  =  (oned_metric->alp[k+1]-oned_metric->alp[k]) / dr *  (r-r1) + oned_metric->alp[k];
        *psi  =  (oned_metric->psi[k+1]-oned_metric->psi[k]) / dr *  (r-r1) + oned_metric->psi[k];
          
    }
    
    if ( r < rmin )  // let's extrapolate to lower values
    {
        *alp =  (oned_metric->alp[1]-oned_metric->alp[0]) / dr *  (r-rmin) + oned_metric->alp[0];
        *psi =  (oned_metric->psi[1]-oned_metric->psi[0]) / dr *  (r-rmin) + oned_metric->psi[0];       
    }
    if ( r > rmax ) // we use the Schwarzschild metric in isotropic coordinates
    {               // but we should be careful that the transitoin is smooth
                    // one day we may introduce a check (todo)
    
    *psi  =  1. + 0.5 * mbh / r;
    
    *alp  =  ( 1. - 0.5 * mbh / r ) / *psi;
        
    }
    
    
    *bx   =  0.0;
    *by   =  0.0;   // todo actually it doens't make much sense to use numerical shift in 1d
    *bz   =  0.0;
    
     /*   if ( !(*alp>0 && *alp <1) )
        {
            printf("cenetr %g %g %g %g\n",xbh,ybh,zbh,mbh); 
            printf("alp %g %g %g %g %g\n",*alp,*psi,r,rmax,rmin);
            printf("\n");            
            int k;
            for(k = 0; k < oned_metric->nsize; k++)
            {
                if (k % 53){
                    printf("%d %g %g %g\n",k,&oned_metric->r[k],&oned_metric->alp[k],&oned_metric->psi[k]);
                }
            }
            printf("\n");
            terminate("gude");
        }*/

    
#endif    
    
}


void get_metric_and_derivs_general_relativity(double x, double y, double z, double *alp, double *psi, double *bx, double *by, double *bz, double dalp[3], double dpsi[3], double dbx[3], double dby[3], double dbz[3])
{

  
#if METRIC_TYPE == 1
  
    // this routine returns the metric of a schwarzschild black hole in isotropic coordinates, i.e. cfc
    
    // possibly one could directly hand over the particle positions; use a black hole mass from the prepo to be sure t
    // that it's the same everywhere, and the same for BH position!!!!
    
    double mbh = BH_Mass;  // maybe use black hole mass from prepro
    
    // we assume that the BH is located 
    
    double xbh = BH_x;  // use preproc or read in BH position
    double ybh = BH_y;
    double zbh = BH_z;
    
    double r = sqrt( (x-xbh)*(x-xbh) + (y-ybh)*(y-ybh) + (z-zbh)*(z-zbh) );
    
    *psi  =  1. + 0.5 * mbh / r;
    
    *alp  =  ( 1. - 0.5 * mbh / r ) / *psi;
    
    *bx   = 0.0;
    *by   = 0.0;
    *bz   = 0.0;
    
    // now we compute the derivatives
    
    dpsi[0] = - 0.5 * mbh * (x-xbh) / pow(r,3.);
    dpsi[1] = - 0.5 * mbh * (y-ybh) / pow(r,3.);
    dpsi[2] = - 0.5 * mbh * (z-zbh) / pow(r,3.);    
    
    dalp[0] = - dpsi[0] * (1.+*alp) / *psi;  // check that this is correct
    dalp[1] = - dpsi[1] * (1.+*alp) / *psi;
    dalp[2] = - dpsi[2] * (1.+*alp) / *psi;
    
    int k;
    
    for(k = 0; k < 3; k++)  // note index gives direction of derivative
       dbx[k] = 0.0;
       
    for(k = 0; k < 3; k++)
       dby[k] = 0.0;

    for(k = 0; k < 3; k++)
       dbz[k] = 0.0;
       
    return;

#endif    
    
#if (METRIC_TYPE == 2) || (METRIC_TYPE == 3) || (METRIC_TYPE == 4)

    // Kerr black hole in cfc and its derivatives
    // since we are lasy and save a bit time, we compute deriavtives numerically

  double lalp, lpsi, lbx, lby, lbz;  // metric; local 

  /*get_metric_general_relativity(x,y,z,&lalp,&lpsi,&lbx,&lby,&lbz);
    
  *psi  = lpsi;
    
  *alp  = lalp;
    
  *bx   = lbx;
  *by   = lby;
  *bz   = lbz;*/

  double dx = BH_Mass * 1e-5;  // characteristic small scale for finite differcing
  double dy = BH_Mass * 1e-5;
  double dz = BH_Mass * 1e-5;

  double lalp2, lpsi2, lbx2, lby2, lbz2;
  
  get_metric_general_relativity(x-dx,y,z,&lalp,&lpsi,&lbx,&lby,&lbz);
  get_metric_general_relativity(x+dx,y,z,&lalp2,&lpsi2,&lbx2,&lby2,&lbz2);

  dpsi[0] = 0.5*(lpsi2 - lpsi) / dx;
  dalp[0] = 0.5*(lalp2 - lalp) / dx;
  dbx[0]  = 0.5*(lbx2  - lbx) / dx;
  dby[0]  = 0.5*(lby2  - lby) / dx;
  dbz[0]  = 0.5*(lbz2  - lbz) / dx;
  
  get_metric_general_relativity(x,y-dy,z,&lalp,&lpsi,&lbx,&lby,&lbz);
  get_metric_general_relativity(x,y+dy,z,&lalp2,&lpsi2,&lbx2,&lby2,&lbz2);

  dpsi[1] = 0.5*(lpsi2 - lpsi) / dy;
  dalp[1] = 0.5*(lalp2 - lalp) / dy;
  dbx[1]  = 0.5*(lbx2  - lbx) / dy;
  dby[1]  = 0.5*(lby2  - lby) / dy;
  dbz[1]  = 0.5*(lbz2  - lbz) / dy;

  get_metric_general_relativity(x,y,z-dz,&lalp,&lpsi,&lbx,&lby,&lbz);
  get_metric_general_relativity(x,y,z+dz,&lalp2,&lpsi2,&lbx2,&lby2,&lbz2);

  dpsi[2] = 0.5*(lpsi2 - lpsi) / dz;
  dalp[2] = 0.5*(lalp2 - lalp) / dz;
  dbx[2]  = 0.5*(lbx2  - lbx) / dz;
  dby[2]  = 0.5*(lby2  - lby) / dz;
  dbz[2]  = 0.5*(lbz2  - lbz) / dz;

  get_metric_general_relativity(x,y,z,&lalp,&lpsi,&lbx,&lby,&lbz);

  /*if ( !gsl_finite(dalp[0]) )
  {
      printf("metric deriv xyz %g %g %g %g %g\n",x,y,z,lalp,dalp[0]);
      terminate("metric deriv x");
  }
    if ( !gsl_finite(dalp[1]))
  {
      printf("metric deriv xyz %g %g %g %g %g\n",x,y,z,lalp,dalp[1]);
      terminate("metric deriv y");
  }
    if ( !gsl_finite(dalp[2]) )
  {
      printf("metric deriv xyz %g %g %g %g %g\n",x,y,z,lalp,dalp[2]);
      terminate("metric deriv z");
  }*/
  
  *psi  = lpsi;
  *alp  = lalp;
  *bx   = lbx;
  *by   = lby;
  *bz   = lbz;
  
#endif
        
}

void compute_source_terms_general_relativity(struct particle_data *P, struct sph_particle_data *SphP, int i, double *srcsx, double *srcsy, double *srcsz, double *srctau)
{
    
    /* this routine computes the general relativistic source term */

    //todo since we evolve absolute quantities and not density we need to multiply by the volume
    
    //this should maybe be done outside this routine or down but then I need to be careful where the volume really cancel where not
    
    
  // Position of Cells
  
  double x = SphP[i].Center[0];
  double y = SphP[i].Center[1];
  double z = SphP[i].Center[2];
    
  // metric quantities in cfc - just declaration
  
  double alp, psi, bx, by, bz;  // metric
  double dalp[3];  // derivatives of metric 
  double dpsi[3];
  double dbx[3];
  double dby[3];
  double dbz[3];
    
  get_metric_and_derivs_general_relativity(x, y, z, &alp, &psi, &bx, &by, &bz, dalp, dpsi, dbx, dby, dbz);
  
  if ( alp < minlapse)
      minlapse = alp;
  
  if ( psi > psimax)
      psimax = psi;
  
  double ps4 = pow(psi,4.);
  
  double det = pow(psi,6.);
  
  double shift[3];
  
  shift[0] = bx;
  shift[1] = by;
  shift[2] = bz;
  
  //ACHTUNG: nicht vergessen durch Volumen zu teilen !!!

  double sii[3]; // this shouldbeome a sij[3,3] for non diagonal metric
  
  int k;
  for(k = 0; k < 3; k++)
    sii[k]  = SphP[i].Momentum[k] / SphP[i].Volume * P[i].Vel[k] + det * SphP[i].Pressure;
  
  *srcsx = 2. * alp * dpsi[0] / psi * ( sii[0] + sii[1] + sii[2] )
          +       SphP[i].Momentum[0] / SphP[i].Volume * dbx[0]
          +       SphP[i].Momentum[1] / SphP[i].Volume * dby[0] 
          +       SphP[i].Momentum[2] / SphP[i].Volume * dbz[0]
          - ( SphP[i].Energy / SphP[i].Volume + P[i].Mass / SphP[i].Volume ) * dalp[0];
  
  *srcsy = 2. * alp * dpsi[1] / psi * ( sii[0] + sii[1] + sii[2] )
          +       SphP[i].Momentum[0] / SphP[i].Volume * dbx[1]
          +       SphP[i].Momentum[1] / SphP[i].Volume * dby[1] 
          +       SphP[i].Momentum[2] / SphP[i].Volume * dbz[1]
          - ( SphP[i].Energy / SphP[i].Volume + P[i].Mass / SphP[i].Volume ) * dalp[1]  ;
  
  *srcsz = 2. * alp * dpsi[2] / psi * ( sii[0] + sii[1] + sii[2] )
          +       SphP[i].Momentum[0] / SphP[i].Volume * dbx[2]
          +       SphP[i].Momentum[1] / SphP[i].Volume * dby[2] 
          +       SphP[i].Momentum[2] / SphP[i].Volume * dbz[2]
          - ( SphP[i].Energy / SphP[i].Volume + P[i].Mass / SphP[i].Volume ) * dalp[2]  ;
    
       double Kxx, Kxy, Kxz, Kyy, Kyz, Kzz;
          
       Kxx  = 0.50 * ( dbx[0] + dbx[0] - 2.0/3.0 * (dbx[0] + dby[1] + dbz[2]) );// / alp(i,j,k)

       Kxy  = 0.50 * ( dbx[1] + dby[0] );// / alp(i,j,k)

       Kxz  = 0.50 * ( dbx[2] + dbz[0] );// / alp(i,j,k)

       Kyy  = 0.50 * ( dby[1] + dby[1] - 2.0/3.0 * (dbx[0] + dby[1] + dbz[2]) );// / alp(i,j,k)

       Kyz  = 0.50 * ( dbz[1] + dby[2] );// / alp(i,j,k)

       Kzz  = 0.50 * ( dbz[2] + dbz[2] - 2.0/3.0 * (dbx[0] + dby[1] + dbz[2]) );// / alp(i,j,k)        
          
    double Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz;
          
      Sxx  = sii[0];
      Sxy  = SphP[i].Momentum[0] / SphP[i].Volume * P[i].Vel[1];   //sx3(i,j,k)*vy3(i,j,k)
      Sxz  = SphP[i].Momentum[0] / SphP[i].Volume * P[i].Vel[2];  //sx3(i,j,k)*vz3(i,j,k)

      Syx  = SphP[i].Momentum[1] / SphP[i].Volume * P[i].Vel[0];  //sy3(i,j,k)*vx3(i,j,k)
      Syy  = sii[1];
      Syz  = SphP[i].Momentum[1] / SphP[i].Volume * P[i].Vel[2];   //sy3(i,j,k)*vz3(i,j,k)

      Szx  = SphP[i].Momentum[2] / SphP[i].Volume * P[i].Vel[0];  //sz3(i,j,k)*vx3(i,j,k)
      Szy  = SphP[i].Momentum[2] / SphP[i].Volume * P[i].Vel[1];  //sz3(i,j,k)*vy3(i,j,k)       
      Szz  = sii[2];    
  
   *srctau =           Sxx * Kxx 
                    + Sxy * Kxy 
                    + Sxz * Kxz 
                    + Syx * Kxy 
                    + Syy * Kyy 
                    + Syz * Kyz 
                    + Szx * Kxz 
                    + Szy * Kyz 
                    + Szz * Kzz  
                    - SphP[i].Momentum[0] / SphP[i].Volume / ps4 * dalp[0] 
                    - SphP[i].Momentum[1] / SphP[i].Volume / ps4 * dalp[1]  
                    - SphP[i].Momentum[2] / SphP[i].Volume / ps4 * dalp[2];   
      
    
  // here we need to be careful when the spatial metric becomes nondiagonal
  
}

#if ADIABATIC_GENERAL_RELATIVITY==1
void update_primitive_variables_general_relativity(struct particle_data *P, struct sph_particle_data *SphP, int i, struct pv_update_data *pvd)
{
  /* here we assume that the EoS is available in this routine either as table or analytically
   for the momemt this routine works only for analytic ideal gas EoS

   general relativstic with polytrope/ideal fluid

   -> simple expressions and method for root finding*/

  
  // Position of Cells
  
  double x = SphP[i].Center[0];
  double y = SphP[i].Center[1];
  double z = SphP[i].Center[2];
    
  // metric quantities in cfc - just declaration
  
  double alp, psi, bx, by, bz;
    
  get_metric_general_relativity(x,y,z,&alp,&psi,&bx,&by,&bz);
   
  double det = pow(psi,6.);

  double ps4 = pow(psi,4.);

  double d = P[i].Mass / SphP[i].Volume / det;
  double sx = SphP[i].Momentum[0] / SphP[i].Volume / det;
  double sy = SphP[i].Momentum[1] / SphP[i].Volume / det;
  double sz = SphP[i].Momentum[2] / SphP[i].Volume / det;
  double t = SphP[i].Energy / SphP[i].Volume / det;

  double p0 = SphP[i].Pressure;

  double     r,p,e,h;         // primitives: rho, pressure, vxyz, internal specific energy
  double     kpoly, gam;         // polytropic constant and index for analytic ideal gas eos

 // internal variables

  int                    iter;               // iterations   
  int                    nmax;            // maximum number of iterations
  double                 s2, tpd, root;     // inetrnal variables for con2primSR
  double                 dr, de;            // derivatives of pressure wrt r, e
  double                 f, df;             // Newton method's function and derivative
  double                 rel, prec;         // relative change in pressure and required precisio

  // EoS settings for analytic Eos only - maybe better taken from some module

  kpoly   =  kpoly_gr;   // <- for the atmosphere treatment

  gam     =  GAMMA;

  // Newton method settings

  nmax = 25;      // <- probably 10 is sufficient

  prec = 1e-10;   // <- may need to be reduce for certain problems, e.g. 1e-8 for sieglers second problem; org 1e-10
  
  double rhofloor = 1.e7 / 6.176e17;   // <- set for the atmosphere treatment - possibly overwritten by rhomax*factior

#ifdef ATMOSPHERE_GENERAL_RELATIVITY

  double pfloor = kpoly * pow(rhofloor,gam);

#if ATMOSPHERE_GENERAL_RELATIVITY == 1
  pfloor = kpoly * pow(rhomax_gr * atm_fct, gam);
  
  rhofloor = rhomax_gr * atm_fct;
#endif
  
#else
  
  double pfloor = 1.e-8;    // <- adapt e.g. for shocktube tests
    
  pfloor = kpoly * pow(1.e8/6.e17,gam);
  
#endif
  
  p  = p0;       // start value of pressure is its previous value

  s2 =  sx * sx + sy * sy + sz * sz;  

  s2 = s2 / ps4;  // because s2 = gamma^ij S_i S_j
  
  // here we will need to be careful when the spatial metric is not diagonal any longer
  
  // estimate rho via p0, could be done also via dd

//  r   = d;

// p = kpoly * pow(r,gam);
 
  r =  pow( p / kpoly, 1./gam);
 
  h = 1. + gam * p / (gam-1.) / r; 

  for (iter = 0; iter < nmax; iter++) {

    double wbar = sqrt( 1. + s2 / (d*h) / (d*h) );

    f    =  r * wbar - d;

    double hprime = kpoly * gam * pow(r,gam-2.) ;

    df   =  wbar - r * s2 * hprime / (wbar*d*d*h*h*h);
    
    // newton step

    rel = fabs(f/df/p);

    r   = r - f/df ;

    //if ( !(r>rhofloor) ) {
    //  r=rhofloor;
    //}
    if ( r < rhofloor){
        rel = 0.1 *prec; // just to exit loop
        r = 0.1*rhofloor; // to enter atmosphere treatment below
    }
    
    p   = kpoly * pow(r, gam);

    h   = 1. + gam * p / (gam-1.) / r; 

    if ( rel < prec ) break;  //  ! required accuracy reached

  }

  /*if ( iter == nmax ) {
    printf( "GENERAL_RELATIVITY: Max iter reached in in con2prim Newton method, rel=%g, prec=%g, f=%g, df=%g, p=%g, nmax=%i, p0=%g, i=%i\n", rel, prec, f, df, p, nmax, p0, i );
    printf( "GENERAL_RELATIVITY: cons: d=%g, s2=%g, t=%g, p0=%g\n",d,s2,t,p0 );
    print_particle_info(i);
    terminate( "general relativity update prim fail" );
  }*/

  // compute all variables based on current pressure a last time

  double w = d / r;

  SphP[i].Density = r ;

  SphP[i].Utherm = p / r / (gam-1.);

  SphP[i].Pressure = p;

  P[i].Vel[0] = sx / d / ps4 / w / h;  // also other numerical expression are possible here
  P[i].Vel[1] = sy / d / ps4 / w / h;  // this is v^i, the latter /ps4 to go from v_j to v^i
  P[i].Vel[2] = sz / d / ps4 / w / h;  // note that the sx ... are divided by the det already above

  // here we will need to be careful when the spatial metric is not diagonal any longer

  // since we evolve adiabatically, we reset tau - one could also let it evolve freely but then it may affect the evolution
  // if it appears in the fluxes or source terms

  SphP[i].Energy = det * ( d * (w * h - 1.) - p ) * SphP[i].Volume; 

  
#ifdef ATMOSPHERE_GENERAL_RELATIVITY  
  
  // simplistic atmosphere treatment

  //if ( !(SphP[i].Density >= rhofloor))
  //if ( !(SphP[i].Density >= 1.e-7))
  
#if ATMOSPHERE_GENERAL_RELATIVITY == 1  
  if ( !gsl_finite( SphP[i].Density ) || (SphP[i].Density < rhomax_gr * atm_fct * thr_fct))
  {
          SphP[i].Density     = rhomax_gr * atm_fct;
      
#else      
  if ( !gsl_finite( SphP[i].Density ) || (SphP[i].Density < rhofloor))
  {
    SphP[i].Density     = rhofloor;
    
#endif
    
    SphP[i].Pressure    = kpoly * pow(SphP[i].Density, gam);
    
    SphP[i].Utherm      = SphP[i].Pressure / (gam-1.) / SphP[i].Density;
    
    P[i].Vel[0] = 0.0;
    P[i].Vel[0] = 0.0; 
    P[i].Vel[0] = 0.0;
      
    // probably we should also update the conserved variables
    
    P[i].Mass = SphP[i].Density * det * SphP[i].Volume;
    SphP[i].Momentum[0] = 0.0;
    SphP[i].Momentum[1] = 0.0; 
    SphP[i].Momentum[2] = 0.0; 
    
    SphP[i].Energy = SphP[i].Density * SphP[i].Utherm * det * SphP[i].Volume;  // is this correct for atmosphere?
    
    
  } 
  
  
#endif  
  
  if(SphP[i].Density < 1.e-18)
    {
      printf("zeug: %g %g %g %g \n",r,p,d,w);
      printf("zeug: %g \n",P[i].Mass);
      printf("callcounter %d\n",GR_counter);
      print_particle_info(i);
      terminate( "shit in adaiabtic con2prim" );
    }
    
    
  // TODO: add scalars, ...
}
#elif ADIABATIC_GENERAL_RELATIVITY==2

// just another way of implementing it

void update_primitive_variables_general_relativity(struct particle_data *P, struct sph_particle_data *SphP, int i, struct pv_update_data *pvd)
{
  /* here we assume that the EoS is available in this routine either as table or analytically
   for the momemt this routine works only for analytic ideal gas EoS

   general relativstic with polytrope/ideal fluid

   -> simple expressions and method for root finding*/

  
  // Position of Cells
  
  double x = SphP[i].Center[0];
  double y = SphP[i].Center[1];
  double z = SphP[i].Center[2];
    
  // metric quantities in cfc - just declaration
  
  double alp, psi, bx, by, bz;
    
  get_metric_general_relativity(x,y,z,&alp,&psi,&bx,&by,&bz);
   
  double det = pow(psi,6.);

  double ps4 = pow(psi,4.);

  double d = P[i].Mass / SphP[i].Volume / det;
  double sx = SphP[i].Momentum[0] / SphP[i].Volume / det;
  double sy = SphP[i].Momentum[1] / SphP[i].Volume / det;
  double sz = SphP[i].Momentum[2] / SphP[i].Volume / det;
  double t = SphP[i].Energy / SphP[i].Volume / det;

  double p0 = SphP[i].Pressure;

  double     r,p,e,h,wbar;     // primitives: rho, pressure, vxyz, internal specific energy
  double     kpoly, gam;         // polytropic constant and index for analytic ideal gas eos

 // internal variables

  int                    iter;               // iterations   
  int                    nmax;            // maximum number of iterations
  double                 s2, tpd, root;     // inetrnal variables for con2primSR
  double                 dr, de;            // derivatives of pressure wrt r, e
  double                 f, df;             // Newton method's function and derivative
  double                 rel, prec;         // relative change in pressure and required precisio

  // EoS settings for analytic Eos only - maybe better taken from some module

  kpoly   =  kpoly_gr;   // <- for the atmosphere treatment

  gam     =  GAMMA;

  // Newton method settings

  nmax = 25;      // <- probably 10 is sufficient

  prec = 1e-10;   // <- may need to be reduce for certain problems, e.g. 1e-8 for sieglers second problem; org 1e-10
  
  double rhofloor = 1.e7 / 6.176e17;   // <- set for the atmosphere treatment - possibly overwriten by rhomax*factior

#ifdef ATMOSPHERE_GENERAL_RELATIVITY

  double pfloor = kpoly * pow(rhofloor,gam);

#if ATMOSPHERE_GENERAL_RELATIVITY == 1
  pfloor = kpoly * pow(rhomax_gr * atm_fct, gam);
  rhofloor = rhomax_gr * atm_fct;
#endif
  
#else
  
  double pfloor = 1.e-8;    // <- adapt e.g. for shocktube tests
    
  pfloor = kpoly * pow(1.e8/6.e17,gam);
  
#endif
  
  p  = p0;       // start value of pressure is its previous value

  s2 =  sx * sx + sy * sy + sz * sz;  

  s2 = s2 / ps4;  // because s2 = gamma^ij S_i S_j
  
  // here we will need to be careful when the spatial metric is not diagonal any longer
  
  r = d;  // as a first guess


  for (iter = 0; iter < nmax; iter++) {

    p = kpoly * pow(r,gam);
    
    e = p / (gam-1.) / r;
    
    h = 1. + e + p/r;
    
    wbar = 1. + s2 / (d*d*h*h);

    wbar = sqrt(wbar);
    
    f = r * wbar - d;
    
    double hprime = kpoly * gam * pow(r,gam-1) / r ;
    
    df = wbar - r * s2 * hprime / (wbar*d*d*h*h*h);
    
    rel=fabs(f/df/r);
    
    r = r -f/df;
    
    if ( r < rhofloor){
        rel = 0.1 *prec; // just to exit loop
        r = 0.1*rhofloor; // to enter atmosphere treatment below
    }
    
    if ( rel < prec ) break;  //  ! required accuracy reached
  }

  /*if ( iter == nmax ) {
    printf( "GENERAL_RELATIVITY: Max iter reached in in con2prim Newton method, rel=%g, prec=%g, f=%g, df=%g, p=%g, nmax=%i, p0=%g, i=%i\n", rel, prec, f, df, p, nmax, p0, i );
    printf( "GENERAL_RELATIVITY: cons: d=%g, s2=%g, t=%g, p0=%g\n",d,s2,t,p0 );
    print_particle_info(i);
    terminate( "general relativity update prim fail" );
  }*/

  if ( !gsl_finite( r ) || (r < rhomax_gr * atm_fct * thr_fct))
  {
          SphP[i].Density     = rhomax_gr * atm_fct; 
          r = SphP[i].Density;
          SphP[i].Pressure = kpoly*pow(r,gam);
          SphP[i].Utherm = SphP[i].Pressure / r / (gam-1.);
          P[i].Vel[0] = 0.0;
          P[i].Vel[1] = 0.0;
          P[i].Vel[2] = 0.0;

      // reset the conserved quantities
              
          P[i].Mass = SphP[i].Density * det * SphP[i].Volume;
          SphP[i].Momentum[0] = 0.0;
          SphP[i].Momentum[1] = 0.0; 
          SphP[i].Momentum[2] = 0.0; 
    
          SphP[i].Energy = SphP[i].Density * SphP[i].Utherm * det * SphP[i].Volume; 
  }        
  else    // compute all variables based on current pressure a last time
  {
        SphP[i].Density = r;    
        SphP[i].Pressure = kpoly*pow(r,gam);    
        SphP[i].Utherm = SphP[i].Pressure / r / (gam-1.);  
        h = 1. + SphP[i].Utherm + SphP[i].Pressure/ SphP[i].Density;
    
        wbar = 1. + s2 / (d*d*h*h);
        
        wbar = sqrt(wbar);
  
        P[i].Vel[0] = sx / d / wbar / h /ps4;
        P[i].Vel[1] = sy / d / wbar / h /ps4;       
        P[i].Vel[2] = sz / d / wbar / h /ps4;

        SphP[i].Energy = det * (d*(h*wbar-1.) -SphP[i].Pressure ) * SphP[i].Volume; 
                  
  }
  
 /* if(p < 0)
    {
      print_particle_info(i);
      terminate( "shit" );
    }*/
    
    
  // TODO: add scalars, ...
}



#else

void update_primitive_variables_general_relativity_org(struct particle_data *P, struct sph_particle_data *SphP, int i, struct pv_update_data *pvd)
{
  /* here we assume that the EoS is available in this routine either as table or analytically
   for the momemt this routine works only for analytic ideal gas EoS

   general relativstic with polytrope/ideal fluid

   -> simple expressions and method for root finding*/

  
  // Position of Cells
  
  double x = SphP[i].Center[0];
  double y = SphP[i].Center[1];
  double z = SphP[i].Center[2];
    
  // metric quantities in cfc - just declaration
  
  double alp, psi, bx, by, bz;
    
  get_metric_general_relativity(x,y,z,&alp,&psi,&bx,&by,&bz);
   
  double det = pow(psi,6.);

  double ps4 = pow(psi,4.);

  double d = P[i].Mass / SphP[i].Volume / det;
  double sx = SphP[i].Momentum[0] / SphP[i].Volume / det;
  double sy = SphP[i].Momentum[1] / SphP[i].Volume / det;
  double sz = SphP[i].Momentum[2] / SphP[i].Volume / det;
  double t = SphP[i].Energy / SphP[i].Volume / det;

  double p0 = SphP[i].Pressure;

  double     r,p,e;     // primitives: rho, pressure, vxyz, internal specific energy
  double     kpoly, gam;         // polytropic constant and index for analytic ideal gas eos

 // internal variables

  int                    iter;               // iterations   
  int                    nmax;            // maximum number of iterations
  double                 s2, tpd, root;     // inetrnal variables for con2primSR
  double                 dr, de;            // derivatives of pressure wrt r, e
  double                 f, df;             // Newton method's function and derivative
  double                 rel, prec;         // relative change in pressure and required precisio

  // EoS settings for analytic Eos only - maybe better taken from some module

  kpoly   =  kpoly_gr;   // <- for the atmosphere treatment

  gam     =  GAMMA;

  // Newton method settings

  nmax = 25;      // <- probably 10 is sufficient

  prec = 1.e-8;//1e-10;   // <- may need to be reduce for certain problems, e.g. 1e-8 for sieglers second problem; org 1e-10
  
  double rhofloor = 1.e7 / 6.176e17;   // <- set for the atmosphere treatment - possibly overwriten by rhomax*factior

#ifdef ATMOSPHERE_GENERAL_RELATIVITY

  double pfloor = kpoly * pow(rhofloor,gam);

#if ATMOSPHERE_GENERAL_RELATIVITY == 1
  pfloor = kpoly * pow(rhomax_gr * atm_fct, gam);
#endif
  
#else
  
  double pfloor = 1.e-8;    // <- adapt e.g. for shocktube tests
    
  pfloor = kpoly * pow(1.e8/6.e17,gam);
  
#endif
  
  p  = p0;       // start value of pressure is its previous value

  s2 =  sx * sx + sy * sy + sz * sz;  

  s2 = s2 / ps4;  // because s2 = gamma^ij S_i S_j
  
  // here we will need to be careful when the spatial metric is not diagonal any longer
  
  
  
  // do some checks: -- maybe one should monitor how often this is imposed

  if ( t < 0.0) { // reset to the previous value by recalculating cons from prim

    double w = 1. / sqrt(1. - (P[i].Vel[0]*P[i].Vel[0] + P[i].Vel[1]*P[i].Vel[1] + P[i].Vel[2]*P[i].Vel[2] ) ); 

    double oldd = w * SphP[i].Density;

    double oldp = SphP[i].Pressure;

    double oldu = SphP[i].Utherm;

    double oldh = 1. + oldu + oldp/SphP[i].Density;

    t =  oldd * ( oldh * w -1. ) - oldp; //thats a bit more

    //printf("tau reset");

    terminate("tau was reset");

    }


  for (iter = 0; iter < nmax; iter++) {

    tpd  =  t + p + d;

    root =  sqrt(tpd*tpd-s2);

    r    =  d / tpd * root ;

    e    =  (root-p*tpd/root-d) / d;

    dr   =  d * s2 / (root*tpd*tpd);

    de   =  p * s2 / (r*root*root*tpd);

    f    =  p - (gam-1.0)*r*e;

    df   =  1.0 - (gam-1.0)*e*dr - (gam-1.0)*r*de;
    
    // newton step
    rel = fabs(f/df/p);

    p   = p - f/df ;

    if ( !(p>pfloor) ) {
      p=pfloor;
    }


   // if ( !(p>0.0) ) {
   //   p = pfloor;
   //  }

    if ( !(p > 0.0) ) {
      printf("  ");
      printf("p<0;p=%g,r=%g,e=%g,f=%g,p0=%g,d=%g,s2=%g  ",p,r,e,f,p0,d,s2);
      printf("t=%g,tpd=%g,root=%g,rel=%g",t,tpd,root,rel);
      printf("  ");
      terminate("Scheisse");
      }


    if ( rel < prec ) break;  //  ! required accuracy reached
  }

  /*if ( iter == nmax ) {
    printf( "GENERAL_RELATIVITY: Max iter reached in in con2prim Newton method, rel=%g, prec=%g, f=%g, df=%g, p=%g, nmax=%i, p0=%g, i=%i\n", rel, prec, f, df, p, nmax, p0, i );
    printf( "GENERAL_RELATIVITY: cons: d=%g, s2=%g, t=%g, p0=%g\n",d,s2,t,p0 );
    print_particle_info(i);
    terminate( "general relativity update prim fail" );
    }*/

  // compute all variables based on current pressure a last time

  tpd  =  t + p + d;

  root =  sqrt(tpd*tpd-s2);

  SphP[i].Density = d / tpd * root ;

  SphP[i].Utherm = (root-p*tpd/root-d) / d;
  SphP[i].Pressure = p;

  P[i].Vel[0] = sx / tpd / ps4;  // also other numerical expression are possible here
  P[i].Vel[1] = sy / tpd / ps4;  // this is v^i, the latter /ps4 to go from v_j to v^i
  P[i].Vel[2] = sz / tpd / ps4;  // note that the sx ... are divided by the det already above
  
  // here we will need to be careful when the spatial metric is not diagonal any longer

  
#ifdef ATMOSPHERE_GENERAL_RELATIVITY  
  
  // simplistic atmosphere treatment

  //if ( !(SphP[i].Density >= rhofloor))
  //if ( !(SphP[i].Density >= 1.e-7))
  
#if ATMOSPHERE_GENERAL_RELATIVITY == 1  
  if ( !gsl_finite( SphP[i].Density ) || (SphP[i].Density < rhomax_gr * atm_fct * thr_fct))
  {
          SphP[i].Density     = rhomax_gr * atm_fct;
      
#else      
  if ( !gsl_finite( SphP[i].Density ) || (SphP[i].Density < rhofloor))
  {
    SphP[i].Density     = rhofloor;
    
#endif
    
    SphP[i].Pressure    = kpoly * pow(SphP[i].Density, gam);
    
    SphP[i].Utherm      = SphP[i].Pressure / (gam-1.) / SphP[i].Density;
    
    P[i].Vel[0] = 0.0;
    P[i].Vel[0] = 0.0; 
    P[i].Vel[0] = 0.0;
      
    // probably we should also update the conserved variables
    
    P[i].Mass = SphP[i].Density * det * SphP[i].Volume;
    SphP[i].Momentum[0] = 0.0;
    SphP[i].Momentum[1] = 0.0; 
    SphP[i].Momentum[2] = 0.0; 
    
    SphP[i].Energy = SphP[i].Density * SphP[i].Utherm * det * SphP[i].Volume;  // is this correct for atmosphere?
    
    
  } 
  
  
#endif  
  
  if(p < 0)
    {
      print_particle_info(i);
      terminate( "shit" );
    }
    
    
  // TODO: add scalars, ...
}


// that's a new implementation for hot evolution closer to my current grid version

void update_primitive_variables_general_relativity(struct particle_data *P, struct sph_particle_data *SphP, int i, struct pv_update_data *pvd)
{
  /* here we assume that the EoS is available in this routine either as table or analytically
   for the momemt this routine works only for analytic ideal gas EoS

   general relativstic with polytrope/ideal fluid

   -> simple expressions and method for root finding*/

  // that's the version that I also use for dirty micphys eoss
  
  // Position of Cells
  
  double x = SphP[i].Center[0];
  double y = SphP[i].Center[1];
  double z = SphP[i].Center[2];
    
  // metric quantities in cfc - just declaration
  
  double alp, psi, bx, by, bz;
    
  get_metric_general_relativity(x,y,z,&alp,&psi,&bx,&by,&bz);
   
  double det = pow(psi,6.);

  double ps4 = pow(psi,4.);

  double dd = P[i].Mass / SphP[i].Volume / det;
  double ssx = SphP[i].Momentum[0] / SphP[i].Volume / det;
  double ssy = SphP[i].Momentum[1] / SphP[i].Volume / det;
  double ssz = SphP[i].Momentum[2] / SphP[i].Volume / det;
  double tt = SphP[i].Energy / SphP[i].Volume / det;

  if (tt<0.0){
      printf("tau smaller 0 %g \n",tt);
      terminate("because of tau");
      }
  
  double rstout,sxout,yeout,tauout,syout,szout;
  
   rstout = dd * det;
   sxout  = ssx * det;
   syout  = ssy * det;
   szout  = ssz * det;
   tauout = tt * det;
  
  double p0 = SphP[i].Pressure;

  double p = p0;
  
  double     r,e,peos,vx,vy,vz;     // primitives: rho, pressure, vxyz, internal specific energy
  double     kpoly, gam, ent;         // polytropic constant and index for analytic ideal gas eos
  
 // internal variables

  int                    iter;               // iterations   
  int                    nmax;            // maximum number of iterations
  double                 s2, tpd, root;     // inetrnal variables for con2primSR
  double                 dr, de;            // derivatives of pressure wrt r, e
  double                 f, df;             // Newton method's function and derivative
  double                 rel, prec;         // relative change in pressure and required precisio

  // EoS settings for analytic Eos only - maybe better taken from some module

  kpoly   =  kpoly_gr;   // <- for the atmosphere treatment

  gam     =  GAMMA;

  // Newton method settings

  nmax = 25;      // <- probably 10 is sufficient

  prec = 1e-10;   // <- may need to be reduce for certain problems, e.g. 1e-8 for sieglers second problem; org 1e-10
  
  
  if ( dd <= rhomax_gr*atm_fct*thr_fct ) {

    r  =  0.1 * rhomax_gr*atm_fct*thr_fct; // i just set a small value that it will enter the atmosphere treatment if clause

    goto ATMLABEL;  // sin, ok, but that's how i did it - to do !!
  }


  
  //if imposewmax==2

 //  check that maximum lorentz factor is bound if root negative

    s2   =  (ssx * ssx + ssy * ssy + ssz * ssz) / ps4;
    tpd  =  tt + p + dd;
    root =  (tpd*tpd-s2);

  if ( !(root > 0.0) && dd < 1.e-5 ) {

    double wmax=2.;
      
    tt = sqrt(s2) * wmax / sqrt( wmax*wmax -1.0 ) - p - dd;

    tauout = tt * det;

    s2   =  (ssx * ssx + ssy * ssy + ssz * ssz) / ps4;
    tpd  =  tt + p + dd;
    root =  (tpd*tpd-s2);

    if ( !(root >= 0.0) ) {

      printf("root still neg. %g %g \n",root,sqrt(1.0 - wmax*wmax));
      printf("%g %g %g %g %g \n",tpd,sqrt(s2),tt, p , dd);
      printf("wmax %g\n",wmax);
      print_particle_info(i);
      terminate("in wamx \n");

    }

  }

 //endif
  
  // start actual loop

  for (iter = 0; iter < nmax; iter++) {
      
    s2   =  (ssx * ssx + ssy * ssy + ssz * ssz) / ps4;
    tpd  =  tt + p + dd;
    root =  sqrt(tpd*tpd-s2);
    
    r    =  dd / tpd * root ;

    e    =  (root-p*tpd/root-dd) / dd;

    dr   =  dd * s2 / (root*tpd*tpd);

    de   =  p * s2 / (r*root*root*tpd);
    
    // call eos:
    
    peos = (gam-1.) * r * e;

//if ( peos < kpoly_gr*pow(r,GAMMA) )
//  peos = kpoly_gr*pow(r,GAMMA); // new as attempt

    f    =  p - peos;

    //determine derivatives numerically! with second eos call (not nice here but valid for general eos)
    double drho = 1.001;  // auxilary to determine deriavteive
    double deps = 1.001;
    
    double p2eos = (gam-1.) * r*drho * e;

//if ( p2eos < kpoly_gr*pow(r,GAMMA) )
//  p2eos = kpoly_gr*pow(r,GAMMA); // new as attempt

    double drhoeos = (p2eos-peos)/((drho-1.0)*r);
          
    p2eos = (gam-1.) * r * e*deps;  

//if ( p2eos < kpoly_gr*pow(r,GAMMA) )
//  p2eos = kpoly_gr*pow(r,GAMMA); // new as attempt
    
    double depseos = (p2eos-peos)/((deps-1.0)*e);
    
    df    =  1.0 -  drhoeos * dr - depseos * de;
    
    //newton step
    
    rel = fabs(f/df/p);

    p   = p - f/df;

    if ( p < 0.0 ) goto ATMLABEL;
    
    if ( rel < prec ) break;
    
  }
        
  // compute all variables based on current pressure a last time

    s2   =  (ssx * ssx + ssy * ssy + ssz * ssz) / ps4;
    tpd  =  tt + p + dd;
    root =  sqrt(tpd*tpd-s2);

    r    =  dd / tpd * root ;

    e    =  (root-p*tpd/root-dd) / dd;
  
    vx    =  ssx / tpd / ps4;        

    vy    =  ssy / tpd / ps4;
    
    vz    =  ssz / tpd / ps4;
    
  // do atmosphere treatment  
ATMLABEL: if ( p < 0.0 || r < rhomax_gr*atm_fct*thr_fct ) {

    r   =  rhomax_gr*atm_fct;
    vx  =  0.0; 
    vy  =  0.0;
    vz  =  0.0;
    
    p   = kpoly *  pow(r,gam);
    
    e   = p/r/(gam-1.);
        
    ent = 1.0 + e + p/r;
    
    rstout = r * det;
    sxout  = 0.0;
    syout  = 0.0;
    szout  = 0.0;
    tauout = r * det *ent - r * det - p * det;
    
    P[i].Mass = r * det * SphP[i].Volume;
    SphP[i].Momentum[0] = 0.0;
    SphP[i].Momentum[1] = 0.0; 
    SphP[i].Momentum[2] = 0.0; 
      
    SphP[i].Energy = P[i].Mass * ent - det * p * SphP[i].Volume - P[i].Mass;

    }
    
  // here actually con2primcheck, currently not used
    
    

  /*if ( iter == nmax ) {
    printf( "GENERAL_RELATIVITY: Max iter reached in in con2prim Newton method, rel=%g, prec=%g, f=%g, df=%g, p=%g, nmax=%i, p0=%g, i=%i\n", rel, prec, f, df, p, nmax, p0, i );
    printf( "GENERAL_RELATIVITY: cons: d=%g, s2=%g, t=%g, p0=%g\n",d,s2,t,p0 );
    print_particle_info(i);
    terminate( "general relativity update prim fail" );
  }*/


  SphP[i].Density  = r;

  SphP[i].Utherm   = e;
  SphP[i].Pressure = p;

  P[i].Vel[0] = vx;
  P[i].Vel[1] = vy;
  P[i].Vel[2] = vz;

      
    // probably we should also update the conserved variables
    
    //  and i could recalculate the cons based on prims
  
  if(p < 1.e-20) //org 0
    {
      print_particle_info(i);
      terminate( "shit in my con2prim" );
    }
    
    
  // TODO: add scalars, ...
}



#endif

double get_cfl_sound_speed_general_relativity(int p)
{
  double vx = P[p].Vel[0];
  double vy = P[p].Vel[1];
  double vz = P[p].Vel[2];
  
  
         
    //todo - correct postitiosn?  check where this is called, probably only in the cells and not on the faces
    
  double x = SphP[p].Center[0];
  double y = SphP[p].Center[1];
  double z = SphP[p].Center[2];
    
  // metric quantities in cfc - just declaration
  
  double alp, psi, bx, by, bz;

  get_metric_general_relativity(x,y,z,&alp,&psi,&bx,&by,&bz); 
  
  double ps4 = pow(psi,4.);
  
  //double bx_turned_to_xface = todo
  
  
  
  double v2 = vx*vx + vy*vy + vz*vz;
  
  v2 = v2 * ps4;
  
  double Enthalpy = 1. + SphP[p].Utherm + SphP[p].Pressure / SphP[p].Density;
  double csnd = sqrt( GAMMA * SphP[p].Pressure / (SphP[p].Density * Enthalpy) );
  double csnd2 = csnd * csnd;
  
  double lambdamax = 0;
  
  double shift[3];
  
  shift[0] = bx;  // todo: this is the the lab frame, so bx is really bx, right? ja
  shift[1] = by;
  shift[2] = bz;
  
  int k;
  for(k = 0; k < 3; k++)
    {
      double lambda0 = fabs( alp * P[p].Vel[k] - shift[k] );
      
      if(lambda0 > lambdamax)
        lambdamax = lambda0;
  
      double t1 = alp * P[p].Vel[k] * (1.-csnd2) / (1.-v2*csnd2);
  
      double t2 = alp * csnd / (1.-v2*csnd2);  // this is just the prefactor
  
      double gammaii = 1./ps4;  // this is gamma^ii; for cfc quite simple, this needs to be adapted if beyond cfc
  
      double t3 = sqrt( (1.-v2) * ( gammaii * (1.-v2*csnd2) - P[p].Vel[k] * P[p].Vel[k] * (1.-csnd2) ) );
  
      double lambdap = fabs( t1 + t2 * t3 - shift[k]);
      double lambdam = fabs( t1 - t2 * t3 - shift[k]);
      
      if(lambdap > lambdamax)
        lambdamax = lambdap;
      
      if(lambdam > lambdamax)
        lambdamax = lambdam;
    }
  
  return lambdamax;
}

void compute_conserved_quantities_from_ICs_general_relativity(int i)
{
    
  // Position of Cells
  
  double x = SphP[i].Center[0];
  double y = SphP[i].Center[1];
  double z = SphP[i].Center[2];
    
  // metric quantities in cfc - just declaration
  
  double alp, psi, bx, by, bz;
    
  get_metric_general_relativity(x,y,z,&alp,&psi,&bx,&by,&bz);  
   
  double det = pow(psi,6.);

  double ps4 = pow(psi,4.);     
    
  double v2 = P[i].Vel[0]*P[i].Vel[0] + P[i].Vel[1]*P[i].Vel[1] + P[i].Vel[2]*P[i].Vel[2];
  
  double w = 1. / sqrt( 1. - v2 * ps4 );
  
#ifdef MESHRELAX_DENSITY_IN_INPUT
  P[i].Mass = det * SphP[i].Density * w * SphP[i].Volume;
#else
  SphP[i].Density = P[i].Mass / (det * w * SphP[i].Volume);
#endif


#if ADIABATIC_GENERAL_RELATIVITY==1 || ADIABATIC_GENERAL_RELATIVITY==2
  SphP[i].Pressure = kpoly_gr * pow(SphP[i].Density, GAMMA);

  double Enthalpy = 1. + SphP[i].Utherm + SphP[i].Pressure / SphP[i].Density;
#else  
  SphP[i].Pressure = (GAMMA-1.) * SphP[i].Utherm * SphP[i].Density;

  double Enthalpy = 1. + GAMMA * SphP[i].Pressure/ SphP[i].Density / (GAMMA-1.); //+ SphP[i].Utherm + SphP[i].Pressure / SphP[i].Density;
  Enthalpy = 1. + SphP[i].Utherm + SphP[i].Pressure / SphP[i].Density;

#endif
  
  int k;
  for(k = 0; k < 3; k++)
    SphP[i].Momentum[k] = P[i].Mass * w * Enthalpy * P[i].Vel[k] * ps4;  // the last ps4 because v_j=ps4*v^j
    
    // here we need to be careful when the spatial metric is not diagonal any longer
  
  SphP[i].Energy = P[i].Mass * w * Enthalpy - det * SphP[i].Pressure * SphP[i].Volume - P[i].Mass;
  
  
}

void compute_conserved_quantities_from_ICs_general_relativity_apply_atmosphere(int i)
{
    
  // Position of Cells
  
  double x = SphP[i].Center[0];
  double y = SphP[i].Center[1];
  double z = SphP[i].Center[2];
    
  // metric quantities in cfc - just declaration
  
  double alp, psi, bx, by, bz;
    
  get_metric_general_relativity(x,y,z,&alp,&psi,&bx,&by,&bz);  
   
  double det = pow(psi,6.);

  double ps4 = pow(psi,4.);     
  
  if (SphP[i].Density < rhomax_gr * atm_fct * thr_fct)
  {
          SphP[i].Density     = rhomax_gr * atm_fct;
          P[i].Vel[0]         = 0.0;
          P[i].Vel[1]         = 0.0;
          P[i].Vel[2]         = 0.0;
          SphP[i].Pressure    = kpoly_gr * pow(SphP[i].Density, GAMMA);
          SphP[i].Utherm      = kpoly_gr * pow(SphP[i].Density, GAMMA-1.0)/(GAMMA-1.);
  }
    
  double v2 = P[i].Vel[0]*P[i].Vel[0] + P[i].Vel[1]*P[i].Vel[1] + P[i].Vel[2]*P[i].Vel[2];
  
  double w = 1. / sqrt( 1. - v2 * ps4 );
  
#ifdef MESHRELAX_DENSITY_IN_INPUT
  P[i].Mass = det * SphP[i].Density * w * SphP[i].Volume;
#else
  SphP[i].Density = P[i].Mass / (det * w * SphP[i].Volume);
#endif


#if ADIABATIC_GENERAL_RELATIVITY==1 || ADIABATIC_GENERAL_RELATIVITY==2
  SphP[i].Pressure = kpoly_gr * pow(SphP[i].Density, GAMMA);

  double Enthalpy = 1. + SphP[i].Utherm + SphP[i].Pressure / SphP[i].Density;
#else  
  SphP[i].Pressure = (GAMMA-1.) * SphP[i].Utherm * SphP[i].Density;

  double Enthalpy = 1. + GAMMA * SphP[i].Pressure/ SphP[i].Density / (GAMMA-1.); //+ SphP[i].Utherm + SphP[i].Pressure / SphP[i].Density;
  Enthalpy = 1. + SphP[i].Utherm + SphP[i].Pressure / SphP[i].Density;

#endif
  
  int k;
  for(k = 0; k < 3; k++)
    SphP[i].Momentum[k] = P[i].Mass * w * Enthalpy * P[i].Vel[k] * ps4;  // the last ps4 because v_j=ps4*v^j
    
    // here we need to be careful when the spatial metric is not diagonal any longer
  
  SphP[i].Energy = P[i].Mass * w * Enthalpy - det * SphP[i].Pressure * SphP[i].Volume - P[i].Mass;
  
  
}

void do_general_relativity_source_terms()
{
  int idx, i;
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      double dt_cell = (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval; /* half the timestep of the cell */

      double sx, sy, sz, tau;
      compute_source_terms_general_relativity( P, SphP, i, &sx, &sy, &sz, &tau);
      
      /*SphP[i].Momentum[0] += sx * dt_cell * SphP[i].Volume;
      SphP[i].Momentum[1] += sy * dt_cell * SphP[i].Volume;
      SphP[i].Momentum[2] += sz * dt_cell * SphP[i].Volume;
      
      SphP[i].Energy += tau * dt_cell * SphP[i].Volume;*/
      
      // store soure terms
      SphP[i].srcMomentum[0] = sx  * dt_cell * SphP[i].Volume;
      SphP[i].srcMomentum[1] = sy  * dt_cell * SphP[i].Volume;
      SphP[i].srcMomentum[2] = sz  * dt_cell * SphP[i].Volume;
      
      SphP[i].srcEnergy      = tau * dt_cell * SphP[i].Volume;
      
    }
  
  // maybe not even needed here because we update later
  //update_primitive_variables();
  
#if ATMOSPHERE_GENERAL_RELATIVITY == 1

 general_relativity_determine_maximum_density();  // this should maybe go to the finite volume solver
  
#endif  
  

}

void apply_general_relativity_source_terms()
{
  int idx, i;
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;
    
      SphP[i].Momentum[0] += SphP[i].srcMomentum[0];
      SphP[i].Momentum[1] += SphP[i].srcMomentum[1];
      SphP[i].Momentum[2] += SphP[i].srcMomentum[2];
      
      SphP[i].Energy += SphP[i].srcEnergy;

    }
  
  // maybe not even needed here because we update later
  //update_primitive_variables(); no, it's done outside
  
}  
  
#if METRIC_TYPE == 3

int read_fixed_numerical_1d_metric(void)
{
    FILE *file;
      
    file = fopen("fixed_1d_metric.dat", "r"); 
    if(file == NULL)
    {
      perror("error opening 1d metric table file");
      fprintf(stderr, "the filname was `%s'\n","fixed_1d_metric.dat" );
      return -1;
    }
    
    oned_metric= malloc(sizeof(struct oned_metric));
    if(oned_metric == NULL)
    {
      perror("could not allocate memory for 1D metric");
      fclose(file);
      return -1;
    }
    
    if(fscanf(file, "%d", &oned_metric->nsize) != 1)  // todo is d correct?
    {
      fprintf(stderr, "error in metric 1d file format");
      free(oned_metric);
      fclose(file);
      return -1;
    }
    
    if ( !(oned_metric->nsize <= 10000) )
    {
        printf(" the metric grid is quite large, are you sure? \n");
        printf(" also for performance reasons if mteric_type=3 %d \n",oned_metric);
        terminate("metric grid size");
    }
    
    oned_metric->r   = malloc(oned_metric->nsize * sizeof(double));
    
    oned_metric->alp = malloc(oned_metric->nsize * sizeof(double));
    oned_metric->psi = malloc(oned_metric->nsize * sizeof(double)); 
    oned_metric->bx  = malloc(oned_metric->nsize * sizeof(double));
    oned_metric->by  = malloc(oned_metric->nsize * sizeof(double));
    oned_metric->bz  = malloc(oned_metric->nsize * sizeof(double));
    
    if ( (oned_metric->r == NULL) || (oned_metric->alp == NULL) || (oned_metric->psi == NULL) || (oned_metric->bx == NULL) || 
       (oned_metric->by == NULL) || (oned_metric->bz == NULL) )
    {      
        perror("error allocating metric table");
        if (!oned_metric->r) 
            free(oned_metric->r);
        if (!oned_metric->alp) 
            free(oned_metric->alp);            
        if (!oned_metric->psi) 
            free(oned_metric->psi);
        if (!oned_metric->bx) 
            free(oned_metric->bx);
        if (!oned_metric->by) 
            free(oned_metric->by);
        if (!oned_metric->bz) 
            free(oned_metric->bz);  
      free(oned_metric);
      fclose(file);
      return -1;
    }
 
    // the numerical shift vector is actually not used because menaingless in 1d
 
    int i;
    for(i = 0; i < oned_metric->nsize; i++)
      {
          
       double  d1,d2,d3,d4,d5,d6;
          
       if ( fscanf(file, "%lf %lf %lf %lf %lf %lf", &d1,&d2,&d3,&d4,&d5,&d6) != 6)
          {
            fprintf(stderr, "error in metric file at element %d\n", i);
            free(oned_metric->r);
            free(oned_metric->alp);
            free(oned_metric->psi);
            free(oned_metric->bx);
            free(oned_metric->by);
            free(oned_metric->bz);
            free(oned_metric);
            fclose(file);
            return -1;
          }
         
       //printf(" ja %g %g %g %g %g %g \n", d1,d2,d3,d4,d5,d6);
           
       //printf(" no %lf %lf %lf %lf %lf %lf \n", d1,d2,d3,d4,d5,d6);
       
       /*if(fscanf(file, "%lf %lf %lf %lf %lf %lf", oned_metric->r + i, oned_metric->alp + i, oned_metric->psi + i, oned_metric->bx + i, oned_metric->by + i, oned_metric->bz + i ) != 6)
          {
            fprintf(stderr, "error in metric file at element %d\n", i);
            free(oned_metric->r);
            free(oned_metric->alp);
            free(oned_metric->psi);
            free(oned_metric->bx);
            free(oned_metric->by);
            free(oned_metric->bz);
            free(oned_metric);
            fclose(file);
            return -1;
          }*/
       
       oned_metric->r[i]   = d1;
       oned_metric->alp[i] = d2;      
       oned_metric->psi[i] = d3;     
       oned_metric->bx[i]  = d4;
       oned_metric->by[i]  = d5;       
       oned_metric->bz[i]  = d6;       
        
       //fscanf(file, "%lf %lf %lf %lf %lf %lf", oned_metric->r + i, oned_metric->alp + i, oned_metric->psi + i, //oned_metric->bx + i, oned_metric->by + i, oned_metric->bz + i )
     
     
       //printf(" was %g %g %g %g %g %g\n", *oned_metric->r, *oned_metric->alp, *oned_metric->psi, *oned_metric->bx, *oned_metric->by, *oned_metric->bz );
      }

    int k;
    for (k=0;k<oned_metric->nsize;k++)
    {
      //printf(" was %g %g %g\n", oned_metric->r[k], oned_metric->alp[k], oned_metric->psi[k]);
    }
//    terminate("file ok?\n");  
      
    double dum1, dum2, dum3;
    if(fscanf(file, "%lf %lf %lf %lf", &BH_Mass, &dum1, &dum2, &dum3) != 4)  // todo is d correct?
    //if(fscanf(file, "%lf %lf %lf %lf", &BH_Mass, &BH_x, &BH_y, &BH_z) != 4)  // todo is d correct?        
    {
      fprintf(stderr, "error in metric 1d source properties");
      free(oned_metric);
      fclose(file);
      return -1;
    }  
    else
    {
        printf("source parameter: %g %g %g %g \n",BH_Mass,BH_x, BH_y, BH_z);
        printf("the center will be overwritten by hand \n");
    }
      
  fclose(file);
  return 0;   
}


void fixed_numerical_1d_metric_deinit(void)  // this should be called somewhere
{
   free(oned_metric->r);
   free(oned_metric->alp);
   free(oned_metric->psi);
   free(oned_metric->bx);
   free(oned_metric->by);
   free(oned_metric->bz);
   free(oned_metric);

}

#endif

#if METRIC_TYPE == 4
    
int read_fixed_numerical_1d_metric(void)
{
    
    FILE *file;
      
    file = fopen("fixed_1d_metric_equigrid.dat", "r"); 
    if(file == NULL)
    {
      perror("error opening 1d metric table file");
      fprintf(stderr, "the filname was `%s'\n","fixed_1d_metric_equigrid.dat" );
      return -1;
    }
    
    oned_metric= malloc(sizeof(struct oned_metric));  // use mymalloc  - todo also above
    if(oned_metric == NULL)
    {
      perror("could not allocate memory for 1D metric");
      fclose(file);
      return -1;
    }
    
    if(fscanf(file, "%d %lf %lf %lf", &oned_metric->nsize, &oned_metric->r0, &oned_metric->dr, &oned_metric->rmax) != 4)  // todo is d correct?
    {
      fprintf(stderr, "error in metric 1d file format");
      free(oned_metric);
      fclose(file);
      return -1;
    }
    
    printf("metric grid parameter: %d %g %g %g\n",oned_metric->nsize, oned_metric->r0, oned_metric->dr, oned_metric->rmax);
    
    //terminate("gude");
    
    oned_metric->alp = malloc(oned_metric->nsize * sizeof(double));
    oned_metric->psi = malloc(oned_metric->nsize * sizeof(double)); 
    
    if ( (oned_metric->alp == NULL) || (oned_metric->psi == NULL) )
    {      
        perror("error allocating metric table");

        if (!oned_metric->alp) 
            free(oned_metric->alp);            
        if (!oned_metric->psi) 
            free(oned_metric->psi);
      free(oned_metric);
      fclose(file);
      return -1;
    }
 
    // the numerical shift vector is actually not used because menaingless in 1d
 
    int i;
    for(i = 0; i < oned_metric->nsize; i++)
      {
          
       double  d1,d2;
          
       if ( fscanf(file, "%lf %lf", &d1,&d2) != 2)
          {
            fprintf(stderr, "error in metric file at element %d\n", i);
            free(oned_metric->alp);
            free(oned_metric->psi);
            free(oned_metric);
            fclose(file);
            return -1;
          }
            
       oned_metric->alp[i] = d1;      
       oned_metric->psi[i] = d2;         
     
       //printf(" was %g %g \n", *oned_metric->alp, *oned_metric->psi );
       
      }

    int k;
    for (k=0;k<oned_metric->nsize;k++)
    {
      //printf(" was %g %g\n", oned_metric->alp[k], oned_metric->psi[k]);
    }
//    terminate("file ok?\n");  
      
    double dum1, dum2, dum3;
    if(fscanf(file, "%lf %lf %lf %lf", &BH_Mass, &dum1, &dum2, &dum3) != 4)  // todo is d correct?
    //if(fscanf(file, "%lf %lf %lf %lf", &BH_Mass, &BH_x, &BH_y, &BH_z) != 4)  // todo is d correct?        
    {
      fprintf(stderr, "error in metric 1d source properties");
      free(oned_metric);
      fclose(file);
      return -1;
    }  
    else
    {
        printf("source parameter: %g %g %g %g \n",BH_Mass,BH_x, BH_y, BH_z);
        printf("the center will be overwritten by hand \n");
    }
      
  fclose(file);
  return 0;   
}


void fixed_numerical_1d_metric_deinit(void)  // this should be called somewhere
{

   free(oned_metric->alp);
   free(oned_metric->psi);
   free(oned_metric);

}

#endif

#ifdef ATMOSPHERE_GENERAL_RELATIVITY
void general_relativity_determine_maximum_density()
{  // for the atmosphere treatment - todo check if it does what it should, it should be called in finite volume solver
    
  double rhomax  = 0.;
      
  int i;
  for(i = 0; i < NumGas; i++)
    {
      if(SphP[i].Density > rhomax)
        rhomax = SphP[i].Density;
    }
  
  double rhomaxAll;
    
  MPI_Allreduce(&rhomax, &rhomaxAll, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  
  rhomax_gr = rhomaxAll;
  
}
#endif

void general_relativity_statistics()
{
  double rhomax  = 0.;
  double rhomin  = 1.e10;
  double umax    = 0.0;
  double mrest   = 0.0;
  double mgrav   = 0.0;
  double eint    = 0.0;
  double vmax    = 0.0;
  double ucenter = 0.0;
  double rmin    = 1.e10;
  
  double taumax  = -100.;
  double taumin  = 100.;
  double tausum = 0.0;
  double rstmax  = 0.0;
  int    icenter;
  
  int i;
  for(i = 0; i < NumGas; i++)
    {
      if(SphP[i].Density > rhomax)
        rhomax = SphP[i].Density;
      
      if(SphP[i].Density < rhomin)
        rhomin = SphP[i].Density;      
      
      if(SphP[i].Utherm > umax)
        umax = SphP[i].Utherm;
      
      if (SphP[i].Energy/SphP[i].Volume > taumax)
          taumax = SphP[i].Energy/SphP[i].Volume;
      
      if (SphP[i].Energy/SphP[i].Volume < taumin)
          taumin = SphP[i].Energy/SphP[i].Volume;      

      if (P[i].Mass/SphP[i].Volume > rstmax)
          rstmax = P[i].Mass/SphP[i].Volume;
      
      double v2 = sqrt( P[i].Vel[0]*P[i].Vel[0] + P[i].Vel[1]*P[i].Vel[1] + P[i].Vel[2]*P[i].Vel[2] );
      
      if(v2 > vmax)
        vmax = v2;     
 
      double r =    (SphP[i].Center[0]-BH_x)*(SphP[i].Center[0]-BH_x) 
                 +  (SphP[i].Center[1]-BH_y)*(SphP[i].Center[1]-BH_y) 
                 +  (SphP[i].Center[2]-BH_z)*(SphP[i].Center[2]-BH_z);
 
      if ( r < rmin )
      {
          rmin     = r;
          ucenter  = SphP[i].Utherm;   
          icenter  = i;
      }
      
      mrest  = mrest + P[i].Mass;
      eint   = eint  + P[i].Mass * SphP[i].Utherm;
      tausum = tausum + SphP[i].Energy/SphP[i].Volume;
      
      // todo - implemt also gravitational mass
            
    }
  
  double rhomaxAll;
  double rhominAll;
  double umaxAll;
  double vmaxAll;
  double rminAll;
  double ucenterAll;
  int    icenterAll;
  double mrestAll;
  double eintAll;
  double minlapseAll;
  double taumaxAll;
  double tauminAll;
  double tausumAll;
  double rstmaxAll;
  double psimaxAll;

  
  MPI_Reduce(&rhomax, &rhomaxAll, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&umax,   &umaxAll,   1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);  
  MPI_Reduce(&vmax,   &vmaxAll,   1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&psimax, &psimaxAll, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); 
  MPI_Reduce(&taumax, &taumaxAll, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&rstmax, &rstmaxAll, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); 
  
  MPI_Reduce(&rhomin, &rhominAll,     1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD); 
  MPI_Reduce(&minlapse, &minlapseAll, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD); 
  MPI_Reduce(&taumin, &tauminAll,     1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

  MPI_Reduce(&tausum,   &tausumAll,1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&mrest,  &mrestAll,   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);  
  MPI_Reduce(&eint,   &eintAll,    1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
   
  // todo check if this works to find central values
  
  struct
  {
      double valr;
      int rank;
  } incenter[1], outcenter[1];
  
  int myrank;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank); // is this ThisTask??
  
  incenter[0].valr = rmin;
  incenter[0].rank = myrank;
  
  MPI_Allreduce(incenter,outcenter,1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
   
  int mintask = outcenter[0].rank;
  
  MPI_Bcast(&ucenter,1,MPI_DOUBLE, mintask, MPI_COMM_WORLD);
  
  if(ThisTask == 0)
    {
        
      //printf("Salve, ich bin root and got %g from %d\n",ucenter,mintask);
      //printf("umax all ist %g\n",umaxAll);  

      fprintf(FdGR, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", All.Time, rhomaxAll, rhominAll, umaxAll, vmaxAll, mrestAll, eintAll, ucenter, minlapseAll, psimaxAll, tausumAll, taumaxAll, tauminAll, rstmaxAll);
      myflush(FdGR);
    }
    
}


#endif




