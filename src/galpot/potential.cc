#include "common.h"
#include "potential.h"
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_expint.h>

// Potential class is the main class

Potential::~Potential() {}

// Constructor for subclass
AxisymmetricMultipolePotential::AxisymmetricMultipolePotential(double G_, const AxisymmetricDensity &rho_,
                                                               int nr_,
                                                               int npoly_, int ngauss_) :
  nr(nr_), npoly(npoly_), ngauss(ngauss_), G(G_), pi(M_PI), rho(&rho_), chatty(0),
  poly(ngauss,npoly), Phil(npoly,nr), dPhil(npoly,nr), ddPhil(npoly,nr)
{
    logr = new double[nr];
    r = new double[nr];
    cosu = new double[ngauss];
    for(int i=0; i<nr; i++) {
        logr[i] = log(rho->rmin) + i*log(rho->rmax/rho->rmin)/(nr-1.0);
        r[i] = exp(logr[i]);
    }
    calculate_multipoles();
}

//Destructor
AxisymmetricMultipolePotential::~AxisymmetricMultipolePotential() {
    delete[] logr;
    delete[] r;
    delete[] cosu;
}

/* Calculate multipole expansion of potential.  See Section 2.1 of thesis.
 The angular integrations in the calculation of rho_l(a) (eq. 2-2) are
 performed using Gauss-Legendre quadrature.  The variable ngauss sets
 the order of this.
 As output, this function returns the quantity inside the square brackets
 of equs 2-1 and 2-4 in the arrays Phil and dPhil respectively.
 These can then be used to interpolate what the potential and forces
 are at any point.
 
 Some points:
 Within the function below, Phil1 and Phil2 are used to store the
 results of integrating rho_l(a)*a^(l+2) and rho_l(a)/a^(l-1) between
 successive grid points.  Function intpot() does the actual integration.
 */

void AxisymmetricMultipolePotential::calculate_multipoles() {
    double wu[ngauss], pol[npoly];
    double Phil1[nr][npoly], Phil2[nr][npoly];
    
    /* Get vertices and weights for the Gauss-Legendre quadrature. */
    gauleg(0.0,1.0,cosu,wu,ngauss);
    for(int i=0;i<ngauss;i++) {
        double pol[npoly];
        wu[i]*=2.0; ///
        legend(pol,cosu[i]);
        for(int np=0;np<npoly;np++)
            poly[i][np]=pol[np]*wu[i];
    }
    
    
    /*  Now do integrations. */
    if(chatty) {
        cout << endl << "Doing multipole expansion";
        cout.flush();
    }
    for(int np=0;np<npoly;np++) {
        if(chatty) cout << "."; cout.flush();
        int l=2*np; ///
        for(int i=0;i<nr;i++) {
            double r1 = (i>0) ? r[i-1] : r[0]*0.4;
            double r2 = (i<nr-1) ? r[i+1] : r[nr-1]*1.4;
            Phil1[i][np]=integrate_rhol(l,l+2,r1,r[i]); // int r^{l+2} rho_l(r) dr between r1 and r[i]
            Phil2[i][np]=integrate_rhol(l,1-l,r[i],r2); // int r^{1-l} rho_l(r) dr between r[i] and r2
        }
        /* Finally, sum up Phil1 and Phil2 to get Phil and dPhil. */
        for(int n=0;n<nr;n++) {
            double p1 = 0.0, p2=0.0;
            double a = r[n];
            for(int i=0;i<=n;i++) p1+=Phil1[i][np];
            for(int i=nr-1;i>=n;i--) p2+=Phil2[i][np];
            Phil[np][n]=-2*pi*G*(p1/pow(a,(double) l+1)+
                                 p2*pow(a,(double) l));
            dPhil[np][n]=-2*pi*G*(-p1/pow(a,(double) l+2)*(l+1)+
                                  p2*pow(a,(double) l-1)*l);
            ddPhil[np][n]=-2*pi*G*( p1*(l+1)*(l+2)/pow(a,(double) l+3)
                                   + p2*l*(l-1)*pow(a,(double) l-2) - (2*l+1)*rhol(l,a));
        }
    }
    
    if(chatty) cout << "done" << endl;
}

/*
 The function legend() returns the values of the
 _even_ Legendre polynomials evaluated at some point.
 dlegend() does the same for the derivatives.
 */
void AxisymmetricMultipolePotential::legend(double *pol,double c) const { ///
    double c2=c*c;
    if(npoly>0) pol[0]=1.0;
    if(npoly>1) pol[1]=1.5*c2-.5;
    for(int np=2;np<npoly;np++) {
        int l=2*(np-1);
        int l2=2*l;
        pol[np]=-pol[np-2]*l*(l-1.)/((l2+1.)*(l2-1.))+
        pol[np-1]*(c2-(l2*l+l2-1.)/((l2-1.)*(l2+3.)));
        pol[np]*=(l2+1.0)*(l2+3.0)/(l+1.0)/(l+2.0);
    }
}

void AxisymmetricMultipolePotential::dlegend(double *dpol,double c) const { ///
    double *pol0 = new double[ max(2*npoly-1,1) ];
    double *dpol0 = new double[ max(2*npoly-1,1) ];
    pol0[0]=1.0; pol0[1]=c;
    dpol0[0]=dpol[0]=0.0; dpol0[1]=1.0;
    for(int l=1;l<=2*npoly-3;l++) {
        pol0[l+1]=2*c*pol0[l]-pol0[l-1]-(c*pol0[l]-pol0[l-1])/(l+1);
        dpol0[l+1]=(l+1)*pol0[l]+c*dpol0[l];
        if(l%2==1) dpol[l/2+1]=dpol0[l+1];
    }
    delete[] dpol0;
    delete[] pol0;
}


/*
 * Calculate rho_l(a) using Gauss-Legendre quadrature
 */
double AxisymmetricMultipolePotential::rhol(int l, double a) const {
    double sum=0.0;
    int np = l/2; ///
    for(int i=0;i<ngauss;i++) {
        double R = a*sqrt(1.0-cosu[i]*cosu[i]);
        double z = a*cosu[i];
        sum += (*rho)(R,z)*poly[i][np];
    }
    return sum;
}

// Do Int_{a1}^{a2} a^m rho_l(a) da
//  = Int_{a1}^{a2} a^{m+1} rho_l(a) d(log a)
double AxisymmetricMultipolePotential::integrate_rhol(int l, int m, double a1, double a2) {
    const int nstep = 200;
    double sum = 0.0;
    double loga1 = log(a1);
    double dloga = log(a2/a1)/nstep;
    
    double one = 1.0;
    
    sum = rhol(l,a1)*pow(a1, m+one);
    for(int i=1;i<nstep;i++) {
        double loga = loga1 + i*dloga;
        double a = exp(loga);
        double fac=(i%2==0) ? 2.0 : 4.0;
        sum += fac*rhol(l,a)*pow(a, m+one);
    }
    sum += rhol(l,a2)*pow(a2, m+one);
    
    return sum*dloga/3.0;
}

double AxisymmetricMultipolePotential::dopot(double logr0,double *pol, const Stupid2DArray &Phil) const {
    double sum=0.0;
    for(int np=npoly-1;np>=0;np--)
        sum += pol[np]*lininterp(nr,logr,Phil[np],logr0);
    return sum;
}


// Finally, the useful public functions
double AxisymmetricMultipolePotential::Phi(double R, double z, double t) const {
    double pol[npoly];
    double a = hypot(R,z);
    double c = z/a;
    double loga0 = log(a);
    legend(pol,c);
    return dopot(loga0,pol,Phil);
}

Vec2 AxisymmetricMultipolePotential::dPhidRz(double R, double z) const {
    double pol[npoly];
    double a = hypot(R,z);
    double c = z/a;
    double loga0 = log(a);
    legend(pol,c);
    double sum1 = dopot(loga0,pol,dPhil);
    dlegend(pol,c);
    double sum2 = dopot(loga0,pol,Phil);
    return Vec2( R/a*sum1 - R*z/(a*a*a)*sum2,
                z/a*sum1 + R*R/(a*a*a)*sum2 );
}

double AxisymmetricMultipolePotential::lininterp(int n,const double *xa, const double *ya,double x) const {
    int klo=0,khi=n-1,k;
    double sgn=(xa[1]>=xa[0]) ? 1 : -1;
    while(khi-klo>1) {
        k=(khi+klo)/2;
        if(sgn*xa[k]>sgn*x) khi=k;
        else klo=k;
    }
    k=klo;
    return (ya[k+1]-ya[k])/(xa[k+1]-xa[k])*(x-xa[k])+ya[k];
}

// Return coeffs for Gauss-Legendre quadrature.
// Stolen from NR.  Modified to output results in
// x[0..n-1], w[0..n-1] instead of x[1..n],w[1..n]
#define EPS 3.0e-11
void AxisymmetricMultipolePotential::gauleg(double x1,double x2,double *x,double *w,int n) {
    --x;
    --w;
    int m,j,i;
    double z1,z,xm,xl,pp,p3,p2,p1;
    
    m=(n+1)/2;
    xm=0.5*(x2+x1);
    xl=0.5*(x2-x1);
    for (i=1;i<=m;i++)  {
        z=cos(pi*(i-0.25)/(n+0.5));
        do {
            p1=1.0;
            p2=0.0;
            for (j=1;j<=n;j++) {
                p3=p2;
                p2=p1;
                p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
            }
            pp=n*(z*p1-p2)/(z*z-1.0);
            z1=z;
            z=z1-p1/pp;
        } while (fabs(z-z1) > EPS);
        x[i]=xm-xl*z;
        x[n+1-i]=xm+xl*z;
        w[i]=2.0*xl/((1.0-z*z)*pp*pp);
        w[n+1-i]=w[i];
    }
}
//#undef EPS

// The functions that I actually use

double AxisymmetricPotential::Phi(const Vec3 &x, double t) {
  return Phi( hypot(x[0],x[1]), x[2] );
}

Vec3 AxisymmetricPotential::dPhidx(const Vec3 &x, double t) {
  double R = hypot(x[0],x[1]), z = x[2];
  Vec2 derivs = dPhidRz(R,z);
  if(R>0.0) return Vec3(derivs[0]*x[0]/R, derivs[0]*x[1]/R, derivs[1]);
  else return Vec3(0.0, 0.0, derivs[1]);
}

// Define a subclass for a particular kind of axysimmetric potential
// The BGSBU potential (see Binney++ 1991)

double BGSBUBarPotential::rho(const Vec3 &x) const {
    double R = sqrt(x[1]*x[1]+x[2]*x[2]);
    double z = x[0];
    return rho_priv(R,z);
}

double BGSBUBarPotential::Phi(const Vec3 &x, double t) {
    return this->AxisymmetricMultipolePotential::Phi(Vec3(x[1],x[2],x[0]));
}

Vec3 BGSBUBarPotential::dPhidx(const Vec3 &x, double t) {
    Vec3 ans = this->AxisymmetricMultipolePotential::dPhidx(Vec3(x[1],x[2],x[0]));
    return Vec3(ans[2],ans[0],ans[1]);
}

RazorDiskPotential::RazorDiskPotential(double Gin, double Sigma0in, double Rdin) :
G(Gin), Sigma0(Sigma0in), Rd(Rdin) {}

double RazorDiskPotential::Phi(const Vec3 &x, double t) {
    double R = sqrt(x[0]*x[0] + x[1]*x[1]);
    double y = R / (2.0 * Rd);
    double ans = 0.0;
    if(R>0 && y<12){ans = - pi*G*Sigma0*R*(gsl_sf_bessel_I0(y)*gsl_sf_bessel_K1(y) - gsl_sf_bessel_I1(y)*gsl_sf_bessel_K0(y));}
    return ans;
};

Vec3 RazorDiskPotential::dPhidx(const Vec3 &x, double t) {
    double R = sqrt(x[0]*x[0] + x[1]*x[1]);
    double y = R / (2.0 * Rd);
    double dPhidR = 0.0;
    if(R>0 && y<12){dPhidR = pi*G*Sigma0*(R/Rd)*(gsl_sf_bessel_I0(y)*gsl_sf_bessel_K0(y) - gsl_sf_bessel_I1(y)*gsl_sf_bessel_K1(y));}
    return Vec3(dPhidR * (x[0]/R),
                dPhidR * (x[1]/R),
                0.0);
};

StupidBarPotential::StupidBarPotential(double Ain, double R0in) :
A(Ain), R0(R0in) {}

double StupidBarPotential::Phi(const Vec3 &x, double t) {
  double R = sqrt(x[0]*x[0]+x[1]*x[1]);
  double theta = atan2(x[1], x[0]);
  double cos2theta = cos(2*theta);
  return -A * (R/R0) * (R/R0) * exp(-(R/R0)) * cos2theta;
}

Vec3 StupidBarPotential::dPhidx(const Vec3 &x, double t) {
  Vec3 xplusdx = x;
  Vec3 xplusdy = x;
  double dx = 1e-5;
  xplusdx[0] += dx;
  xplusdy[1] += dx;
  double dphidx = (Phi(xplusdx) - Phi(x)) / dx;
  double dphidy = (Phi(xplusdy) - Phi(x)) / dx;
  return Vec3(dphidx,
              dphidy,
              0.0);
};

LogPotential::LogPotential(double v0in, double Rcin, double qin) :
v0(v0in), Rc(Rcin), q(qin) {}

double LogPotential::Phi(const Vec3 &x, double t) {
    double r2 = x[0]*x[0]+x[1]*x[1];
    double Rc2 = Rc * Rc, q2 = q * q;
    return 0.5 * v0*v0 * log(Rc2 + x[0]*x[0] + x[1]*x[1]/q2);
}

Vec3 LogPotential::dPhidx(const Vec3 &x, double t) {
    double rq2=Rc*Rc + x[0]*x[0] + x[1]*x[1]/(q*q);
    return Vec3(v0*v0 * x[0]/rq2,
                v0*v0 * x[1]/(q*q*rq2),
                0.0);
};

/*Sormani, Binney, Magorrian 2015c quadrupole*/

double PhiTilde(double x) {
  double x2 = x*x, x3 = x2*x, x4 = x3*x, x5 = x4*x;
  double ExpIntegralEi = gsl_sf_expint_Ei(-x);
  double F1 = (24. + 24.*x + 12.*x2 + 4.*x3 + x4);
  return (exp(-x)*F1 - 24. + x5*ExpIntegralEi) / (5.*x3);
}

double dPhiTildedx(double x) {
  double x2 = x*x, x3 = x2*x, x4 = x3*x, x5 = x4*x;
  double ExpIntegralEi = gsl_sf_expint_Ei(-x);
  double F1 = (24. + 24.*x + 12.*x2 + 4.*x3 + x4);
  return -(exp(-x)*3.*F1 - 72. - 2.*x5*ExpIntegralEi) / (5.*x4);
}

SBM2015Quadrupole::SBM2015Quadrupole(double Ain, double v0in, double Rqin) :
A{Ain}, v0{v0in}, Rq{Rqin} {}

double SBM2015Quadrupole::Phi(const Vec3 &x, double t) {
    double R = sqrt(x[0]*x[0] + x[1]*x[1]);
    double phi = atan2(x[1], x[0]);
    if (R == 0.) return 0.;
    double K = A/4. * v0*v0;
    double Phi2 = K * PhiTilde(2. * R/Rq);
    return  Phi2 * cos(2.*phi);
}

Vec3 SBM2015Quadrupole::dPhidx(const Vec3 &x, double t) {
  
    double eps = 0.001;

    double R = sqrt(x[0]*x[0] + x[1]*x[1]);
    double phi = atan2(x[1], x[0]);

    if (R < eps) return Vec3(0., 0., 0.);

    double K = A/4. * v0*v0;
    double Phi2 = K * PhiTilde(2. * R/Rq);
    double dPhi2dR = K * dPhiTildedx(2. * R/Rq) * (2./Rq);

    double dPhidR = dPhi2dR*cos(2.*phi);
    double dPhidphi = Phi2*2.*sin(2.*phi);

    double dPhidx_x = dPhidR*cos(phi) + dPhidphi*sin(phi)/R;
    double dPhidx_y = dPhidR*sin(phi) - dPhidphi*cos(phi)/R;    
    return Vec3(dPhidx_x,
                dPhidx_y,
                0.0);
};

/*Potential generated by a power law density distribution*/
PowerLawPotential::PowerLawPotential(double v0in, double alphain) :
v0{v0in}, alpha{alphain} {}

double PowerLawPotential::Phi(const Vec3 &x, double t) {
  return  0.;
}

Vec3 PowerLawPotential::dPhidx(const Vec3 &x, double t) {

  double R = sqrt(x[0]*x[0] + x[1]*x[1]);
  double phi = atan2(x[1], x[0]);

  double dPhidR = v0*v0 * pow(R,1.-alpha);

  double dPhidx_x = dPhidR*cos(phi);
  double dPhidx_y = dPhidR*sin(phi);

  return Vec3(dPhidx_x,
              dPhidx_y,
              0.0);
};


/*Exponential disk*/
double ExpDiskPotential::rho(const Vec3 &x) const {
    double R = sqrt(x[0]*x[0]+x[1]*x[1]);
    double z = x[2];
    return rho_priv(R,z);
}

double ExpDiskPotential::Phi(const Vec3 &x, double t) {
    return this->AxisymmetricMultipolePotential::Phi(Vec3(x[0],x[1],x[2]));
}

Vec3 ExpDiskPotential::dPhidx(const Vec3 &x, double t) {
    Vec3 ans = this->AxisymmetricMultipolePotential::dPhidx(Vec3(x[0],x[1],x[2]));
    return Vec3(ans[0],ans[1],ans[2]);
}

/*Double power law profile*/
double DPLPotential::rho(const Vec3 &x) const {
    double R = sqrt(x[0]*x[0]+x[1]*x[1]);
    double z = x[2];
    return rho_priv(R,z);
}

double DPLPotential::Phi(const Vec3 &x, double t) {
    return this->AxisymmetricMultipolePotential::Phi(Vec3(x[0],x[1],x[2]));
}

Vec3 DPLPotential::dPhidx(const Vec3 &x, double t) {
    Vec3 ans = this->AxisymmetricMultipolePotential::dPhidx(Vec3(x[0],x[1],x[2]));
    return Vec3(ans[0],ans[1],ans[2]);
}

/*McMillan Bulge */
double McMillanBulgePotential::rho(const Vec3 &x) const {
    double R = sqrt(x[0]*x[0]+x[1]*x[1]);
    double z = x[2];
    return rho_priv(R,z);
}

double McMillanBulgePotential::Phi(const Vec3 &x, double t) {
    return this->AxisymmetricMultipolePotential::Phi(Vec3(x[0],x[1],x[2]));
}

Vec3 McMillanBulgePotential::dPhidx(const Vec3 &x, double t) {
    Vec3 ans = this->AxisymmetricMultipolePotential::dPhidx(Vec3(x[0],x[1],x[2]));
    return Vec3(ans[0],ans[1],ans[2]);
}

double ExpBarPotential::rho(const Vec3 &x) const {
    double R = sqrt(x[1]*x[1]+x[2]*x[2]);
    double z = x[0];
    return rho_priv(R,z);
}

double ExpBarPotential::Phi(const Vec3 &x, double t) {
    return this->AxisymmetricMultipolePotential::Phi(Vec3(x[1],x[2],x[0]));
}

Vec3 ExpBarPotential::dPhidx(const Vec3 &x, double t) {
    Vec3 ans = this->AxisymmetricMultipolePotential::dPhidx(Vec3(x[1],x[2],x[0]));
    return Vec3(ans[2],ans[0],ans[1]);
}

double ModifiedMcMillanBulgePotential::rho(const Vec3 &x) const {
  double R = sqrt(x[0]*x[0]+x[1]*x[1]);
  double z = x[2];
  return rho_priv(R,z);
}

double ModifiedMcMillanBulgePotential::Phi(const Vec3 &x, double t) {
  return this->AxisymmetricMultipolePotential::Phi(Vec3(x[0],x[1],x[2]));
}

Vec3 ModifiedMcMillanBulgePotential::dPhidx(const Vec3 &x, double t) {
  Vec3 ans = this->AxisymmetricMultipolePotential::dPhidx(Vec3(x[0],x[1],x[2]));
  return Vec3(ans[0],ans[1],ans[2]);
}

/*Miyamoto & Nagai potential*/
MiyamotoNagaiPotential::MiyamotoNagaiPotential(double Gin, double ain, double bin, double Min) :
G{Gin}, a{ain}, b{bin}, M{Min} {}

double MiyamotoNagaiPotential::Phi(const Vec3 &x, double t) {
  double R = sqrt(x[0]*x[0] + x[1]*x[1]);
  return  -(G * M) / sqrt(R*R + pow((a + sqrt(x[2]*x[2] + b*b)),2) );
}

Vec3 MiyamotoNagaiPotential::dPhidx(const Vec3 &x, double t) {
  
  double factor1 = a + sqrt(x[2]*x[2] + b*b); 
  double factor2 = x[0]*x[0] + x[1]*x[1] + factor1*factor1;

  double dPhidx_x = G*M * pow(factor2,-3./2.) * x[0];
  double dPhidx_y = G*M * pow(factor2,-3./2.) * x[1];
  double dPhidx_z;
  if (b == 0.)
    dPhidx_z = G*M * pow(factor2,-3./2.) * factor1 * (x[2] > 0) ? 1 : ((x[2] < 0) ? -1 : 0);
  else
    dPhidx_z = G*M * pow(factor2,-3./2.) * factor1 * x[2] / sqrt(x[2]*x[2] + b*b);

  return Vec3(dPhidx_x,
              dPhidx_y,
              dPhidx_z);
}

/*Spiral arm potential perturbation given by Cox & Gomez (2002)*/
SpiralArmPotential::SpiralArmPotential(double Sr0in, 
                                       double St0in, double SNarmsin, double Sphirin, double Salphain, 
                                       double SHzin, double Sp0in, double Srsin, double Gin) :
Sr0{Sr0in}, St0{St0in}, SNarms{SNarmsin}, 
Sphir{Sphirin}, Salpha{Salphain}, SHz{SHzin}, Sp0{Sp0in}, Srs{Srsin}, G{Gin} {}

double SpiralArmPotential::Phi(const Vec3 &x, double t) {

  // Convert to polar co-ords
  double R = sqrt(x[0]*x[0] + x[1]*x[1]);
  double phi = atan2(x[1], x[0]);

  double gamma, sum;
  double Kn, Bn, Dn;
  double np1float;
  int n;

  double Cz[3];

  Cz[0] = 8. / (3. * M_PI);
  Cz[1] = 0.5;
  Cz[2] = 8. / (15. * M_PI);
  
  // Calculate the potential 
  gamma = SNarms * (phi + Sphir * (St0 + t) - log(R / Sr0) / tan(Salpha));

  sum = 0.;
  for(n = 0; n < 3; n++)
    {
      np1float = 1. + (float) n;
      Kn = np1float * SNarms / (R * sin(Salpha));
      Bn = Kn * SHz * (1. + 0.4 * Kn * SHz);
      Dn = (1. + Kn * SHz + 0.3 * pow(Kn * SHz, 2)) / (1. + 0.3 * Kn * SHz);

      sum += (Cz[n] / (Dn * Kn)) * cos(np1float * gamma) * pow(1. / cosh((Kn * x[2]) / Bn), Bn);
    }

  return -4. * M_PI * G * SHz * Sp0 * exp(-(R - Sr0) / Srs) * sum;
}

Vec3 SpiralArmPotential::dPhidx(const Vec3 &x, double t) {

  Vec3 xplusdx = x;
  Vec3 xplusdy = x;
  Vec3 xplusdz = x;

  double dx = 1e-3;

  xplusdx[0] += dx;
  xplusdy[1] += dx;
  xplusdz[2] += dx;

  double phix = Phi(x, t);

  double dPhidx_x = -(Phi(xplusdx, t) - phix) / dx;
  double dPhidx_y = -(Phi(xplusdy, t) - phix) / dx;
  double dPhidx_z = -(Phi(xplusdz, t) - phix) / dx;

 return Vec3(dPhidx_x,
             dPhidx_y,
             dPhidx_z);
}

/*NFW potential moving along a straight line at given velocity*/
double MovingNFWPotential::rho(const Vec3 &x, double t) const {
    Vec3 center = Vec3(x0[0] + v[0]*t,
                       x0[1] + v[1]*t,
                       x0[2] + v[2]*t);
    
    Vec3 x1 = Vec3(x[0] - center[0],
                   x[1] - center[1],
                   x[2] - center[2]);

    double R = sqrt(x1[0]*x1[0]+x1[1]*x1[1]);
    double z = x1[2];
    return rho_priv(R,z);
}

double MovingNFWPotential::Phi(const Vec3 &x, double t) {
    Vec3 center = Vec3(x0[0] + v[0]*t,
                       x0[1] + v[1]*t,
                       x0[2] + v[2]*t);
    
    Vec3 x1 = Vec3(x[0] - center[0],
                   x[1] - center[1],
                   x[2] - center[2]);

    return this->AxisymmetricMultipolePotential::Phi(Vec3(x1[0],x1[1],x1[2]));
}

Vec3 MovingNFWPotential::dPhidx(const Vec3 &x, double t) {
    Vec3 center = Vec3(x0[0] + v[0]*t,
                       x0[1] + v[1]*t,
                       x0[2] + v[2]*t);
 
    Vec3 x1 = Vec3(x[0] - center[0],
                   x[1] - center[1],
                   x[2] - center[2]);
    
    Vec3 ans = this->AxisymmetricMultipolePotential::dPhidx(Vec3(x1[0],x1[1],x1[2]));
    return Vec3(ans[0],ans[1],ans[2]);
}

