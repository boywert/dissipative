class Potential {
 public:
  virtual ~Potential()=0;
  virtual double Phi(const Vec3 &x, double t=0.0)=0;
  virtual Vec3 dPhidx(const Vec3 &x, double t=0.0)=0;
};

class Stupid2DArray {
 public:
  Stupid2DArray(int n1_, int n2_) : n1(n1_), n2(n2_) {
      arr = new double[n1*n2];
  }
  ~Stupid2DArray() {
      delete[] arr;
  }
  inline double* operator[](int i) { return (arr + (n2 * i)); }
  inline double const* const operator[](int i) const { return (arr + (n2 * i)); }
 private:
  const int n1, n2;
  double *arr;
};


class AxisymmetricPotential : public Potential {
 public:
  virtual double Phi(double R, double z, double t=0.0) const=0;
  virtual double Phi(const Vec3 &x, double t=0.0);
  virtual Vec3 dPhidx(const Vec3 &x, double t=0.0);
  virtual Vec2 dPhidRz(double R, double z) const=0;
};


class AxisymmetricDensity {
 public:
  const double rmin, rmax;
  AxisymmetricDensity(double r1, double r2) : rmin(r1), rmax(r2) {}
  virtual double operator()(double, double) const = 0;
  virtual ~AxisymmetricDensity() = 0;
};
inline AxisymmetricDensity::~AxisymmetricDensity() {}

/*
*
* The main event.  Class to represent (axisymmetric) potential calculated
* by multipole expansion.
* This version assumes reflection symmetry in z=0.
*
*/
class AxisymmetricMultipolePotential : public AxisymmetricPotential {

public:
  const double G;
  const double pi;
  AxisymmetricMultipolePotential(double G, const AxisymmetricDensity &, int nr,
                                 int npoly=6, int ngauss=8);
  ~AxisymmetricMultipolePotential();
  virtual double rhol(int, double) const;
  virtual double Phi(double R, double z, double t=0.0) const;
  virtual Vec2 dPhidRz(double R, double z) const;
  virtual double Phi(const Vec3 &x, double t=0.0) {
    return AxisymmetricPotential::Phi(x);
  }
  virtual Vec3 dPhidx(const Vec3 &x, double t=0.0) {
    return AxisymmetricPotential::dPhidx(x);
  }
    
private:
  int chatty;
  const int nr, npoly, ngauss;
  double *logr, *r, *cosu;
  const AxisymmetricDensity *rho;
  Stupid2DArray Phil, dPhil, ddPhil, poly;
  void calculate_multipoles();
  void init();
  void legend(double *pol, double c) const;
  void dlegend(double *dpol, double c) const;
  double integrate_rhol(int, int, double, double);
  double dopot(double ga,double *pol, const Stupid2DArray &Phil) const;
  void gauleg(double x1, double x2, double *x, double *w, int n);
  double lininterp(int, const double *, const double *, double) const;
};


/*
 *
 * Class used to represent BGBSU91 (Binney++91) bar potential.
 * The density is given by equation (2) of that paper.
 * The potential is calculated using an axisymmetric multipole expansion,
 * taking the x-axis to be the symmetry axis.
 *
 */

// rminmax(rmin,a0_*1e-3), rminmax(rmax,a0_*1e2):
// these two contains the inner and outer limits of the potential calculation.
// nr=100 is the radial number of points where the potential is calculated
// npoly is the parameter that controls the order of the multipole expansion.
// The maximum value of l included in the expansion is 2*(npoly-1) [I think]

class BGSBUBarDensity : public AxisymmetricDensity {
public:
  const double alf, bet, a0, qm;
  const double rho0;
  BGSBUBarDensity(double rho0_=0.069,
                  double alf_=1.75, double bet_=3.5,
                  double a0_=1.2, double qm_=0.75, double rmin=0, double rmax=0)
  : AxisymmetricDensity(rminmax(rmin,a0_*1e-3), rminmax(rmax,a0_*1e2)),
  alf(alf_), bet(bet_), a0(a0_), qm(qm_), rho0(rho0_) {} 
  double operator()(double R, double z) const {
      double eps=0.01;
      double a = sqrt(eps*eps + z*z + (R*R)/(qm*qm));
      if(a<a0) return rho0*pow(a/a0,-alf);
      return rho0*pow(a/a0,-bet);
  }
private:
  double rminmax(double choice1, double choice2) {
      return choice1!=0 ? choice1 : choice2;
  }
};

class BGSBUBarPotential : public AxisymmetricMultipolePotential {
public:
    BGSBUBarPotential(double G, const BGSBUBarDensity &rho_in,
                      int nr=100, int npoly=6, int ngauss=8)
    : AxisymmetricMultipolePotential(G, rho_in, nr, npoly, ngauss),
    rho_priv(rho_in) {}
    double rho(const Vec3 &x) const;
    double Phi(const Vec3 &x, double t=0.0);
    Vec3 dPhidx(const Vec3 &x, double t=0.0);
private:
    BGSBUBarDensity rho_priv;
};

/*3D bar potential (Sormani, Binney & Magorrian 2015)*/
class ExpBarDensity : public AxisymmetricDensity {
public:
  const double rho0, a0, q;
  ExpBarDensity(double rho0_=0.2,
                 double a0_=1.5,
                 double q_=0.5,
                 double rmin=0, 
                 double rmax=0)
  : AxisymmetricDensity(rminmax(rmin,a0_*1e-3), rminmax(rmax,a0_*1e2)),
  rho0(rho0_), a0(a0_), q(q_) {} 
  double operator()(double R, double z) const {
  double a = sqrt(z*z + (R*R)/(q*q));
  return rho0*exp(-a/a0);
  }
private:
  double rminmax(double choice1, double choice2) {
    return choice1!=0 ? choice1 : choice2;
  }
};

class ExpBarPotential : public AxisymmetricMultipolePotential {
public:
  ExpBarPotential(double G, const ExpBarDensity &rho_in,
                   int nr=100, int npoly=6, int ngauss=8)
  : AxisymmetricMultipolePotential(G, rho_in, nr, npoly, ngauss),
  rho_priv(rho_in) {}
  double rho(const Vec3 &x) const;
  double Phi(const Vec3 &x, double t=0.0);
  Vec3 dPhidx(const Vec3 &x, double t=0.0);
private:
  ExpBarDensity rho_priv;
};

/*exponential disk potential rho(R,Z) = (Sigma0/(2*z_d)) * exp(-|z|/z_d - R/R_d)  */
class ExpDiskDensity : public AxisymmetricDensity {
public:
  const double Sigma0, zd, Rd;
  ExpDiskDensity(double Sigma0_=1.83e8 , double zd_=9., double Rd_=36., double rmin=0, double rmax=0)
  : AxisymmetricDensity(rminmax(rmin,Rd_*1e-3), rminmax(rmax,Rd_*1e2)),
  Sigma0(Sigma0_), Rd(Rd_), zd(zd_) {}
  double operator()(double R, double z) const {
      return (Sigma0/(2.*zd)) * exp(-abs(z)/zd - R/Rd);
  }
private:
  double rminmax(double choice1, double choice2) {
      return choice1!=0 ? choice1 : choice2;
  }
};

class ExpDiskPotential : public AxisymmetricMultipolePotential {
public:
    ExpDiskPotential(double G, const ExpDiskDensity &rho_in,
                      int nr=100, int npoly=6, int ngauss=8)
    : AxisymmetricMultipolePotential(G, rho_in, nr, npoly, ngauss),
    rho_priv(rho_in) {}
    double rho(const Vec3 &x) const;
    double Phi(const Vec3 &x, double t=0.0);
    Vec3 dPhidx(const Vec3 &x, double t=0.0);
private:
    ExpDiskDensity rho_priv;
};

/*Double power law profile (see eq. 2.64 of BT08)  */
class DPLDensity : public AxisymmetricDensity {
public:
  const double rho0, a, alpha, beta; 
  DPLDensity(double rho0_=8.54e3, double a_=196., double alpha_=1., double beta_=3., double rmin=0, double rmax=0)
  : AxisymmetricDensity(rminmax(rmin,a_*1e-3), rminmax(rmax,a_*1e2)),
  rho0(rho0_), a(a_), alpha(alpha_), beta(beta_) {}
  double operator()(double R, double z) const {
      double r = sqrt(R*R + z*z);
      return rho0 / ( pow(r/a,alpha) * pow((1.+r/a),(beta-alpha)) );
  }
private:
  double rminmax(double choice1, double choice2) {
      return choice1!=0 ? choice1 : choice2;
  }
};

/*McMillan bulge profile*/
class McMillanBulgeDensity : public AxisymmetricDensity {
public:
  const double alf, a0, acut, qm;
  const double rho0;
  McMillanBulgeDensity(double rho0_=9.84e7,
                       double alf_=1.8, double acut_=21.,
                       double a0_=0.75, double qm_=0.5, double rmin=0., double rmax=0.)
  : AxisymmetricDensity(rminmax(rmin,acut_*1e-3), rminmax(rmax,acut_*1e2)),
  alf(alf_), acut(acut_), a0(a0_), qm(qm_), rho0(rho0_) {}
  double operator()(double R, double z) const {
      double eps=0.01;
      double a = sqrt(eps*eps + R*R + (z*z)/(qm*qm));
      return rho0*exp(-pow(a/acut,2.))*pow(1.+a/a0,-alf);
  }
private:
  double rminmax(double choice1, double choice2) {
      return choice1!=0 ? choice1 : choice2;
  }
};

class McMillanBulgePotential : public AxisymmetricMultipolePotential {
public:
    McMillanBulgePotential(double G, const McMillanBulgeDensity &rho_in,
                      int nr=100, int npoly=6, int ngauss=8)
    : AxisymmetricMultipolePotential(G, rho_in, nr, npoly, ngauss),
    rho_priv(rho_in) {}
    double rho(const Vec3 &x) const;
    double Phi(const Vec3 &x, double t=0.0);
    Vec3 dPhidx(const Vec3 &x, double t=0.0);
private:
    McMillanBulgeDensity rho_priv;
};

/*Modified McMillan bulge profile*/
class ModifiedMcMillanBulgeDensity : public AxisymmetricDensity {
public:
  const double alf, a0, acut, qm;
  const double rho0;
  ModifiedMcMillanBulgeDensity(double rho0_=9.84e7,
                       double alf_=1.8, double a0_=1.0, double acut_=21., 
                       double qm_=0.5, double rmin=0., double rmax=0.)
  : AxisymmetricDensity(rminmax(rmin,acut_*1e-3), rminmax(rmax,acut_*1e2)),
  alf(alf_), a0(a0_), acut(acut_), qm(qm_), rho0(rho0_) {}
  double operator()(double R, double z) const {
    double eps=0.001;
    double a = sqrt(eps*eps + R*R + (z*z)/(qm*qm));
    return rho0*exp(-pow(a/acut,2.))*pow(a/a0,-alf);
  }
private:
  double rminmax(double choice1, double choice2) {
    return choice1!=0 ? choice1 : choice2;
  }
};

class ModifiedMcMillanBulgePotential : public AxisymmetricMultipolePotential {
public:
  ModifiedMcMillanBulgePotential(double G, const ModifiedMcMillanBulgeDensity &rho_in,
                         int nr=100, int npoly=6, int ngauss=8)
  : AxisymmetricMultipolePotential(G, rho_in, nr, npoly, ngauss),
  rho_priv(rho_in) {}
  double rho(const Vec3 &x) const;
  double Phi(const Vec3 &x, double t=0.0);
  Vec3 dPhidx(const Vec3 &x, double t=0.0);
private:
  ModifiedMcMillanBulgeDensity rho_priv;
};

class DPLPotential : public AxisymmetricMultipolePotential {
public:
    DPLPotential(double G, const DPLDensity &rho_in,
                      int nr=100, int npoly=6, int ngauss=8)
    : AxisymmetricMultipolePotential(G, rho_in, nr, npoly, ngauss),
    rho_priv(rho_in) {}
    double rho(const Vec3 &x) const;
    double Phi(const Vec3 &x, double t=0.0);
    Vec3 dPhidx(const Vec3 &x, double t=0.0);
private:
    DPLDensity rho_priv;
};


//*** class to calculate exponential disk potential. See BT 2.164 and 2.165
class RazorDiskPotential : public Potential {
public:
    RazorDiskPotential(double Gin, double Sigma0in, double Rdin);
    double Phi(const Vec3 &x, double t=0.0);
    Vec3 dPhidx(const Vec3 &x,double t=0.0);
    const double Sigma0, Rd, G;
};

class StupidBarPotential : public Potential {
public:
  StupidBarPotential(double Ain, double R0in);
  double Phi(const Vec3 &x, double t=0.0);
  Vec3 dPhidx(const Vec3 &x,double t=0.0);
  const double A, R0;
};

class LogPotential : public Potential {
public:
    LogPotential(double v0in, double Rcin, double qin);
    double Phi(const Vec3 &x, double t=0.0);
    Vec3 dPhidx(const Vec3 &x,double t=0.0);
    const double v0, Rc, q;
};

/*Only quadrupole potential i.e. the quadrupole is generated by
 * a density rho_0 exp(-2R/Rq) sin^2(theta) cos(2 phi) */
class SBM2015Quadrupole : public Potential {
public:
    SBM2015Quadrupole(double Ain, double v0in, double Rqin);
    double Phi(const Vec3 &x, double t=0.0);
    Vec3 dPhidx(const Vec3 &x, double t=0.0);
    const double v0, Rq, A;
};

/*Potential generated by a power law density distribution*/
class PowerLawPotential : public Potential {
public:
    PowerLawPotential(double v0in, double alphain);
    double Phi(const Vec3 &x, double t=0.0);
    Vec3 dPhidx(const Vec3 &x, double t=0.0);
    const double v0, alpha;
};

/*Miyamoto & Nagai potential*/
class MiyamotoNagaiPotential : public Potential {
public:
    MiyamotoNagaiPotential(double Gin, double ain, double bin, double Min);
    double Phi(const Vec3 &x, double t=0.0);
    Vec3 dPhidx(const Vec3 &x, double t=0.0);
    const double G, a, b, M;
};

/*Spiral arm potential perturbation given by Cox & Gomez (2002)*/
class SpiralArmPotential : public Potential {
public:
    SpiralArmPotential(double Sr0in,
                       double St0in, double SNarmsin, double Sphirin, double Salphain,
                       double SHzin, double Sp0in, double Srsin, double Gin);
    double Phi(const Vec3 &x, double t);
    Vec3 dPhidx(const Vec3 &x, double t);
    const double Sr0, St0, SNarms, Sphir, Salpha, SHz, Sp0, Srs, G;
};

/*NFW potential moving along a straight line at given velocity*/
class MovingNFWPotential : public AxisymmetricMultipolePotential {
public:
    MovingNFWPotential(double G, const DPLDensity &rho_in,
                      int nr, int npoly, int ngauss,
                      const Vec3 x0_, const Vec3 v_)
    : AxisymmetricMultipolePotential(G, rho_in, nr, npoly, ngauss),
    rho_priv(rho_in), x0(x0_), v(v_) {}
    double rho(const Vec3 &x, double t) const;
    double Phi(const Vec3 &x, double t);
    Vec3 dPhidx(const Vec3 &x, double t);
    const Vec3 x0, v;
private:
    DPLDensity rho_priv;
};
