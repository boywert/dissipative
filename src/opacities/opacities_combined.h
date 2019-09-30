#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>

/* data structure for fa opacities */
typedef struct {
  int nr;
  int nt;
  double *logR, *logT;
  double *logkappa;
  gsl_interp2d * inter;
  gsl_interp_accel *xacc;
  gsl_interp_accel *yacc;
} opacitydata;

/* functions for combined opacities */
float opacities(float rho, float t, opacitydata *data);
void initialize_opacities(opacitydata *data, char *dir);
void finalize_opacities(opacitydata *data);

/* functions for OPAL opacities */
void init_opal(char * path);
void opacity_opal(float * xh, float * rho, float * t, float * kappa);
