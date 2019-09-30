#include "opacities_combined.h"

#include <stdio.h>
#include <string.h>


void initialize_fa_opacity(opacitydata *data, char *dir)
{
  int i, j, status;
  float tmp;
  FILE* f;
  //char * fname;
  char fname[200];

  data->nr = 19;
  data->nt = 85;

  data->logR = (double*)malloc(data->nr * sizeof(double));
  data->logT = (double*)malloc(data->nt * sizeof(double));
  data->logkappa = (double*)malloc(data->nt * data->nr * sizeof(double));

  strcpy(fname, dir);
  strcat(fname, "/data_logR");
  f = fopen(fname, "r");
  for (i = 0; i < data->nr; i++)
    {
      fscanf(f, "%f ", &tmp);
      data->logR[i] = (double)tmp;
    }
  fclose(f);

  strcpy(fname, dir);
  strcat(fname, "/data_logT");
  f = fopen(fname, "r");
  for (i = 0; i < data->nt; i++)
    {
      fscanf(f, "%f ", &tmp);
      data->logT[data->nt - i - 1] = (double)tmp;
    }
  fclose(f);

  strcpy(fname, dir);
  strcat(fname, "/data_logkappa");
  f = fopen(fname, "r");
  for (i = 0; i < data->nt; i++)
    {
      for (j = 0; j < data->nr; j++)
        {
          fscanf(f, "%f ", &tmp);
          //data->logkappa[data->nt - i - 1][j] = (double)tmp;
          data->logkappa[(data->nt - i - 1)*data->nr + j] = (double)tmp;
        }
      fscanf(f, "\n");
    }
  fclose(f);

  data->inter = gsl_interp2d_alloc (gsl_interp2d_bilinear, data->nr, data->nt);
  data->xacc = gsl_interp_accel_alloc();
  data->yacc = gsl_interp_accel_alloc();

  status = gsl_interp2d_init(data->inter, (double *)data->logR, (double *)data->logT, 
      (double *)data->logkappa, data->nr, data->nt);
}

double get_fa_opacity(double rho, double t, opacitydata *data)
{
  double result;
  int status;
  double logt, t6, logr;

  logt = log10(t);
  t6 = t/1e6;
  logr = log10(rho) - 3*log10(t6);


  if (logr < data->logR[0] || logr > data->logR[data->nr-1] || 
      logt < data->logT[0] || logt > data->logT[data->nt-1])
      result = 0.0;
  else
    {
      int err;
      err = gsl_interp2d_eval_e(data->inter, (double *)data->logR, (double *)data->logT, 
          (double *)data->logkappa, logr, logt, data->xacc, data->yacc, &result);
      result = pow(10, result);
    }

  return result;
}

void finalize_opacities(opacitydata *data)
{
  free(data->logR);
  free(data->logT);
  free(data->logkappa);

  gsl_interp2d_free(data->inter);
  gsl_interp_accel_free(data->xacc);
  gsl_interp_accel_free(data->yacc);
}

void initialize_opacities(opacitydata *data, char *dir)
{
  init_opal(dir);

  initialize_fa_opacity(data, dir);
}

float opacities(float rho, float t, opacitydata *data)
{
  float xh = 0.7;
  float kappa;

  if (log10(t) > 4)
  {
    opacity_opal(&xh, &rho, &t, &kappa);
    if (kappa != kappa)
      {
        opacity_opal(&xh, &rho, &t, &kappa);
      }
    if (kappa != kappa)
      {
        printf("Error! Opacity is NaN! Setting to zero.\n");
        kappa = 0;
      }
  }
  else
    kappa = get_fa_opacity(rho, t, data);

  return kappa;
}
