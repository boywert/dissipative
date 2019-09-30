/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dg/dg_legendre.c
 * \date        10/2014
 * \author		Kevin Schaal
 * \brief		Scaled Legendre polynomials and their derivatives
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

//3d done (not yet interface parts, refinement matrices)

#include "../allvars.h"
#include "../proto.h"

#ifdef DG

#define HORNER_SCHEME

/*!
 * Hard-coded scaled Legendre Polynomials
 * for testing/debugging purposes
 */

double P_0(double x)
{
  return 1;
}

double P_1(double x)
{
  return sqrt(3) * x;
}

double P_2(double x)
{
#ifdef HORNER_SCHEME
  return sqrt(5.) * 0.5 * (-1. + x * (x * 3.));
#else
  return sqrt(5.) * 0.5 * (3. * x * x - 1.);
#endif
}

double P_3(double x)
{
#ifdef HORNER_SCHEME
  return sqrt(7.) * 0.5 * (x * (-3. + x * (x * 5.)));
#else
  return sqrt(7.) * 0.5 * (5. * x * x * x - 3. * x);
#endif
}

double P_4(double x)
{
#ifdef HORNER_SCHEME
  return 3. / 8. * (3. + x * (x * (-30. + x * (x * 35.))));
#else
  return 3. / 8. * (35. * x * x * x * x - 30. * x * x + 3.);
#endif
}

double P_5(double x)
{
#ifdef HORNER_SCHEME
  return sqrt(11.) / 8. * (x * (15 + x * (x * (-70 + x * (x * 63)))));
#else
  return sqrt(11.) / 8. * (63. * x * x * x * x * x - 70. * x * x * x + 15 * x);
#endif
}

/*!
 * Hard-coded derivatives of th scaled Legendre Polynomials
 * for testing/debugging purposes
 */

double dP_0(double x)
{
  return 0;
}

double dP_1(double x)
{
  return sqrt(3.);
}

double dP_2(double x)
{
  return sqrt(5.) * 0.5 * 6. * x;
}

double dP_3(double x)
{
#ifdef HORNER_SCHEME
  return sqrt(7.) * 0.5 * (-3. + x * (x * 15.));
#else
  return sqrt(7.) * 0.5 * (15. * x * x - 3.);
#endif
}

double dP_4(double x)
{
#ifdef HORNER_SCHEME
  return sqrt(9.) / 8. * (x * (-60. + x * (0 + x * 140.)));
#else
  return sqrt(9.) / 8. * (140 * x * x * x - 60 * x);
#endif
}

double dP_5(double x)
{
#ifdef HORNER_SCHEME
  return sqrt(11.) / 8. * (15. + x * (x * (-210. + x * (x * 315.))));
#else
  return sqrt(11.) / 8. * (315. * x * x * x * x - 210 * x * x + 15.);
#endif
}

/*!
 * Calculate the value of Legendre Polynom n (n=0,...) at position x in [-1,1]
 */
double Legendre(int n, double x)
{
  if(n == 0)
    {
      return 1;
    }
  else if(n == 1)
    {
      return x;
    }
  else
    {
      return 1. / n * ((2. * n - 1) * x * Legendre(n - 1, x) - (n - 1.) * Legendre(n - 2, x));
    }
}

/*!
 * Calculate the value of the scaled Legendre Polynom n (n=0,...) at position x in [-1,1]
 */
double scaled_Legendre(int n, double x)
{
  return sqrt(2. * n + 1) * Legendre(n, x);
}

/*!
 * Calculate the value of the derivative of the Legendre Polynom n (n=0,...) at position x in [-1,1]
 */
double deriv_Legendre(int n, double x)
{
  if(n == 0)
    {
      return 0;
    }
  else if(n == 1)
    {
      return 1;
    }
  else
    {
      return 1. / n * ((2. * n - 1) * Legendre(n - 1, x) + (2 * n - 1) * x * deriv_Legendre(n - 1, x) - (n - 1.) * deriv_Legendre(n - 2, x));
    }
}


/*!
 * Calculate the value of the derivative of the scaled Legendre Polynom n (n=0,...) at position x in [-1,1]
 */
double deriv_scaled_Legendre(int n, double x)
{
  return sqrt(2. * n + 1) * deriv_Legendre(n, x);
}

/*!
 * Calculate the root (y=0) of the Legendre polynomial n (n=0,...) close to position x0
 */
static double newton_raphson_legendre(int n, double x0)
{
  double tol = 1e-12;

  double x_n = x0;
  double x_n_plus_1;
  double change = 1;
  int counter = 0;

  while(change > tol)
    {
      x_n_plus_1 = x_n - Legendre(n, x_n) / deriv_Legendre(n, x_n);

      change = 2 * fabs((x_n_plus_1 - x_n) / (x_n_plus_1 + x_n));

      counter++;

      if(counter > 100)
        {
          terminate("Error in newton_raphson_legendre: Newton-Raphson method did not converge!\n");
        }

      x_n = x_n_plus_1;
    }

  return x_n_plus_1;
}

/*!
 * calculate the kth root (k=1,2,...) of the Legendre polynomial n (n>=1, since P0=const).
 * the roots correspond to the 1d quadrature points
 */
double legendre_root(int n, int k)
{
  mpi_fprintf(FdInfoDG, "\nCalculating legendre root %d for polynomial %d\n", k, n);

  double approx_root;
  double root;

  approx_root = (1 - 1. / (8. * n * n) + 1. / (8 * n * n * n)) * cos(M_PI * (4 * k - 1) / (4 * n + 2));

  mpi_fprintf(FdInfoDG, "\tapproximated root:%f\n", approx_root);

  root = newton_raphson_legendre(n, approx_root);

  mpi_fprintf(FdInfoDG, "\troot:%f\n", root);

  return root;
}

/*!
 * calculate the gauss weight for a quadrature point
 */
double gauss_weight(int n, double root)
{
  double denominator_sqrt;

  denominator_sqrt = n * Legendre(n - 1, root);

  return 2 * (1 - root * root) / (denominator_sqrt * denominator_sqrt);
}

/*!
 * calculate the base function indices from the 1d index
 *
 * \param k 1d index for accessing base functions
 * \param Px output parameter, base function in x-direction
 * \param Py output parameter, base function in y-direction
 * \param Pz output parameter, base function in z-direction
 */

#if (NUMDIMS==3)
void index_to_base_function(int k, int *Px, int *Py, int *Pz)
{
  int u, v, w, deg_k;

  int counter = 0;

  for(deg_k = 0; deg_k <= DEGREE_K; deg_k++)
    {
      for(u = 0; u <= deg_k; u++)
        {
          for(v = 0; v <= deg_k - u; v++)
            {
              for(w = 0; w <= deg_k - u - v; w++)
                {
                  if(u + v + w == deg_k)
                    {
                      if(counter == k)
                        {
                          *Px = w;
                          *Py = v;
                          *Pz = u;
                          // printf("\t base function: %d, Px=%d, Py=%d, Pz=%d\n", k, *Px, *Py, *Pz); //kundo
                          return;
                        }
                      else
                        {
                          counter++;
                        }
                    }
                }
            }
        }
    }

  *Px = 0;
  *Py = 0;
  *Pz = 0;

  terminate("Error in index_to_base_function!\n");
  return;
}

#endif

#if (NUMDIMS==2)
void index_to_base_function(int k, int *Px, int *Py, int *Pz)
{
  int degree = 0;
  int counter = 0;

  while(1)
    {
      if(k <= counter)
        {
          break;
        }
      else
        {
          degree++;
          counter += degree + 1;
        }
    }

  *Px = 0;
  *Py = degree;
  *Pz = 0;

  while(k != counter)
    {
      counter--;
      *Px = *Px + 1;
      *Py = *Py - 1;
      *Pz = 0;
    }

  return;
}

#endif

/*!
 * get the value of the base function k (k=0, ... , Nof_base_functions-1) at (x,y)
 */
double base_function(int k, double x, double y, double z)
{
  int Px, Py, Pz;

  index_to_base_function(k, &Px, &Py, &Pz);

#ifdef TWODIMS
  assert(scaled_Legendre(Pz, z) == 1);
#endif

  return scaled_Legendre(Px, x) * scaled_Legendre(Py, y) * scaled_Legendre(Pz, z);
}

/*!
 * get the value of the x-derivative of the base function k (k=0, ... , Nof_base_functions-1) at (x,y,z)
 */
double deriv_x_base_function(int k, double x, double y, double z)
{
  int Px, Py, Pz;

  index_to_base_function(k, &Px, &Py, &Pz);

#ifdef TWODIMS
  assert(scaled_Legendre(Pz, z) == 1);
#endif

  return deriv_scaled_Legendre(Px, x) * scaled_Legendre(Py, y) * scaled_Legendre(Pz, z);
}

/*!
 * get the value of the y-derivative of the base function k (k=0, ... , Nof_base_functions-1) at (x,y,z)
 */
double deriv_y_base_function(int k, double x, double y, double z)
{
  int Px, Py, Pz;

  index_to_base_function(k, &Px, &Py, &Pz);

#ifdef TWODIMS
  assert(scaled_Legendre(Pz, z) == 1);
#endif

  return scaled_Legendre(Px, x) * deriv_scaled_Legendre(Py, y) * scaled_Legendre(Pz, z);
}

/*!
 * get the value of the z-derivative of the base function k (k=0, ... , Nof_base_functions-1) at (x,y,z)
 */
double deriv_z_base_function(int k, double x, double y, double z)
{
  int Px, Py, Pz;

  index_to_base_function(k, &Px, &Py, &Pz);

  return scaled_Legendre(Px, x) * scaled_Legendre(Py, y) * deriv_scaled_Legendre(Pz, z);
}

/*!
 * Calculate the absolute position of an inner quadrature point
 */
void abs_pos_inner_quad_points(int q, double cell_dl, double cell_center[3], double pos_abs[3])
{
  pos_abs[0] = 0.5 * GET_Inner_quad_points_x(q) * cell_dl + cell_center[0];
  pos_abs[1] = 0.5 * GET_Inner_quad_points_y(q) * cell_dl + cell_center[1];
  pos_abs[2] = 0.5 * GET_Inner_quad_points_z(q) * cell_dl + cell_center[2];
}

/*!
 * Calculate the absolute position of an outer quadrature point
 */
void abs_pos_outer_quad_points(int e, int q, double cell_dl, double cell_center[3], double pos_abs[3])
{
  pos_abs[0] = 0.5 * GET_Outer_quad_points_x(e, q) * cell_dl + cell_center[0];
  pos_abs[1] = 0.5 * GET_Outer_quad_points_y(e, q) * cell_dl + cell_center[1];
  pos_abs[2] = 0.5 * GET_Outer_quad_points_z(e, q) * cell_dl + cell_center[2];
}


/*!
 * Initialize the quadrature points (positions) and weights for one dimension
 */
void ini_1d_quadrature_points()
{

  int i;

#ifndef CALC_QUADRATURE_DATA    //hard-coded positions and weights

  switch (Nof_quad_points_1d)
    {
    case 1:
      SET_Quad_points_1d(0, 0);
      SET_Quad_points_weights_1d(0, 2);
      break;

    case 2:
      SET_Quad_points_1d(0, -sqrt(1. / 3.));
      SET_Quad_points_1d(1, sqrt(1. / 3.));
      SET_Quad_points_weights_1d(0, 1);
      SET_Quad_points_weights_1d(1, 1);
      break;

    case 3:
      SET_Quad_points_1d(0, -sqrt(3. / 5.));
      SET_Quad_points_1d(1, 0);
      SET_Quad_points_1d(2, sqrt(3. / 5.));
      SET_Quad_points_weights_1d(0, 5. / 9.);
      SET_Quad_points_weights_1d(1, 8. / 9.);
      SET_Quad_points_weights_1d(2, 5. / 9.);
      break;

    default:
      terminate("Error in ini_1d_quadrature_points: Quadrature data for degree k=%d is not in the table!\n", DEGREE_K);
      break;
    }
#else //calculate the positions and weights


  for(i = 0; i < Nof_quad_points_1d; i++)
    {
      SET_Quad_points_1d(Nof_quad_points_1d - 1 - i, legendre_root(Nof_quad_points_1d, i + 1));
      SET_Quad_points_weights_1d(i, gauss_weight(Nof_quad_points_1d, GET_Quad_points_1d(Nof_quad_points_1d - 1 - i)));
    }

#endif

  //initialize the 1d Lobatto quadrature points and weights (needed for the positivity limiter)
  for(i = 0; i < Nof_lobatto_points_1d; i++)
    {
      switch (Nof_lobatto_points_1d)
        {
        case 2:
          SET_Lobatto_points_1d(0, -1);
          SET_Lobatto_points_1d(1, 1);
          SET_Lobatto_points_weights_1d(0, 1);
          SET_Lobatto_points_weights_1d(1, 1);
          break;

        case 3:
          SET_Lobatto_points_1d(0, -1);
          SET_Lobatto_points_1d(1, 0);
          SET_Lobatto_points_1d(2, 1);
          SET_Lobatto_points_weights_1d(0, 1. / 3.);
          SET_Lobatto_points_weights_1d(1, 4. / 3.);
          SET_Lobatto_points_weights_1d(2, 1. / 3.);
          break;

        case 4:
          SET_Lobatto_points_1d(0, -1);
          SET_Lobatto_points_1d(1, -sqrt(1. / 5.));
          SET_Lobatto_points_1d(2, sqrt(1. / 5.));
          SET_Lobatto_points_1d(3, 1);
          SET_Lobatto_points_weights_1d(0, 1. / 6.);
          SET_Lobatto_points_weights_1d(1, 5. / 6.);
          SET_Lobatto_points_weights_1d(2, 5. / 6.);
          SET_Lobatto_points_weights_1d(3, 1. / 6.);
          break;

        case 5:
          SET_Lobatto_points_1d(0, -1);
          SET_Lobatto_points_1d(1, -sqrt(3. / 7.));
          SET_Lobatto_points_1d(2, 0);
          SET_Lobatto_points_1d(3, sqrt(3. / 7.));
          SET_Lobatto_points_1d(4, 1);
          SET_Lobatto_points_weights_1d(0, 1. / 10.);
          SET_Lobatto_points_weights_1d(1, 49. / 90.);
          SET_Lobatto_points_weights_1d(2, 32. / 45.);
          SET_Lobatto_points_weights_1d(3, 49. / 90.);
          SET_Lobatto_points_weights_1d(4, 1. / 10.);
          break;

        default:
          terminate("Error in ini_1d_quadrature_points: Lobatto quadrature data for degree k=%d is not in the table!\n", DEGREE_K);
          break;
        }
    }
}

void delete_duplicated_union_quad_points()
{
#ifdef OLD_VERSION
  return;
#endif

  mpi_fprintf(FdInfoDG, "\nNumber of union quadrature points: %d\n", Nof_union_quad_points);

  double Union_quad_points_x_temp[Nof_union_quad_points];
  double Union_quad_points_y_temp[Nof_union_quad_points];
  double Union_quad_points_z_temp[Nof_union_quad_points];

  int Nof_union_quad_points_temp = 0;
  int duplicate = 0;

  int counter_duplicates = 0;

  int k, l;

  for(k = 0; k < Nof_union_quad_points; k++)
    {
      //check whether the union quadrature point already exists
      duplicate = 0;

      for(l = 0; l < Nof_union_quad_points_temp; l++)
        {
          if(Union_quad_points_x_temp[l] == GET_Union_quad_points_x(k) && Union_quad_points_y_temp[l] == GET_Union_quad_points_y(k) && Union_quad_points_z_temp[l] == GET_Union_quad_points_z(k))
            {
              counter_duplicates++;
              duplicate = 1;

              mpi_fprintf(FdInfoDG, "Deleting duplicate union quadrature point (%f, %f, %f)\n", GET_Union_quad_points_x(k), GET_Union_quad_points_y(k), GET_Union_quad_points_z(k));
            }
        }

      if(!duplicate)
        {
          Union_quad_points_x_temp[Nof_union_quad_points_temp] = GET_Union_quad_points_x(k);
          Union_quad_points_y_temp[Nof_union_quad_points_temp] = GET_Union_quad_points_y(k);
          Union_quad_points_z_temp[Nof_union_quad_points_temp] = GET_Union_quad_points_z(k);

          Nof_union_quad_points_temp++;
        }
    }



  //new union quadrature points (non-duplicates)
  Nof_union_quad_points = Nof_union_quad_points_temp;

  myfree(Union_base_values);
  myfree(Union_quad_points_z);
  myfree(Union_quad_points_y);
  myfree(Union_quad_points_x);

  Union_quad_points_x = (double *) mymalloc("Union_quad_points_x", Nof_union_quad_points * sizeof(double));
  Union_quad_points_y = (double *) mymalloc("Union_quad_points_y", Nof_union_quad_points * sizeof(double));
  Union_quad_points_z = (double *) mymalloc("Union_quad_points_z", Nof_union_quad_points * sizeof(double));
  Union_base_values = (double (*)[NOF_BASE_FUNCTIONS]) mymalloc("Union_base_values", Nof_base_functions * Nof_union_quad_points * sizeof(double));

  for(k = 0; k < Nof_union_quad_points; k++)
    {
      SET_Union_quad_points_x(k, Union_quad_points_x_temp[k]);
      SET_Union_quad_points_y(k, Union_quad_points_y_temp[k]);
      SET_Union_quad_points_z(k, Union_quad_points_z_temp[k]);
    }

  mpi_fprintf(FdInfoDG, "Found %d duplicates!, New number of union quadrature points: %d\n", counter_duplicates, Nof_union_quad_points);
}

/*!
 * Calculate the full quadrature data with the 1d quadrature points and weights
 */
void ini_multi_d_quadrature_points()
{

  int i, j, k;

  DG3D(int l;) k = 0;

  //union quadrature data
  //


#if (NUMDIMS==2)
  for(i = 0; i < Nof_lobatto_points_1d; i++)
    {
      for(j = 0; j < Nof_quad_points_1d; j++)
        {
          SET_Union_quad_points_x(k, GET_Lobatto_points_1d(i));
          SET_Union_quad_points_y(k, GET_Quad_points_1d(j));
          SET_Union_quad_points_z(k, 0);

          k++;
        }
    }

  for(i = 0; i < Nof_lobatto_points_1d; i++)
    {
      for(j = 0; j < Nof_quad_points_1d; j++)
        {
          SET_Union_quad_points_x(k, GET_Quad_points_1d(j));
          SET_Union_quad_points_y(k, GET_Lobatto_points_1d(i));
          SET_Union_quad_points_z(k, 0);

          k++;
        }
    }
#else
  for(i = 0; i < Nof_lobatto_points_1d; i++)
    {
      for(j = 0; j < Nof_quad_points_1d; j++)
        {
          for(l = 0; l < Nof_quad_points_1d; l++)
            {
              SET_Union_quad_points_x(k, GET_Lobatto_points_1d(i));
              SET_Union_quad_points_y(k, GET_Quad_points_1d(j));
              SET_Union_quad_points_z(k, GET_Quad_points_1d(l));

              k++;
            }
        }
    }

  for(i = 0; i < Nof_quad_points_1d; i++)
    {
      for(j = 0; j < Nof_lobatto_points_1d; j++)
        {
          for(l = 0; l < Nof_quad_points_1d; l++)
            {
              SET_Union_quad_points_x(k, GET_Quad_points_1d(i));
              SET_Union_quad_points_y(k, GET_Lobatto_points_1d(j));
              SET_Union_quad_points_z(k, GET_Quad_points_1d(l));

              k++;
            }
        }
    }

  for(i = 0; i < Nof_quad_points_1d; i++)
    {
      for(j = 0; j < Nof_quad_points_1d; j++)
        {
          for(l = 0; l < Nof_lobatto_points_1d; l++)
            {
              SET_Union_quad_points_x(k, GET_Quad_points_1d(i));
              SET_Union_quad_points_y(k, GET_Quad_points_1d(j));
              SET_Union_quad_points_z(k, GET_Lobatto_points_1d(l));

              k++;
            }
        }
    }
#endif

  assert(k == Nof_union_quad_points);

  delete_duplicated_union_quad_points();

  //quadrature data within the cell
  //

#if (NUMDIMS==2)
  for(i = 0; i < Nof_inner_quad_points; i++)
    {
      SET_Inner_quad_points_x(i, GET_Quad_points_1d(i % (DEGREE_K + 1)));
      SET_Inner_quad_points_y(i, GET_Quad_points_1d(i / (DEGREE_K + 1)));
      SET_Inner_quad_points_z(i, 0);

      SET_Inner_quad_points_weights(i, GET_Quad_points_weights_1d(i % (DEGREE_K + 1)) * GET_Quad_points_weights_1d(i / (DEGREE_K + 1)));
    }
#else
  k = 0;

  for(i = 0; i < Nof_quad_points_1d; i++)
    {
      for(j = 0; j < Nof_quad_points_1d; j++)
        {
          for(l = 0; l < Nof_quad_points_1d; l++)
            {

              SET_Inner_quad_points_x(k, GET_Quad_points_1d(i));
              SET_Inner_quad_points_y(k, GET_Quad_points_1d(j));
              SET_Inner_quad_points_z(k, GET_Quad_points_1d(l));

              SET_Inner_quad_points_weights(k, GET_Quad_points_weights_1d(i) * GET_Quad_points_weights_1d(j) * GET_Quad_points_weights_1d(l));

              k++;
            }
        }
    }
#endif


  //quadrature data at the cell interfaces
  //

#if (NUMDIMS==2)
  for(i = 0; i < 4; i++)        //loop over cell interfaces
    {
      for(j = 0; j < (Nof_quad_points_1d); j++) //loop over 1d quadrature points
        {
          SET_Outer_quad_points_weights(i, j, GET_Quad_points_weights_1d(j));

          switch (i)
            {
            case FRONT:
              SET_Outer_quad_points_y(FRONT, j, 1.);
              SET_Outer_quad_points_x(FRONT, j, GET_Quad_points_1d(j));
              break;

            case BACK:
              SET_Outer_quad_points_y(BACK, j, -1.);
              SET_Outer_quad_points_x(BACK, j, GET_Quad_points_1d(j));
              break;

            case LEFT:
              SET_Outer_quad_points_x(LEFT, j, -1.);
              SET_Outer_quad_points_y(LEFT, j, GET_Quad_points_1d(j));
              break;

            case RIGHT:
              SET_Outer_quad_points_x(RIGHT, j, 1.);
              SET_Outer_quad_points_y(RIGHT, j, GET_Quad_points_1d(j));
              break;

            default:
              terminate("Error in ini_multi_d_quadrature_points: Index out of bounds!\n");
              break;
            }
        }
    }
#else

  for(i = 0; i < 6; i++)        //loop over cell interfaces
    {
      k = 0;

      for(j = 0; j < Nof_quad_points_1d; j++)
        {
          for(l = 0; l < Nof_quad_points_1d; l++)
            {
              SET_Outer_quad_points_weights(i, k, GET_Quad_points_weights_1d(j) * GET_Quad_points_weights_1d(l));

              switch (i)
                {
                case FRONT:
                  SET_Outer_quad_points_y(FRONT, k, 1.);
                  SET_Outer_quad_points_x(FRONT, k, GET_Quad_points_1d(j));
                  SET_Outer_quad_points_z(FRONT, k, GET_Quad_points_1d(l));
                  break;

                case BACK:
                  SET_Outer_quad_points_y(BACK, k, -1.);
                  SET_Outer_quad_points_x(BACK, k, GET_Quad_points_1d(j));
                  SET_Outer_quad_points_z(BACK, k, GET_Quad_points_1d(l));
                  break;

                case LEFT:
                  SET_Outer_quad_points_x(LEFT, k, -1.);
                  SET_Outer_quad_points_y(LEFT, k, GET_Quad_points_1d(j));
                  SET_Outer_quad_points_z(LEFT, k, GET_Quad_points_1d(l));
                  break;

                case RIGHT:
                  SET_Outer_quad_points_x(RIGHT, k, 1.);
                  SET_Outer_quad_points_y(RIGHT, k, GET_Quad_points_1d(j));
                  SET_Outer_quad_points_z(RIGHT, k, GET_Quad_points_1d(l));
                  break;

                case TOP:
                  SET_Outer_quad_points_z(TOP, k, 1.);
                  SET_Outer_quad_points_x(TOP, k, GET_Quad_points_1d(j));
                  SET_Outer_quad_points_y(TOP, k, GET_Quad_points_1d(l));
                  break;

                case BOTTOM:
                  SET_Outer_quad_points_z(BOTTOM, k, -1.);
                  SET_Outer_quad_points_x(BOTTOM, k, GET_Quad_points_1d(j));
                  SET_Outer_quad_points_y(BOTTOM, k, GET_Quad_points_1d(l));
                  break;

                default:
                  terminate("Error in ini_multi_d_quadrature_points: Index out of bounds!\n");
                  break;
                }

              k++;
            }
        }

      assert(k == Nof_outer_quad_points);
    }
#endif

  //fine quadrature data for AMR boundary cells
  //
}


/*
 * Calculate the values of the base functions at the quadrature points
 */
void ini_base_function_values()
{
  int i, j;

  //set base values at the Lobatto-Legendre union quadrature points
  for(j = 0; j < Nof_union_quad_points; j++)
    {
      for(i = 0; i < Nof_base_functions; i++)
        {
          SET_Union_base_values(j, i, base_function(i, GET_Union_quad_points_x(j), GET_Union_quad_points_y(j), GET_Union_quad_points_z(j)));
        }
    }


  //set base values within the cell
  //

  for(j = 0; j < Nof_inner_quad_points; j++)
    {
      for(i = 0; i < Nof_base_functions; i++)
        {
          SET_Inner_base_values(j, i, base_function(i, GET_Inner_quad_points_x(j), GET_Inner_quad_points_y(j), GET_Inner_quad_points_z(j)));
          SET_Inner_base_dx_values(j, i, deriv_x_base_function(i, GET_Inner_quad_points_x(j), GET_Inner_quad_points_y(j), GET_Inner_quad_points_z(j)));
          SET_Inner_base_dy_values(j, i, deriv_y_base_function(i, GET_Inner_quad_points_x(j), GET_Inner_quad_points_y(j), GET_Inner_quad_points_z(j)));
#ifdef TWODIMS
          SET_Inner_base_dz_values(j, i, 0);
#else
          SET_Inner_base_dz_values(j, i, deriv_z_base_function(i, GET_Inner_quad_points_x(j), GET_Inner_quad_points_y(j), GET_Inner_quad_points_z(j)));
#endif
        }
    }

  //set base values at the interfaces
  //

  int k;


  for(k = 0; k < Nof_interfaces; k++)   //loop over interfaces
    {
      for(j = 0; j < Nof_outer_quad_points; j++)        //loop over outer quadrature points
        {
          for(i = 0; i < Nof_base_functions; i++)
            {
              SET_Outer_base_values(k, j, i, base_function(i, GET_Outer_quad_points_x(k, j), GET_Outer_quad_points_y(k, j), GET_Outer_quad_points_z(k, j)));
            }
        }
    }

  //set base function values at the fine quadrature points
  //

  double x = 0;
  double y = 0;
  double z = 0;

  double s1, s2;

  int q, q2, s, l;

  int outer_fine_quad_point;

  int nof_fine_sets = pow(2, NUMDIMS - 1);

  for(i = 0; i < Nof_interfaces; i++)
    {
      outer_fine_quad_point = 0;

      for(s = 0; s < nof_fine_sets; s++)
        {
          for(q = 0; q < Nof_quad_points_1d; q++)
            {
#ifdef TWODIMS
              q2 = q;
#else
              for(q2 = 0; q2 < Nof_quad_points_1d; q2++)
#endif
                {
                  switch (s)
                    {
                    case 0:
                      s1 = -1;
                      s2 = 1;
                      break;

                    case 1:
                      s1 = 1;
                      s2 = -1;
                      break;

                    case 2:
                      s1 = 1;
                      s2 = 1;
                      break;

                    case 3:
                      s1 = -1;
                      s2 = -1;
                      break;

                    default:
                      terminate("bad fine quadrature point set!\n");
                    }

                  switch (i)
                    {
                    case LEFT:
                      x = -1;
                      y = (GET_Quad_points_1d(q) + s1) * 0.5;
                      z = (GET_Quad_points_1d(q2) + s2) * 0.5;
                      break;

                    case RIGHT:
                      x = 1;
                      y = (GET_Quad_points_1d(q) + s1) * 0.5;
                      z = (GET_Quad_points_1d(q2) + s2) * 0.5;
                      break;

                    case BACK:
                      x = (GET_Quad_points_1d(q) + s1) * 0.5;
                      y = -1;
                      z = (GET_Quad_points_1d(q2) + s2) * 0.5;
                      break;

                    case FRONT:
                      x = (GET_Quad_points_1d(q) + s1) * 0.5;
                      y = 1;
                      z = (GET_Quad_points_1d(q2) + s2) * 0.5;
                      break;

                    case BOTTOM:
                      x = (GET_Quad_points_1d(q) + s1) * 0.5;
                      y = (GET_Quad_points_1d(q2) + s2) * 0.5;
                      z = -1;
                      break;

                    case TOP:
                      x = (GET_Quad_points_1d(q) + s1) * 0.5;
                      y = (GET_Quad_points_1d(q2) + s2) * 0.5;
                      z = 1;
                      break;

                    default:
                      terminate("bad interface!\n");
                    }

#ifdef OUTER_FINE_QUAD_POINT_INFO
#ifdef TWODIMS
                  mpi_fprintf(FdInfoDG, "fine quad point (%f, %f)\n", x, y);
#else
                  mpi_fprintf(FdInfoDG, "fine quad point (%f, %f, %f)\n", x, y, z);
#endif
#endif

                  for(l = 0; l < Nof_base_functions; l++)
                    {
                      SET_Outer_fine_base_values(i, outer_fine_quad_point, l, base_function(l, x, y, z));
#ifdef OUTER_FINE_QUAD_POINT_BASE_INFO
                      mpi_fprintf(FdInfoDG, "\t base fct %d value: %f\n", l, base_function(l, x, y, z));
#endif
                    }

                  outer_fine_quad_point++;
                }
            }                   //end quad point loop
        }                       //end set loop

      assert(outer_fine_quad_point == Nof_outer_fine_quad_points);

    }                           //end interface loop
}

/*!
 * Calculate the values for the refinement/derefinement matrices
 */
void ini_refinement_matrices()
{
  int j, l;

  int u, v;
  int w = 0;

  double sum_a, sum_b, sum_c, sum_d DG3D(, sum_e, sum_f, sum_g, sum_h);

  for(j = 0; j < Nof_base_functions; j++)
    {
      for(l = 0; l < Nof_base_functions; l++)
        {
          sum_a = 0;
          sum_b = 0;
          sum_c = 0;
          sum_d = 0;

#ifndef TWODIMS
          sum_e = 0;
          sum_f = 0;
          sum_g = 0;
          sum_h = 0;
#endif

          //gaussian quadrature
          for(u = 0; u < Nof_quad_points_1d; u++)
            {
              for(v = 0; v < Nof_quad_points_1d; v++)
                {
                  DG3D(for(w = 0; w < Nof_quad_points_1d; w++))
                    {
                      sum_a += base_function(j, (GET_Quad_points_1d(u) - 1) * 0.5, (GET_Quad_points_1d(v) - 1) * 0.5, (GET_Quad_points_1d(w) - 1) * 0.5)
                        * base_function(l, GET_Quad_points_1d(u), GET_Quad_points_1d(v), GET_Quad_points_1d(w))
                        * GET_Quad_points_weights_1d(u) * GET_Quad_points_weights_1d(v) DG3D(*GET_Quad_points_weights_1d(w));

                      sum_b += base_function(j, (GET_Quad_points_1d(u) + 1) * 0.5, (GET_Quad_points_1d(v) - 1) * 0.5, (GET_Quad_points_1d(w) - 1) * 0.5)
                        * base_function(l, GET_Quad_points_1d(u), GET_Quad_points_1d(v), GET_Quad_points_1d(w))
                        * GET_Quad_points_weights_1d(u) * GET_Quad_points_weights_1d(v) DG3D(*GET_Quad_points_weights_1d(w));

                      sum_c += base_function(j, (GET_Quad_points_1d(u) - 1) * 0.5, (GET_Quad_points_1d(v) + 1) * 0.5, (GET_Quad_points_1d(w) - 1) * 0.5)
                        * base_function(l, GET_Quad_points_1d(u), GET_Quad_points_1d(v), GET_Quad_points_1d(w))
                        * GET_Quad_points_weights_1d(u) * GET_Quad_points_weights_1d(v) DG3D(*GET_Quad_points_weights_1d(w));

                      sum_d += base_function(j, (GET_Quad_points_1d(u) + 1) * 0.5, (GET_Quad_points_1d(v) + 1) * 0.5, (GET_Quad_points_1d(w) - 1) * 0.5)
                        * base_function(l, GET_Quad_points_1d(u), GET_Quad_points_1d(v), GET_Quad_points_1d(w))
                        * GET_Quad_points_weights_1d(u) * GET_Quad_points_weights_1d(v) DG3D(*GET_Quad_points_weights_1d(w));

#ifndef TWODIMS
                      sum_e += base_function(j, (GET_Quad_points_1d(u) - 1) * 0.5, (GET_Quad_points_1d(v) - 1) * 0.5, (GET_Quad_points_1d(w) + 1) * 0.5)
                        * base_function(l, GET_Quad_points_1d(u), GET_Quad_points_1d(v), GET_Quad_points_1d(w))
                        * GET_Quad_points_weights_1d(u) * GET_Quad_points_weights_1d(v) DG3D(*GET_Quad_points_weights_1d(w));

                      sum_f += base_function(j, (GET_Quad_points_1d(u) + 1) * 0.5, (GET_Quad_points_1d(v) - 1) * 0.5, (GET_Quad_points_1d(w) + 1) * 0.5)
                        * base_function(l, GET_Quad_points_1d(u), GET_Quad_points_1d(v), GET_Quad_points_1d(w))
                        * GET_Quad_points_weights_1d(u) * GET_Quad_points_weights_1d(v) DG3D(*GET_Quad_points_weights_1d(w));

                      sum_g += base_function(j, (GET_Quad_points_1d(u) - 1) * 0.5, (GET_Quad_points_1d(v) + 1) * 0.5, (GET_Quad_points_1d(w) + 1) * 0.5)
                        * base_function(l, GET_Quad_points_1d(u), GET_Quad_points_1d(v), GET_Quad_points_1d(w))
                        * GET_Quad_points_weights_1d(u) * GET_Quad_points_weights_1d(v) DG3D(*GET_Quad_points_weights_1d(w));

                      sum_h += base_function(j, (GET_Quad_points_1d(u) + 1) * 0.5, (GET_Quad_points_1d(v) + 1) * 0.5, (GET_Quad_points_1d(w) + 1) * 0.5)
                        * base_function(l, GET_Quad_points_1d(u), GET_Quad_points_1d(v), GET_Quad_points_1d(w))
                        * GET_Quad_points_weights_1d(u) * GET_Quad_points_weights_1d(v) DG3D(*GET_Quad_points_weights_1d(w));
#endif
                    }
                }
            }

#ifdef OLD_VERSION
          sum_a *= 0.25;
          sum_b *= 0.25;
          sum_c *= 0.25;
          sum_d *= 0.25;
#ifndef TWODIMS
          sum_e *= 0.25;
          sum_f *= 0.25;
          sum_g *= 0.25;
          sum_h *= 0.25;
#endif
#else
          sum_a /= DG_PROJ_NORM;
          sum_b /= DG_PROJ_NORM;
          sum_c /= DG_PROJ_NORM;
          sum_d /= DG_PROJ_NORM;
#ifndef TWODIMS
          sum_e /= DG_PROJ_NORM;
          sum_f /= DG_PROJ_NORM;
          sum_g /= DG_PROJ_NORM;
          sum_h /= DG_PROJ_NORM;
#endif
#endif

          SET_P_A(j, l, sum_a);
          SET_P_B(j, l, sum_b);
          SET_P_C(j, l, sum_c);
          SET_P_D(j, l, sum_d);

#ifndef TWODIMS
          SET_P_E(j, l, sum_e);
          SET_P_F(j, l, sum_f);
          SET_P_G(j, l, sum_g);
          SET_P_H(j, l, sum_h);
#endif
        }
    }
}


/*!
 *  Print the quadrature point and weights
 */
void print_quad_info()
{
  int i, j;

  mpi_fprintf(FdInfoDG, "\n=== Discontinuous Galerkin quadrature info ===\n\n");

  mpi_fprintf(FdInfoDG, "1d Labatto-quadrature points and weights:\n\n");

  for(j = 0; j < Nof_lobatto_points_1d; j++)
    {
      mpi_fprintf(FdInfoDG, "x=%f\tw=%f\n", GET_Lobatto_points_1d(j), GET_Lobatto_points_weights_1d(j));
    }

  mpi_fprintf(FdInfoDG, "\n1d Gauss quadrature points and weights:\n\n");

  for(j = 0; j < Nof_quad_points_1d; j++)
    {
      mpi_fprintf(FdInfoDG, "x=%f\tw=%f\n", GET_Quad_points_1d(j), GET_Quad_points_weights_1d(j));
    }

  mpi_fprintf(FdInfoDG, "\nGauss U Lobatto quadrature points (union):\n\n");

  for(j = 0; j < Nof_union_quad_points; j++)
    {
#ifdef TWODIMS
      mpi_fprintf(FdInfoDG, "x=%f\ty=%f\n", GET_Union_quad_points_x(j), GET_Union_quad_points_y(j));
#else
      mpi_fprintf(FdInfoDG, "x=%f\ty=%f\tz=%f\n", GET_Union_quad_points_x(j), GET_Union_quad_points_y(j), GET_Union_quad_points_z(j));
#endif
    }

  mpi_fprintf(FdInfoDG, "\nGauss quadrature points and weights:\n\n");

  for(j = 0; j < Nof_inner_quad_points; j++)
    {
#ifdef TWODIMS
      mpi_fprintf(FdInfoDG, "x=%f\ty=%f\tw=%f\n", GET_Inner_quad_points_x(j), GET_Inner_quad_points_y(j), GET_Inner_quad_points_weights(j));
#else
      mpi_fprintf(FdInfoDG, "x=%f\ty=%f\tz=%f\tw=%f\n", GET_Inner_quad_points_x(j), GET_Inner_quad_points_y(j), GET_Inner_quad_points_z(j), GET_Inner_quad_points_weights(j));
#endif
    }

  mpi_fprintf(FdInfoDG, "\ninterface Gauss quadrature points and weights:\n\n");

  for(i = 0; i < Nof_interfaces; i++)
    {
      for(j = 0; j < Nof_outer_quad_points; j++)
        {
#ifdef TWODIMS
          mpi_fprintf(FdInfoDG, "x=%f\ty=%f\tw=%f\n", GET_Outer_quad_points_x(i, j), GET_Outer_quad_points_y(i, j), GET_Outer_quad_points_weights(i, j));
#else
          mpi_fprintf(FdInfoDG, "x=%f\ty=%f\tz=%f\tw=%f\n", GET_Outer_quad_points_x(i, j), GET_Outer_quad_points_y(i, j), GET_Outer_quad_points_z(i, j), GET_Outer_quad_points_weights(i, j));
#endif
        }
    }
}

/*!
 * Print the values of the base functions at the quadrature points
 */
void print_base_info()
{
  mpi_fprintf(FdInfoDG, "\n=== base function info ===\n\n");

  int i, j, k;

  int Px, Py, Pz;

  mpi_fprintf(FdInfoDG, "Lobatto values:\n\n");

  for(j = 0; j < Nof_union_quad_points; j++)
    {
      for(i = 0; i < Nof_base_functions; i++)
        {
          index_to_base_function(i, &Px, &Py, &Pz);

#ifdef TWODIMS
          mpi_fprintf(FdInfoDG, "base function: %d: P%d(x)*P%d(y), quad_point: (%f, %f), value: %f\n", i, Px, Py, GET_Union_quad_points_x(j), GET_Union_quad_points_y(j), GET_Union_base_values(j, i));
#else
          mpi_fprintf(FdInfoDG, "base function: %d: P%d(x)*P%d(y)*P%d(z), quad_point: (%f, %f, %f), value: %f\n", i, Px, Py, Pz, GET_Union_quad_points_x(j), GET_Union_quad_points_y(j),
                      GET_Union_quad_points_z(j), GET_Union_base_values(j, i));
#endif
        }
    }


  mpi_fprintf(FdInfoDG, "\nCell values:\n\n");

  for(j = 0; j < Nof_inner_quad_points; j++)
    {
      for(i = 0; i < Nof_base_functions; i++)
        {

          index_to_base_function(i, &Px, &Py, &Pz);

#ifdef TWODIMS
          mpi_fprintf(FdInfoDG, "base function: %d: P%d(x)*P%d(y), quad_point: (%f, %f), value: %f, d/dx value: %f, d/dy value: %f\n",
                      i, Px, Py, GET_Inner_quad_points_x(j), GET_Inner_quad_points_y(j), GET_Inner_base_values(j, i), GET_Inner_base_dx_values(j, i), GET_Inner_base_dy_values(j, i));
#else
          mpi_fprintf(FdInfoDG, "base function: %d: P%d(x)*P%d(y)*P%d(z), quad_point: (%f, %f, %f), value: %f, d/dx value: %f, d/dy value: %f, d/dz value: %f\n",
                      i, Px, Py, Pz, GET_Inner_quad_points_x(j), GET_Inner_quad_points_y(j), GET_Inner_quad_points_z(j), GET_Inner_base_values(j, i), GET_Inner_base_dx_values(j, i),
                      GET_Inner_base_dy_values(j, i), GET_Inner_base_dz_values(j, i));
#endif
        }
    }


  mpi_fprintf(FdInfoDG, "\nInterface values:\n\n");


  for(k = 0; k < Nof_interfaces; k++)   //loop over interfaces
    {
      for(j = 0; j < (Nof_outer_quad_points); j++)
        {
          for(i = 0; i < Nof_base_functions; i++)
            {
              index_to_base_function(i, &Px, &Py, &Pz);

#ifdef TWODIMS
              mpi_fprintf(FdInfoDG, "base function: %d: P%d(x)*P%d(y), quad_point: (%f, %f), value: %f\n",
                          i, Px, Py, GET_Outer_quad_points_x(k, j), GET_Outer_quad_points_y(k, j), GET_Outer_base_values(k, j, i));
#else
              mpi_fprintf(FdInfoDG, "base function: %d: P%d(x)*P%d(y)*P%d(z), quad_point: (%f, %f, %f), value: %f\n",
                          i, Px, Py, Pz, GET_Outer_quad_points_x(k, j), GET_Outer_quad_points_y(k, j), GET_Outer_quad_points_z(k, j), GET_Outer_base_values(k, j, i));
#endif
            }
        }
    }
}

/*!
 * Print the values of the refinement/derefinement matrices
 */
void print_ref_matrix_info()
{
  mpi_fprintf(FdInfoDG, "\n=== refinement matrix info ===\n\n");

  int j, l;

  for(j = 0; j < Nof_base_functions; j++)
    {
      for(l = 0; l < Nof_base_functions; l++)
        {
          mpi_fprintf(FdInfoDG, "P_A[%d][%d]=%f\n", j, l, GET_P_A(j, l));
        }
    }

  mpi_fprintf(FdInfoDG, "\n");

  for(j = 0; j < Nof_base_functions; j++)
    {
      for(l = 0; l < Nof_base_functions; l++)
        {
          mpi_fprintf(FdInfoDG, "P_B[%d][%d]=%f\n", j, l, GET_P_B(j, l));
        }
    }

  mpi_fprintf(FdInfoDG, "\n");

  for(j = 0; j < Nof_base_functions; j++)
    {
      for(l = 0; l < Nof_base_functions; l++)
        {
          mpi_fprintf(FdInfoDG, "P_C[%d][%d]=%f\n", j, l, GET_P_C(j, l));
        }
    }

  mpi_fprintf(FdInfoDG, "\n");

  for(j = 0; j < Nof_base_functions; j++)
    {
      for(l = 0; l < Nof_base_functions; l++)
        {
          mpi_fprintf(FdInfoDG, "P_D[%d][%d]=%f\n", j, l, GET_P_D(j, l));
        }
    }

#ifndef TWODIMS
  mpi_fprintf(FdInfoDG, "\n");

  for(j = 0; j < Nof_base_functions; j++)
    {
      for(l = 0; l < Nof_base_functions; l++)
        {
          mpi_fprintf(FdInfoDG, "P_E[%d][%d]=%f\n", j, l, GET_P_E(j, l));
        }
    }

  mpi_fprintf(FdInfoDG, "\n");

  for(j = 0; j < Nof_base_functions; j++)
    {
      for(l = 0; l < Nof_base_functions; l++)
        {
          mpi_fprintf(FdInfoDG, "P_F[%d][%d]=%f\n", j, l, GET_P_F(j, l));
        }
    }

  mpi_fprintf(FdInfoDG, "\n");

  for(j = 0; j < Nof_base_functions; j++)
    {
      for(l = 0; l < Nof_base_functions; l++)
        {
          mpi_fprintf(FdInfoDG, "P_G[%d][%d]=%f\n", j, l, GET_P_G(j, l));
        }
    }

  mpi_fprintf(FdInfoDG, "\n");

  for(j = 0; j < Nof_base_functions; j++)
    {
      for(l = 0; l < Nof_base_functions; l++)
        {
          mpi_fprintf(FdInfoDG, "P_H[%d][%d]=%f\n", j, l, GET_P_H(j, l));
        }
    }
#endif
}


#endif /* DG */
