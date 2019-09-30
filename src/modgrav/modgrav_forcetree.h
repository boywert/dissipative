/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/forcetree.h
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


int modgrav_force_construct_full_tree();
void modgrav_tree_fifth_force(double pos_x, double pos_y, double pos_z, int no, double asmthfac, double hmin, float * shortrange_table,
                              double * acc_x, double * acc_y, double * acc_z
#ifdef EVALPOTENTIAL
                 , MyLongDouble * pot
#endif
  );
void modgrav_update_node(int no, int father, int * suns);
