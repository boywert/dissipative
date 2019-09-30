/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/shock_finder_fields.h
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

#ifndef SHOCK_FINDER_FIELDS_H
#define SHOCK_FINDER_FIELDS_H

//MyFloat GravFraction; /** < The gravity fraction of the jump 

//fields which are always present
#define SHOCK_FINDER_FIELDS \
  MyFloat ShockSurfaceArea; /**< The shock surface area */ \
  MyFloat ShockDir[3]; /**< Shock direction */ \
  MyFloat Divvel; /** < The velocity divergence of the cell */ \
  int ShockZone;  /**< Flags whether the cell is in a shock zone */ \
  int ShockSurface; /**< Flags whether the cell is part of the shock surface in the shock zone */ \

//fields which are present for shock finding in ideal hydro
#define SHOCK_FINDER_IDEAL_HYDRO_FIELDS \
    MyFloat Temperature; /**< The temperature of the cell */ \
    MyFloat Tgrad[3]; /**< The gradient of the temperature of the cell */ \
    MyFloat VpostShock[3]; /**< the post shock velocity */ \
    MyFloat VpreShock[3]; /**< the pre shock velocity */ \
    MyFloat CpreShock; /**< pre shock sound speed */ \
    MyFloat RhopreShock; /**< pre shock density */ \
    MyFloat PpreShock; /**< pre shock pressure */ \
    MyFloat TpreShock; /**< pre shock temperature */ \
    MyFloat RhopostShock; /**< the post shock density */ \
    MyFloat PpostShock; /**< the post shock pressure */ \

//fields which are present for shock finding with cosmic rays
#define SHOCK_FINDER_CR_FIELDS \
    MyFloat CRpseudoTemperature; /** < The effective temperature */ \
    MyFloat CRpseudoTgrad[3]; /** The gradient of the effective temperature */ \
    MyFloat VpostShock[3]; /**< the post shock velocity */ \
    MyFloat VpreShock[3]; /**< the pre shock velocity */ \
    MyFloat RhopreShock; /**< pre shock density */ \
    MyFloat RhopostShock; /**< the post shock density */ \

//fields which are present for SHOCK_FINDER_BEFORE_OUTPUT_MORE
#define SHOCK_FINDER_FIELDS_MORE \
    MyFloat ShockSurfaceArea; \
    MyFloat ShockDir[3]; \
    MyFloat Divvel; \
    int ZoneFlag; \

#endif
