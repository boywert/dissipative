/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/domain_checks.c
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

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"
#include "domain.h"
#include "voronoi.h"

#ifdef TRACER_MC
static long long check_total_tracers;
static long long check_total_bh_tracers;
static long long check_total_star_tracers;
#endif

void domain_prechecks(void)
{
#ifdef TRACER_MC
  /* tracer counts by particle type should be conserved during domain decomposition */
  check_total_tracers = get_total_number_of_tracers(-1);
  check_total_star_tracers = get_total_number_of_tracers(4);
  check_total_bh_tracers = get_total_number_of_tracers(5);
#endif
}


void domain_post_checks(void)
{
#ifdef TRACER_MC
  long long new_global_tr_count = get_total_number_of_tracers(-1);
  if(check_total_tracers != new_global_tr_count)
    terminate("TRACER_MC: unconserved tracers during domain decomposition: total number of tracers BEFORE = %lld, AFTER = %lld\n", check_total_tracers, new_global_tr_count);

  if(check_total_star_tracers != get_total_number_of_tracers(4))
    terminate("TRACER_MC: strange, global number of STAR tracers decreased in domain, BEFORE = %lld, AFTER = %lld\n", check_total_star_tracers, get_total_number_of_tracers(4));

  if(check_total_bh_tracers != get_total_number_of_tracers(5))
    terminate("TRACER_MC: strange, global number of BH tracers decreased in domain, BEFORE = %lld, AFTER = %lld\n", check_total_bh_tracers, get_total_number_of_tracers(5));
#endif
}
