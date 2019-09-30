/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        
 * \date        
 * \author      
 * \brief
 * \details     Parallel logging to single file
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#ifndef PARALLEL_LOGS_H
#define PARALLEL_LOGS_H

#include "allvars.h"
#include "proto.h"

#define PLOG_LOCAL_FEEDBACK  1

struct plog_data {
  char *text;
  int *text_lengths;
  int text_count;
  int offset;
  int size;
  int count;
  int extendable;
};

void parallel_log_allocate( struct plog_data *data, int extendable );
void pl_fprintf( struct plog_data *data, const char *fmt, ... );
void parallel_log_write( struct plog_data *data, char *text );
void parallel_log_clear( struct plog_data *data, int type );
void parallel_log_clear_bunch( struct plog_data *data, int type );
void parallel_log_free( struct plog_data *data );

#endif

