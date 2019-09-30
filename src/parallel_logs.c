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

#include <stdarg.h>
#include "allvars.h"
#include "proto.h"
#include "parallel_logs.h"

#define INITIAL_TEXT_SIZE 100000
#define INITIAL_TEXT_COUNT 1000

void parallel_log_allocate( struct plog_data *data, int extendable )
{
  data->offset = 0;
  data->extendable = extendable;
  data->size = INITIAL_TEXT_SIZE;
  data->count = INITIAL_TEXT_COUNT;
  data->text_count = 0;
  
  if(data->extendable)
    {
      data->text = mymalloc_movable( &data->text, "pLogText", data->size );
      data->text_lengths = mymalloc_movable( &data->text_lengths, "pLogTextLengths", data->count * sizeof(int) );
    }
  else
    {
      data->text = mymalloc( "pLogText", data->size );
      data->text_lengths = mymalloc( "pLogTextLengths", data->count * sizeof(int) );
    }
}

void parallel_log_write( struct plog_data *data, char *text )
{
  int len = strlen( text );
  
  if(data->offset+len > data->size || data->text_count == data->count)
    {
      if(!data->extendable)
        {
          printf( "Maximum log size reached, use extendable log, clear more often, or increase constants in parallel_logs.c\n" );
        }
      else
        {
          data->size += INITIAL_TEXT_SIZE;
          data->count += INITIAL_TEXT_COUNT;
          
          data->text = myrealloc_movable( data->text, data->size );
          data->text_lengths = myrealloc_movable( data->text_lengths, data->count * sizeof(int) );
        }
    }
  
  strncpy( &data->text[data->offset], text, len );
  data->text_lengths[data->text_count] = len;
  data->text_count++;
  data->offset += len;
}

void pl_fprintf( struct plog_data *data, const char *format, ... )
{
  char buf[1000];
  va_list args;
  va_start( args, format );
  vsprintf( buf, format, args );
  va_end( args );
  
  parallel_log_write( data, buf );
}

void parallel_log_clear( struct plog_data *data, int type )
{
  if(ThisTask == 0)
    {
      parallel_log_clear_bunch( data, type );
      for(int task=1; task<NTask; task++)
        {
          MPI_Status status;
          
          int text_size, text_count;
          MPI_Recv( &text_size, 1, MPI_INT, task, TAG_PLOGS, MPI_COMM_WORLD, &status );
          
          if(text_size > 0)
            {
              MPI_Recv( &text_count, 1, MPI_INT, task, TAG_PLOGS, MPI_COMM_WORLD, &status );
          
              struct plog_data temp;
              temp.text = mymalloc( "text", text_size );
              temp.text_lengths = mymalloc( "textlengths", text_count * sizeof(int) );
              temp.text_count = text_count;
          
              MPI_Recv( temp.text, text_size, MPI_BYTE, task, TAG_PLOGS, MPI_COMM_WORLD, &status );
              MPI_Recv( temp.text_lengths, text_count, MPI_INT, task, TAG_PLOGS, MPI_COMM_WORLD, &status );
              
              parallel_log_clear_bunch( &temp, type );
              
              myfree( temp.text_lengths );
              myfree( temp.text );
            }
        }
    }
  else
    {
      MPI_Send( &data->offset, 1, MPI_INT, 0, TAG_PLOGS, MPI_COMM_WORLD );
      
      if(data->offset > 0)
        {
          MPI_Send( &data->text_count, 1, MPI_INT, 0, TAG_PLOGS, MPI_COMM_WORLD );
          MPI_Send( data->text, data->offset, MPI_BYTE, 0, TAG_PLOGS, MPI_COMM_WORLD );
          MPI_Send( data->text_lengths, data->text_count, MPI_INT, 0, TAG_PLOGS, MPI_COMM_WORLD );
        }
    }
  
  data->text_count = 0;
  data->offset = 0;
}

void parallel_log_clear_bunch( struct plog_data *data, int type )
{
  FILE *fp;
  switch(type)
    {
#ifdef LOCAL_FEEDBACK
      case PLOG_LOCAL_FEEDBACK:
        fp = FdLocalFeedback;
        break;
#endif
      default:
        terminate( "Type not implemented." );
        break;
    }
  
  char line[1001];
  
  int offset = 0;
  for(int i=0; i<data->text_count; i++)
    {
      int len = imin( data->text_lengths[i], 1000 );
      strncpy( line, &data->text[offset], len );
      line[len] = 0;
      fprintf( fp, line );
      offset += data->text_lengths[i];
    }
}

void parallel_log_free( struct plog_data *data )
{
  myfree( data->text_lengths );
  myfree( data->text );
}

