/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/main.c
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
#include <math.h>
#include <gsl/gsl_math.h>


#include "allvars.h"
#include "proto.h"

#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

/*! \file main.c
 *  \brief Start of the program
 */
/*! \brief The entry point of the program.
 *
 *  This function initializes the MPI communication packages, and sets
 *  cpu-time counters to 0. Then begrun1() is called, which sets up
 *  the simulation. Then either IC's or restart files are loaded. In
 *  case of IC's init() is called which prepares the IC's for the run.
 *  A call to begrun2() finishes the initialization. Finally, run() is
 *  started, the main simulation loop, which iterates over the timesteps.
 */
int main(int argc, char **argv)
{
#ifdef IMPOSE_PINNING
  detect_topology();
  get_core_set();
#endif

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);

  /* output a welcome message */
  hello();

#ifdef EXPLICIT_VECTORIZATION
  if(instrset_detect() < 7)
    {
      mpi_printf("You compiled with explicit vectorization in terms of AVX instructions, but this CPU does not support AVX or higher.\n\n");
      endrun();
    }
#endif

  /* initialize CPU-time/Wallclock-time measurement */
  init_cpu_log();

  determine_compute_nodes();

#ifdef IMPOSE_PINNING
  /* pin the MPI ranks to the available core set */
  pin_to_core_set();
#endif

#if (NUM_THREADS > 1)
  char *p = getenv("OMP_NUM_THREADS");
  int num_threads;
  if(p)
    num_threads = atoi(p);
  else
    num_threads = 1;

  if(num_threads != NUM_THREADS)
    terminate("\n\nYou have activated NUM_THREADS, but the value of the environment\n" "variable OMP_NUM_THREADS=%d is not equal to NUM_THREADS=%d!\n\n", num_threads, NUM_THREADS);

  omp_set_num_threads(NUM_THREADS);
  omp_set_dynamic(0);           /* explicitly turn off dynamic threads */
  mpi_printf("\nUsing %d OpenMP threads\n", omp_get_max_threads());
#ifdef IMPOSE_PINNING
  set_pinning_openmp_threads();
  report_pinning_openmp_threads();
#endif
#else
  mpi_printf("\nPINNING: We are not using OpenMP.\n\n");
#ifdef IMPOSE_PINNING
  report_pinning();
#endif
#endif

#ifdef HOST_MEMORY_REPORTING
  mpi_report_committable_memory();
#endif

  Argc = argc;
  Argv = argv;

  for(PTask = 0; NTask > (1 << PTask); PTask++);

  begrun0();

  if(argc < 2)
    {
      if(ThisTask == 0)
        {
          printf("\nParameters are missing. \n");
          printf("Call with <ParameterFile> [<RestartFlag>] [<RestartSnapNum>] [<SpecialOptions>]\n");
          printf("\n");
          printf("   RestartFlag    Action\n");
          printf("       0          Read initial conditions and start simulation\n");
          printf("       1          Read restart files and resume simulation\n");
          printf("       2          Restart from specified snapshot dump and resume simulation\n");
          printf("       3          Run FOF and optionally SUBFIND: [<SubboxSnapNum> for SUBBOX_SNAPSHOTS]\n");
          printf("       4          Make an image slice:    <SnapNum> <pixelsX> <pixelsY> <axisX> <axisY> <axisZ> <xmin> <xmax> <ymin> <ymax> <zval>\n");
          printf("       5          Make a projected image: <SnapNum> <pixelsX> <pixelsY> <axisX> <axisY> <axisZ> <xmin> <xmax> <ymin> <ymax> <zmin> <zmax> [<SubboxSnapNum> for SUBBOX_SNAPSHOTS]\n");
          printf("       6          Convert snapshot file to different format [input=ICFormat  output=SnapFormat   NOTE: derived quantities have round-off errors!\n");
          printf("       7          Calculate a velocity power spectrum for the gas cells\n");
          printf("       8          Make a grid projection: <SnapNum> <pixelsX> <pixelsY> <pixelsZ> 0 0  <xmin> <xmax> <ymin> <ymax> <zmin> <zmax>\n");
          printf("       9          Make a projection along an arbitrary axis: <SnapNum> <pixelsX> <pixelsY> <centerX> <centerY> <centerZ> <dirX> <dirY> <dirZ> <boxX> <boxY> <boxZ>\n");
          printf("      10          Make a perspective camera projection: <SnapNum> <pixelsX> <pixelsY> <filename of camera file> [<SubboxSnapNum> for SUBBOX_SNAPSHOTS]\n");
          printf("      11          Calculate power spectra of various quantities for TRACER_PARTICLEs\n");
          printf("      12          Calculate two-point correlation function: <SnapNum> <parttype bitmask> [output path]\n");
          printf("      13          Calculate power spectrum: <SnapNum> <parttype bitmask> [output path]\n");
          printf("      14          Write out the Voronoi mesh: <SnapNum>\n");
          printf("      15          Run the post-processing shock finder: <SnapNum> [<SubboxNum> for SUBBOX_SNAPSHOTS]\n");
          printf("      16          Write out a two-dimensional slice of the Voronoi mesh: <SnapNum> <center_x> <center_y> <center_z> <normal_x> <normal_y> <normal_z>\n");
          printf("      17          Write out snapshot dump with measured gradients\n");
          printf("      18          Recalculate gravitational potential values for specified snaphot dump: <snapnum>\n");
          printf("      19          Calculate additional quantities from a snapshot dump: <snapnum>\n");
          printf("\n");
        }
      endrun();
    }

  strcpy(ParameterFile, argv[1]);

  if(argc >= 3)
    RestartFlag = atoi(argv[2]);
  else
    RestartFlag = 0;

  if(argc >= 4)
    RestartSnapNum = atoi(argv[3]);
  else
    RestartSnapNum = -1;

  // Do minimal validation of arguments here rather than in random places in the code
  if((RestartFlag == 3 || RestartFlag == 4 || RestartFlag == 5 || RestartFlag == 6 ||
      RestartFlag == 7 || RestartFlag == 9 || RestartFlag == 11 || RestartFlag == 12 || RestartFlag == 14 || RestartFlag == 15 || RestartFlag == 16 || RestartFlag == 17 || RestartFlag == 18 || RestartFlag == 19) && RestartSnapNum < 0)
    {
      mpi_printf("Need to give the snapshot number\n");
      return (0);
    }

#ifndef RECOMPUTE_POTENTIAL_IN_SNAPSHOT
 if(RestartFlag == 18)
   {
     mpi_printf("Need RECOMPUTE_POTENTIAL_IN_SNAPSHOT for this option\n");
     return (0);
   }
#endif

#if (GFM_STELLAR_EVOLUTION == 2)
  test_stellar_evolution();     /* call only the stellar evolution test routine */
  MPI_Finalize();               /* clean up & finalize MPI */
  return 0;
#endif

#ifdef TEST_COOLING_METAL
  test_cooling_function();      /* call only the metal cooling test routine */
  test_cooling();
  MPI_Finalize();               /* clean up & finalize MPI */
  return 0;
#endif


#ifdef TEST_SFR
  test_sfr();
  MPI_Finalize();               /* clean up & finalize MPI */
  return 0;
#endif

#ifdef RUNNING_SAFETY_FILE
  /* do not run if 'running' safety file exists */
  int runningflag = 0;
  if(ThisTask == 0)
    {
      FILE *fd;
      char runningfname[MAXLEN_PATH];

      sprintf(runningfname, "./running");
      if((fd = fopen(runningfname, "r")))       /* Is the running-file present? If yes, interrupt the run. */
        {
          fclose(fd);
          printf("running-file detected. stopping.\n");
          runningflag = 1;
        }
    }
  MPI_Bcast(&runningflag, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if(runningflag)
    {
      MPI_Finalize(); /* do not call endrun() */
      return 0;
    }
  else
    {
      /* touch a running safety file */
      if(ThisTask == 0)
        {
          FILE *fd;
          char runningfname[MAXLEN_PATH];

          sprintf(runningfname, "./running");
          if((fd = fopen(runningfname, "w")))
            {
              fclose(fd);
              printf("touching a running-file: %s \n", runningfname);
            }
          else
            terminate("could not touch a running-file: %s\n", runningfname);
        }
    }
#endif

  begrun1();                    /* set-up run  */

  /* see if we are loading a restart file or an IC file */
  if(RestartFlag == 1)
    loadrestart();
  else
    {
      /* We're reading an IC file. Is it a snapshot or really an IC? */
      char fname[MAXLEN_PATH];

      if(RestartFlag >= 2 && RestartSnapNum >= 0)
        {

          if(All.NumFilesPerSnapshot > 1)
            sprintf(fname, "%s/snapdir_%03d/%s_%03d", All.OutputDir, RestartSnapNum, All.SnapshotFileBase, RestartSnapNum);
          else
            sprintf(fname, "%s%s_%03d", All.OutputDir, All.SnapshotFileBase, RestartSnapNum);

#ifdef SUBBOX_SNAPSHOTS
          if((RestartFlag == 3) || (RestartFlag == 5) || (RestartFlag == 10) || (RestartFlag == 15))
            {
              if((RestartFlag == 3) && (Argc < 5))
                {
                  mpi_printf("subbox missing: %d\n", Argc);
                  endrun();
                }
              if((RestartFlag == 5) && (Argc < 16))
                {
                  mpi_printf("subbox missing: %d\n", Argc);
                  endrun();
                }
              if((RestartFlag == 10) && (Argc != 8))
                {
                  mpi_printf("subbox missing: %d\n", Argc);
                  endrun();
                }
              else
                {
                  if(RestartFlag == 3)
                   {
                     All.SubboxNumber = atoi(Argv[4]);
                     All.NumFilesPerSnapshot = All.SubboxNumFilesPerSnapshot;
                     All.NumFilesWrittenInParallel = All.SubboxNumFilesWrittenInParallel;
                   }
                  if(RestartFlag == 5)
                    All.SubboxNumber = atoi(Argv[15]);
                  if(RestartFlag == 10)
                    All.SubboxNumber = atoi(Argv[7]);
                  if(RestartFlag == 15)
                    {
                      if(Argv[4] == 0)
                        {
                          terminate("<SubboxNum> is missing!\n");
                        }
                      All.SubboxNumber = atoi(Argv[4]);
                    }
                }
              if(All.SubboxNumFilesPerSnapshot > 1)
                sprintf(fname, "%s/snapdir_subbox%d_%03d/%s_subbox%d_%03d", All.OutputDir, All.SubboxNumber, RestartSnapNum, All.SnapshotFileBase, All.SubboxNumber, RestartSnapNum);
              else
                sprintf(fname, "%s/%s_subbox%d_%03d", All.OutputDir, All.SnapshotFileBase, All.SubboxNumber, RestartSnapNum);
            }
#endif
        }
      else
        strcpy(fname, All.InitCondFile);

#ifdef SHOCK_FINDER_POST_PROCESSING
      sprintf(All.OutputDir, "%s/", All.OutputDirShockFinder);
#endif

#ifdef POWERSPECTRUM_IN_POSTPROCESSING
      if(RestartFlag == 13)
        {
          int typemask = atoi(Argv[4]);

#ifdef POWERSPECTRUM_IN_POSTPROCESSING_ICS
          strcpy(All.InputFileName, All.InitCondFile);
#else
          strcpy(All.InputFileName, fname);
#endif

          long long n_type[NTYPES];
          for(int type=0; type<NTYPES; type++)
            {
              if(typemask & (0x1 << type))
                n_type[type] = 1;
              else
                n_type[type] = 0;
            }

          calculate_power_spectra(RestartSnapNum, n_type);
        
          endrun();
          
          return 0;
        }
      else
        {
          mpi_printf( "POWERSPECTRUM_IN_POSTPROCESSING only allows RestartFlag 13\n" );
          endrun();
          return 0;
        }
#endif

      /* now we can load the file */
      if((RestartFlag == 12) || (RestartFlag == 13))
        {
          int typemask = atoi(Argv[4]);
          read_ic(fname, typemask);
        }
      else /* readTypes=0x01: just load the gas particles. readTypes=LOAD_TYPES: load all particles */
        {
#ifdef READ_DM_AS_GAS
          read_ic(fname, ((RestartFlag == 4) || (RestartFlag == 5) || (RestartFlag == 10) || (RestartFlag == 14) || (RestartFlag == 16)) ? 0x02 : LOAD_TYPES);
#else
          read_ic(fname, ((RestartFlag == 4) || (RestartFlag == 5) || (RestartFlag == 10) || (RestartFlag == 14) || (RestartFlag == 15) || (RestartFlag == 16)) ? 0x01 : LOAD_TYPES);
#endif
        }

      /* If we are supposed to just convert the file, write and exit here. */
      if(RestartFlag == 6)
        {
#ifdef GFM_COOLING_METAL
          read_cooling_tables_current_time();
#endif
          /* important for proper functioning of FOF+SUBFIND */
          if(All.ComovingIntegrationOn) /* change to new velocity variable */
            {
              int i, j;
              for(i = 0; i < NumPart; i++)
                for(j = 0; j < 3; j++)
                  P[i].Vel[j] *= sqrt(All.Time) * All.Time;
            }
          set_softenings();
          All.TopNodeAllocFactor = 0.08;
          All.TreeAllocFactor = 0.7;
          All.NgbTreeAllocFactor = 0.7;

          sprintf(All.SnapshotFileBase, "%s_converted", All.SnapshotFileBase);
          mpi_printf("Start writing file %s\nRestartSnapNum %d\n", All.SnapshotFileBase, RestartSnapNum);
          savepositions(RestartSnapNum, 0);
          endrun();
        }

      /* init returns a status code, where a value of >=0 means that endrun() should be called. */
      int status = init();

      if(status >= 0)
        {
          if(status > 0)
            mpi_printf("init() returned with %d\n", status);

          endrun();
        }
    }

  begrun2();

  run();                        /* main simulation loop */

  endrun();                     /* clean up & finalize MPI */

  return 0;
}

/** \brief This function aborts the simulations.
 *
 * This method has to be called by all processes. It should be used only
 * if the simulation ends without a errors or a an error message is already printed.
 * Otherwise terminate() should be used instead.
 */
void endrun()
{
  mpi_printf("Code run for %f seconds!\n", timediff(StartOfRun, second()));
  mpi_printf("endrun called, calling MPI_Finalize()\nbye!\n\n");
  fflush(stdout);

#ifdef HAVE_HDF5
  /*The hdf5 library will sometimes register an atexit() handler that calls its error handler.
   * In arepo this is set to my_hdf_error_handler, which calls MPI_Abort.
   * Calling MPI_Abort after MPI_Finalize is not allowed.
   * Hence unset the HDF error handler here*/
  H5Eset_auto(NULL, NULL);
#endif

#ifdef CUDA
  cuda_finish();
#endif

#ifdef SNE_FEEDBACK
  sne_destroy();
#endif

#ifdef RUNNING_SAFETY_FILE
  if(All.Ti_Current < TIMEBASE) /* simulation has not reached the final time */
  {
    char running_fname[MAXLEN_PATH], running_done_fname[MAXLEN_PATH];
    sprintf(running_fname, "./running");
    sprintf(running_done_fname, "./running_done");
    rename(running_fname, running_done_fname);
    mpi_printf("moved ./running file to ./running_done, job can now restart.\n");
  }
  else
    mpi_printf("leaving ./running file in place since run is complete to prevent any restarts.\n");
#endif

  MPI_Finalize();
  exit(0);
}
