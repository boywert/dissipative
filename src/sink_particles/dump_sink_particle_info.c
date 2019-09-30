#include "../allvars.h"
#include "../proto.h"
#include <sys/stat.h>
#include <sys/types.h>

/* PCC - 15.11.2014
        This function writes out the properties of all the sink particles that are present
        in the simulation at EACH timestep. At the moment, we're writing a bitstream. 
        Only Task 0 calls this function.
*/
void dump_sink_particle_info(int success)
{
  int i;

  if(ThisTask != 0)             /* only the root processors writes to the log files */
    return;

  /* Check that it's time to write a dump. If not, return...
 *   */
  printf("Checking if it is time to write sink properties to file %g %g %g \n", SinksLastEvolutionDumpTime, All.Time, SinkEvolutionDumpRateYears);
  if( ((All.Time - SinksLastEvolutionDumpTime) > SinkEvolutionDumpRateYears) || (success == 1))
    {
      SinksLastEvolutionDumpTime = All.Time;
      printf("SINK_PARTICLES: Writing the sink properities to the evolution file \n");

      fwrite(&All.Time, sizeof(double), 1, FdSinkPart); 
      fwrite(&NSinksAllTasks, sizeof(int), 1, FdSinkPart);
      fwrite(&SinkP[0], sizeof(struct global_sink_particle_data), NSinksAllTasks, FdSinkPart);
    }
}

/* PCC - 20.07.2015
        The function takes care of opening and (re)naming the files that store the
        time evolution of the sink particles. The idea here is that each 'run' (i.e
        restart) of the code will get a new file holding the sink particle info.
        Although this is a little clumsy, it's also the most flexible, and means that
        the code is backwards compatible with previous versions -- it doesn't care if
        new blocks are added in the future, just that the user reads the files sensibly!

        Empty files will be reused. In some (obscure) cases this might result in the
        file number being at odds with the temporal evolution, so it's probably best
        to read in all the files when using them in analysis. 

        The file FdSinkPart is opened below, but it is closed in close_logfiles(void)
        along with the other log files. This function is CALLED from open_logfiles(void)    

        Only Task 0 calls this function.
*/
int open_sink_particle_evolution_file(void )
{

  int i;
  int file_no;
  struct stat st;
  int max_file_no = 999;
  char sink_dir[MAXLEN_PATH];
  char filename[MAXLEN_PATH];
  char msg[100];

  /* Set the directory for the sink particle information
  */
  sprintf(sink_dir, "%s/sink_particle_info/", All.OutputDir);

  /* Make sure the sink directory is present 
  */
  if( stat(sink_dir, &st) == -1)
    mkdir(sink_dir,02755);
  else
    printf("SINK_PARTICLES: sink dump directory already exists\n");

  /* Now find the first available file name
  */
  file_no = 0;
  while (file_no < max_file_no)
    {
      sprintf(filename, "%ssink_particle_info/sink_particle_evolution_%03d", All.OutputDir, file_no);
      printf("SINK_PARTICLES: Checking for file %s \n", filename);
      if( (stat(filename, &st) == -1) || (stat(filename, &st) == 0 && st.st_size == 0))
        {
          /* set the sink particle file info and open the file */
         printf("SINK_PARTICLES: Sink particle evolution will be writen to file sink_particle_evolution_%03d \n", file_no);
         if(!(FdSinkPart = fopen(filename, "w")))
           {
             sprintf(msg, "error in opening file '%s'\n", filename);
             terminate(msg);
           }
         /*return with success!*/
         return(1);
        }
      printf("SINK_PARTICLES: Sink file ending %03d already exists \n", file_no);
      file_no++;
    }
  return(0);
}

/* This function is called if there are sinks present in the particle structure on startup.
   The function reads the sink_info files to get the last value of the SinkP that was written
   before the code stopped. This retrieves the data that's only stored on the SinkP array, such
   as the time the sink formed, etc, which gives us information on the sink particle histories.

   The function scans the exising sink dump file to find the first entry for which
     1) time <= start time in snapshot
     2) NSinks = NSinks in snapshot
   and terminates if nothing is found...

   FUTURE: might want to make this find the last dump before the snapshot was written
*/
int open_old_sink_file_and_read_sink_data(void )  
{
  int numsinks;
  int file_no;
  struct stat st;
  int max_file_no = 999;
  char sink_dir[MAXLEN_PATH];
  char filename[MAXLEN_PATH];
  char msg[100];
  double ScanTime;
  int ScanNSinks;  

  /* Find the number of sinks residing this snapshot's particle date
  */
  numsinks = get_all_sink_particle_info(0);
  mpi_printf("SINK PARTICLES READ -- numsinks %d NSinksAllTasks %d \n", numsinks, NSinksAllTasks);
  if(NSinksAllTasks < 1)
    return(0);

  /* Set the directory for the sink particle information
  */
  sprintf(sink_dir, "%s/sink_particle_info/", All.OutputDir);

  /* Make sure the sink directory is present 
  */
  if( stat(sink_dir, &st) == -1)
    { 
      sprintf("SINK PARTICLES READ -- problem! Sinks present but can't find sink directory %s\n", sink_dir);
      terminate(msg); 
    }
  else
    {
      mpi_printf("SINK_PARTICLES READ -- sink dump directory already exists \n");
      mpi_printf("Scaning through the files...");
    }

  file_no = 0;
  /* Loop around files
  */
  while (file_no < max_file_no)
    {
      sprintf(filename, "%ssink_particle_info/sink_particle_evolution_%03d", All.OutputDir, file_no);
      mpi_printf("SINK_PARTICLES: Checking for file %s \n", filename);
      if( stat(filename, &st) == 0  &&  st.st_size > 0)
        {
         /* open the file 
         */
         mpi_printf("SINK_PARTICLES: File present. scanning file...\n");
         if(!(FdSinkPart = fopen(filename, "r")))
           {
             sprintf(msg, "error in opening file '%s'\n", filename);
             terminate(msg);
           }
         /* Scan the file entries... 
         */
         while(1)
           {
             my_fread(&ScanTime, sizeof(double), 1, FdSinkPart);
             my_fread(&ScanNSinks, sizeof(int), 1, FdSinkPart);
             my_fread(&SinkP[0], sizeof(struct global_sink_particle_data), ScanNSinks, FdSinkPart);
             if(ScanTime <= All.Time  &&  ScanNSinks == NSinksAllTasks)
               {
                 /* return with success! 
                 */
                 mpi_printf("Found the sink record we want \n");
                 return(1);
               }
           }
        }
      mpi_printf("SINK_PARTICLES: Sink file ending %03d didn't help... Moving to next file... \n", file_no);
      file_no++;
    }

  mpi_printf("SINK_PARTICLES: No appropriate sink information from the dump files was found...\n");
  return(0);
}  
