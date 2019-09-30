#ifdef MYIBARRIER

#include <strings.h>
#include "myIBarrier.h"

void myIBarrier(MPI_Comm comm, struct sMyIBarrier *barrier)
{
  barrier->comm = comm;
  MPI_Comm_rank(comm, &barrier->rank);
  MPI_Comm_size(comm, &barrier->nTasks);

  barrier->nLevels = fls(barrier->rank - 1);
  barrier->LevelDone = mymalloc("myIBarrier", barrier->nLevels);
  memset(barrier->LevelDone, 0, barrier->nLevels);

  /* find messages we would expect from nonexisting tasks */
  for(level = 0; level < barrier->nLevels; level++)
    if((barrier->rank & (1 << level) == 0) && (barrier->rank + (1 << level) >= barrier->nTasks))
      barrier->LevelDone[level] = 1;

  /* find out if we have to send or wait */
  int level = 0;
  while(level < barrier->nLevels)
    {
      if(barrier->rank & (1 << level))
        {
          /* we need to send our result */
          int target = barrier->rank - (1 << level);
          int level = barrier->nLevels;
          MPI_Isend(&level, 1, MPI_INT, target, MPI_TAG_IBARRIER, barrier->comm);
          break;
        }
      else
        {
          /* check if there is something to recieve in which case we have to wait, otherwise go down one level */
          if(barrier->rank + (1 << level) < barrier->nTasks)
            {
              barrier->levelDone[level] = 1;
              break;
            }
          else
            level++;
        }
    }
}

void myIBarrierTest(struct sMyIBarrier *barrier, int *flag, MPI_Status * unused)
{
  flag = 0;

  int rflag;
  MPI_Status status;

  MPI_Iprobe(MPI_ANY_SOURCE, MPI_TAG_IBARRIER, barrier->comm, &rflag, &status);

  if(rflag)
    {
      int source = status.MPI_SOURCE;

      int level;
      MPI_Recv(&level, 1, MPI_INT, source, MPI_TAG_IBARRIER, barrier->comm, MPI_STATUS_IGNORE);

      if(source > barrier->rank)
        {
          /* we got another result, so lets check if we can send out further */
          while((level < barrier->nLevels) && barrier->LevelDone[level])
            level++;

          if(level == barrier->nLevels)
            {
              if(barrier->rank != 0)
                terminate("fail");
              /* ok, the barrier resolved, tell everyone */

              for(level = 0; level < barrier->nLevels; level++)
                {
                  if(barrier->rank & (1 << level) == 0)
                    {
                      int target = barrier->rank + (1 << level);
                      if(target < barrier->nTasks)
                        MPI_Isend(&level, 1, MPI_INT, target, MPI_TAG_IBARRIER, barrier->comm);
                    }
                  else
                    break;
                }

              flag = 1;
            }
          else
            {
              if(barrier->rank & (1 << level))
                {
                  /* we need to send our result */
                  int target = barrier->rank - (1 << level);
                  int level = barrier->nLevels;
                  MPI_Isend(&level, 1, MPI_INT, target, MPI_TAG_IBARRIER, barrier->comm);
                }
              else
                {
                  barrier->LevelDone[level] = 1;
                }
            }
        }
      else
        {
          for(; level < barrier->nLevels; level++)
            {
              if(barrier->rank & (1 << level) == 0)
                {
                  int target = barrier->rank + (1 << level);
                  if(target < barrier->nTasks)
                    MPI_Isend(&level, 1, MPI_INT, target, MPI_TAG_IBARRIER, barrier->comm);
                }
              else
                break;
            }

          flag = 1;
        }
    }
}

#endif
