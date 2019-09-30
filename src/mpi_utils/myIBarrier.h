#ifndef MYIBARRIER_H
#define MYIBARRIER_H

#ifdef MYIBARRIER
#define MPI_TAG_IBARRIER 0x666

struct sMyIBarrier {
  MPI_Comm comm;
  int rank;
  int nTasks;
  int nLevels;
  char *LevelDone;
};

void myIBarrier( MPI_Comm comm, struct sMyIBarrier * barrier );
void myIBarrierTest( struct sMyIBarrier * barrier, int *flag, MPI_Status *unused );
#endif

#endif /* MYIBARRIER_H */
