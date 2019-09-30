/* PCC - 19.12.2013
	This header holds the external structures for the sink 
        particle functions.
*/


extern int NSinksThisTask;
extern int NSinksAllTasks;
extern int NSinksAllTasksOld;
extern int NSinkBufferSize;
extern int SinksFormedSinceLastDomain;
extern double SinkCreationDensityCodeUnits;
extern double SinkAccretionRadius;
extern double SinkFormationRadiusSquared;
extern double SinkFormationRadius;
extern double SinkSofteningLength;
extern double SinkTargetAccretionMass;
extern double SinkTestsGiveupDensityCodeUnits;
extern double SinkEvolutionDumpRateYears;
extern double SinksLastEvolutionDumpTime;


/* This structure is the one that we'll use throughout Arepo */
extern struct global_sink_particle_data
{
  double Pos[3];
  double Vel[3];
  double Accel[3];
  double Mass;
  double FormationMass;
  double FormationTime;
  int ID;
  int HomeTask;
  int Index;
  int FormationOrder;
} *SinkP,
  *export_SinkP;

/*The following keeps track of which Task just formed a sink
*/
extern int LastTaskNewSink;


