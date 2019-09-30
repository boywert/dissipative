/* PCC - 02.01.2014
	Contains the prototypes for the ITA-style sink particles 
*/

 void sink_particles(void);
 void init_sink_particles(int mode); 
 int get_all_sink_particle_info(int mode);
 double accrete_onto_sink_particles(void);
 int create_sink_particles(void);
 void add_sink_to_particle_structure(int idx, int candidate_index, int candidate_task);
 void dump_sink_particle_info(int success);
 int open_sink_particle_evolution_file(void);
 int open_old_sink_file_and_read_sink_data(void);

