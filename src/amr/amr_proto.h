/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/amr/amr_proto.h
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

#ifndef AMR_PROTO_H
#define AMR_PROTO_H

#ifdef AMR

/* general, amr.c */
void amr_init();
int amr_check_domain_decomposition();
void amr_update_oldMass();
void amr_reset_node_data(int no);
int amr_get_subnode(int node, int father);
void amr_accumulate_node_data(int no, int p);

/* amr_exchange.c */
void amr_export_node(int p, int dir, int ngb);
void mesh_setup_exchange();
int compare_list_export_data(const void *a, const void *b);
void amr_exchange_ghost_nodes(tessellation * T);
void amr_exchange_ghost_cells(tessellation * T);

void amr_update_primitive_variables_nodes();
void amr_exchange_primitive_variables();
void amr_calculate_gradients();
void amr_exchange_primitive_variables_and_gradients();



/* mesh generation, amr_mesh.c */
void create_mesh();

void amr_initmesh();

void free_mesh(void);

int amr_create_faces(tessellation * T);
void amr_create_face(int p1, int ngb, int direction, tessellation * T);
int amr_copy_dp(int i, tessellation * T);

void amr_set_cell_data(int cell, int father, int lists, tessellation * T);
void amr_set_node_data(int target, int parent, int subnode, int lists, tessellation * T);

void amr_free_mesh(tessellation * T);
void amr_allocate_mesh(tessellation * T);

void amr_treemodifylength(int delta_NgbMaxPart);
void amr_treerealloc(int delta_Nodes);


/* amr_ngb.c */
void amr_link_ngb_node(int new_node, int node, int subnode);
void amr_link_ngb_node_toplevel(int new_node, int node, int subnode);
void amr_link_ngb_particle(int particle, int node, int subnode);

void amr_link_ngb_particle_ghost(int particle, int node, int subnode);
void amr_link_ngb_node_ghost(int new_node, int th, int subnode);
void amr_link_ngb(int no);

void amr_link_toptree(int no);


void amr_update_nodes();

/* amr_search.c */
void amr_point_location(tessellation * T, double x, double y, double z, int* task, int* cell);

/* gradient calculation, amr_gradients.c */
void calculate_gradients(void);
void limit_gradient(double *d, double phi, double min_phi, double max_phi, MySingle *dphi);
void exchange_primitive_variables(void);
void exchange_primitive_variables_and_gradients(void);
void exchange_gradients(void);
void exchange_node_data(void);

int amr_should_this_node_be_split(int no);

void amr_validate_mesh(int node, int parent, int subnode);
void amr_validate_lists();


void amr_generate_gas_in_ics();

int amr_treefind_single_threads(MyDouble searchcenter[3], int target, int mode, int thread_id, int numnodes, int *firstnode);

/* image generation for 2D */
void do_special_dump(int num, int gradients_flag);

#endif /* AMR */ 

#endif /* AMR_PROTO_H */
