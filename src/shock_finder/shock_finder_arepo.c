/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/shock_finder/shock_finder_arepo.c
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

#include "../allvars.h"
#include "../proto.h"



#ifdef SHOCK_FINDER_AREPO

#define MOVEPRINTF(...)         //printf(__VA_ARGS__)


/*
 * shock finder for Arepo
*/


static void set_shock_zone();
static int is_in_shock_zone(int k);
static int calculate_shock_direction(int k, double *dir);
static void ray_action();
static void rays_ini();
static inline void insert_ray(double x, double y, double z, double dir_x, double dir_y, double dir_z, double surface_value, int current_cell, char preshock_dir);
static inline shock_ray *insert_ray_empty();
static inline void insert_export_ray(shock_ray * sr, unsigned int to_task);
static int ray_wrong_direction(MyFloat * shock_dir, shock_ray * ray, int this_task);
static void set_post_shock_quant(shock_ray * ray, int task, int original_index);
static void calculate_machnumber(shock_ray * ray, int this_task);
static void move_ray(shock_ray * ray);
static void copy_ray(shock_ray * origin, shock_ray * destination);
static void sf_exchange_rays(unsigned int *imported, unsigned int *exported);
static void rays_finish();
static void set_shock_surface();
static void set_machnumber();   //and pre-shock quantities
static void set_additional_info();

#ifdef SHOCK_FINDER_FILTER_INCONSISTENT
static int inconsistent_jump(shock_ray * ray);
#endif

static unsigned int *rays_send_count;   //!< rays_send_count[i]: number of elements which have to be send to task i
static unsigned int *rays_receive_count;        //!< rays_receive_count[i]: number of elements which have to be received from task i


//counters
static long long int machnumber_calculated = 0;
static long long int particles_in_shock_zone = 0;
static long long int particles_in_shock_surface = 0;
static unsigned int average_steps_to_preshock = 0;
static unsigned int wrong_jump_direction = 0;


//rays
static unsigned int MaxRays = 0;        //!< maximum number of elements in the array
static shock_ray *Rays = NULL;  //!< start of the array
static shock_ray *NewRay = NULL;        //!< first ray in array which is not in use

static shock_ray **RayHoles = NULL;     //!< save the holes in the Rays array
static unsigned int NofHoles = 0;       //!< number of elements in RayHoles array
static unsigned int MaxNofHoles = 0;    //!< maximum number of elements in NofHoles

static void reset_static_vars()
{
  machnumber_calculated = 0;
  particles_in_shock_zone = 0;
  particles_in_shock_surface = 0;
  average_steps_to_preshock = 0;
  wrong_jump_direction = 0;

  MaxRays = 0;
  Rays = NULL;
  NewRay = NULL;

  RayHoles = NULL;
  NofHoles = 0;
  MaxNofHoles = 0;
}

//rays to export
static struct
{
  unsigned int max_elements;    //!< maximum number of elements in elements
  shock_ray **elements;         //!< shock ray elements to send (array of adresses)
  unsigned int *to_task;        //!< send elements[i] to to_task[i]
  unsigned int counter;         //!< the first element in elements which is not yet used

} export_rays;

void shock_finder_arepo()
{
  if(NumGas != TimeBinsHydro.NActiveParticles)
    {
      mpi_printf("SHOCKS: Running shock finder for a local time step\n");
    }

  reset_static_vars();

  TIMER_START(CPU_SF_ZONE);

  mpi_printf("SHOCKS: Searching for shock zone\n");

  set_shock_zone();

  TIMER_STOP(CPU_SF_ZONE);


  TIMER_START(CPU_SF_SURFACE);

  exchange_shock_data();

#ifdef SHOCK_FINDER_ON_THE_FLY
  exchange_primitive_variables_and_gradients();
#endif

  rays_ini();

  mpi_printf("SHOCKS: Sending rays\n");

  ray_action();

  mpi_printf("SHOCKS: Setting shock surface\n");

  set_shock_surface();

  mpi_printf("SHOCKS: Assigning mach numbers\n");

  set_machnumber();

  rays_finish();

  TIMER_STOP(CPU_SF_SURFACE);

#ifdef SHOCK_FINDER_POST_PROCESSING
  TIMER_START(CPU_SF_MISC);

  set_additional_info();

  TIMER_STOP(CPU_SF_MISC);
#endif
}

static void set_shock_zone()
{
  int idx, i;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      //check whether the particle is in the shock zone and set flag correspondingly
      SDATA(i).ShockZone = is_in_shock_zone(i);
    }
#ifdef SHOCK_FINDER_POST_PROCESSING
  mpi_printf_each("\tTask %d: Particles in shock zone: %d\n", ThisTask, particles_in_shock_zone);
#endif
}

static int is_in_shock_zone(int k)
{
  VPRINTF("particle: (%f, %f, %f)\n", P[k].Pos[0], P[k].Pos[1], P[k].Pos[2]);

  if((P[k].Mass == 0) && (P[k].ID == 0))
    {
      return 0;
    }

#ifdef SHOCK_FINDER_FILTER_SFR
  //exclude cells on the effective equation of state (Springel & Hernquist 2003)
  if(SphP[k].Sfr > 0)
    {
      return 0;
    }
#endif


#ifdef SKIP_BORDER
#ifdef REFLECTIVE_X
  if(distance_to_border_x(k) < All.DistToBorder)
    {
      return 0;
    }
#endif

#ifdef REFLECTIVE_Y
  if(distance_to_border_y(k) < All.DistToBorder)
    {
      return 0;
    }
#endif

#ifdef REFLECTIVE_Z
  if(distance_to_border_z(k) < All.DistToBorder)
    {
      return 0;
    }
#endif

#endif



  //div v < 0 => maybe a shock zone (a factor not needed)
  if(divergence(k) >= 0)
    {
      return 0;
    }

  //grad(T)*grad(rho) > 0 => maybe a shock zone (a factors not needed)
#ifdef COSMIC_RAYS
  if(SDATA(k).CRpseudoTgrad[0] * SphP[k].Grad.drho[0] + SDATA(k).CRpseudoTgrad[1] * SphP[k].Grad.drho[1] + SDATA(k).CRpseudoTgrad[2] * SphP[k].Grad.drho[2] <= 0)
    {
      return 0;
    }
#else
  if(SDATA(k).Tgrad[0] * SphP[k].Grad.drho[0] + SDATA(k).Tgrad[1] * SphP[k].Grad.drho[1] + SDATA(k).Tgrad[2] * SphP[k].Grad.drho[2] <= 0)
    {
      return 0;
    }
#endif

  //variable for storing the shock direction
  double dir[3];


  //calculate the shock direction
  if(!calculate_shock_direction(k, dir))        //gradient has zero length => no shock zone
    {
      return 0;
    }

  VPRINTF("shock direction: (%f, %f, %f)\n", SphP[k].ShockDir[0], SphP[k].ShockDir[1], SphP[k].ShockDir[2]);


  //vector antiparallel to shock direction
  double asdir[3];

  asdir[0] = -1 * dir[0];
  asdir[1] = -1 * dir[1];
  asdir[2] = -1 * dir[2];

  //find neighbours
  int task;
  int edge;
  int boundary_reached;
  int original_cell;

  int preshockcell = find_neighbouring_cell_para(k, dir, &task, &original_cell, &edge, &boundary_reached);

  int this_task = (task == ThisTask);

#ifdef SKIP_BORDER
  int reflective_boundary_reached = 0;
  if(boundary_reached)
    {
#ifdef REFLECTIVE_X
      if(boundary_reached & REFL_X_FLAGS)
        {
          reflective_boundary_reached = 1;
        }
#endif
#ifdef REFLECTIVE_Y
      if(boundary_reached & REFL_Y_FLAGS)
        {
          reflective_boundary_reached = 1;
        }
#endif
#ifdef REFLECTIVE_Z
      if(boundary_reached & REFL_Z_FLAGS)
        {
          reflective_boundary_reached = 1;
        }
#endif

      if(reflective_boundary_reached)
        {
          terminate("ray reached the border, increase DistToBorder in param.txt");
        }
    }
#endif

#ifdef ZONE_JUMP_CR
  double pth1, pth2, pcr1, pcr2, rho1, rho2, gamma1, gamma2, e1, e2;
#endif

#ifdef ZONE_JUMP_V
  double v_preshock;
  double v_postshock;
  double c_preshock;
#endif

#ifdef ZONE_JUMP_P
  double p_preshock;
  double p_postshock;
#endif

#ifdef ZONE_JUMP_S
  double s_preshock;
  double s_postshock;
#endif

#ifdef ZONE_JUMP_T
  double T_preshock;
  double T_postshock;
#endif



  if(this_task)
    {
#ifdef ZONE_JUMP_CR
      pth1 = SphP[preshockcell].Pressure;
      pcr1 = SphP[preshockcell].CR_Pressure;
      rho1 = SphP[preshockcell].Density;
      e1 = SphP[preshockcell].Pressure / GAMMA_MINUS1 + SphP[preshockcell].CR_Pressure / (All.GammaCR - 1.);
      gamma1 = (pth1 + pcr1) / e1 + 1;
#endif
#ifdef ZONE_JUMP_V
      v_preshock = P[preshockcell].Vel[0] * asdir[0] + P[preshockcell].Vel[1] * asdir[1] + P[preshockcell].Vel[2] * asdir[2];

      myassert(SphP[preshockcell].Density != 0);
      c_preshock = sqrt(GAMMA * SphP[preshockcell].Pressure / SphP[preshockcell].Density);
#endif

#ifdef ZONE_JUMP_P
      p_preshock = SphP[preshockcell].Pressure;
#endif

#ifdef ZONE_JUMP_S
      s_preshock = SphP[preshockcell].Pressure / pow(SphP[preshockcell].Density, GAMMA);
#endif

#ifdef ZONE_JUMP_T
      T_preshock = SDATA(preshockcell).Temperature;
#endif

    }
  else                          //find neighbour function crossed the domain
    {
#ifdef ZONE_JUMP_CR
      pth1 = PrimExch[preshockcell].Pressure;
      pcr1 = PrimExch[preshockcell].CR_Pressure;
      rho1 = PrimExch[preshockcell].Density;
      e1 = PrimExch[preshockcell].Pressure / GAMMA_MINUS1 + PrimExch[preshockcell].CR_Pressure / (All.GammaCR - 1.);
      gamma1 = (pth1 + pcr1) / e1 + 1;
#endif
#ifdef ZONE_JUMP_V
      v_preshock = PrimExch[preshockcell].VelGas[0] * asdir[0] + PrimExch[preshockcell].VelGas[1] * asdir[1] + PrimExch[preshockcell].VelGas[2] * asdir[2];

      myassert(PrimExch[preshockcell].Density != 0);
      c_preshock = sqrt(GAMMA * PrimExch[preshockcell].Pressure / PrimExch[preshockcell].Density);
#endif

#ifdef ZONE_JUMP_P
      p_preshock = PrimExch[preshockcell].Pressure;
#endif

#ifdef ZONE_JUMP_S
      s_preshock = PrimExch[preshockcell].Pressure / pow(PrimExch[preshockcell].Density, GAMMA);
#endif

#ifdef ZONE_JUMP_T
      T_preshock = SDATAEXCH(preshockcell).Temperature;
#endif
    }

  VPRINTF("\tneighbour preshock: (%f, %f)\n", Mesh.DP[DC[edge].dp_index].x, Mesh.DP[DC[edge].dp_index].y);

  int postshockcell = find_neighbouring_cell_para(k, asdir, &task, &original_cell, &edge, &boundary_reached);

  this_task = (task == ThisTask);

#ifdef SKIP_BORDER
  reflective_boundary_reached = 0;
  if(boundary_reached)
    {
#ifdef REFLECTIVE_X
      if(boundary_reached & REFL_X_FLAGS)
        {
          reflective_boundary_reached = 1;
        }
#endif
#ifdef REFLECTIVE_Y
      if(boundary_reached & REFL_Y_FLAGS)
        {
          reflective_boundary_reached = 1;
        }
#endif
#ifdef REFLECTIVE_Z
      if(boundary_reached & REFL_Z_FLAGS)
        {
          reflective_boundary_reached = 1;
        }
#endif

      if(reflective_boundary_reached)
        {
          terminate("ray reached the border, increase DistToBorder in param.txt");
        }
    }
#endif

  if(this_task)
    {
#ifdef ZONE_JUMP_CR
      pth2 = SphP[postshockcell].Pressure;
      pcr2 = SphP[postshockcell].CR_Pressure;
      rho2 = SphP[postshockcell].Density;
      e2 = SphP[postshockcell].Pressure / GAMMA_MINUS1 + SphP[postshockcell].CR_Pressure / (All.GammaCR - 1.);
      gamma2 = (pth2 + pcr2) / e2 + 1;
#endif
#ifdef ZONE_JUMP_V
      v_postshock = P[postshockcell].Vel[0] * asdir[0] + P[postshockcell].Vel[1] * asdir[1] + P[postshockcell].Vel[2] * asdir[2];
#endif

#ifdef ZONE_JUMP_P
      p_postshock = SphP[postshockcell].Pressure;
#endif

#ifdef ZONE_JUMP_S
      s_postshock = SphP[postshockcell].Pressure / pow(SphP[postshockcell].Density, GAMMA);
#endif

#ifdef ZONE_JUMP_T
      T_postshock = SDATA(postshockcell).Temperature;
#endif
    }
  else                          //find neighbour function crossed the domain
    {
#ifdef ZONE_JUMP_CR
      pth2 = PrimExch[postshockcell].Pressure;
      pcr2 = PrimExch[postshockcell].CR_Pressure;
      rho2 = PrimExch[postshockcell].Density;
      e2 = PrimExch[postshockcell].Pressure / GAMMA_MINUS1 + PrimExch[postshockcell].CR_Pressure / (All.GammaCR - 1.);
      gamma2 = (pth2 + pcr2) / e2 + 1;
#endif
#ifdef ZONE_JUMP_V
      v_postshock = PrimExch[postshockcell].VelGas[0] * asdir[0] + PrimExch[postshockcell].VelGas[1] * asdir[1] + PrimExch[postshockcell].VelGas[2] * asdir[2];
#endif

#ifdef ZONE_JUMP_P
      p_postshock = PrimExch[postshockcell].Pressure;
#endif

#ifdef ZONE_JUMP_S
      s_postshock = PrimExch[postshockcell].Pressure / pow(PrimExch[postshockcell].Density, GAMMA);
#endif

#ifdef ZONE_JUMP_T
      T_postshock = SDATAEXCH(postshockcell).Temperature;
#endif
    }

  VPRINTF("\tneighbour postshock: (%f, %f)\n", Mesh.DP[DC[edge].dp_index].x, Mesh.DP[DC[edge].dp_index].y);

#ifdef ZONE_JUMP_CR

  if(rho2 < rho1 || (pth2 + pcr2) < (pth1 + pcr1))      //jump in wrong direction
    {
      return 0;
    }

  if(2 * fabs(gamma1 - gamma2) / fabs(gamma1 + gamma2) < 0.01)  //fallback to thermal criterion
    {
      double gamma_mean = 0.5 * (gamma1 + gamma2);

      double delta_log_p = log(pth2 + pcr2) - log(pth1 + pcr1);

      if(delta_log_p < 0)       //pressure jump in wrong direction
        {
          return 0;
        }

      if(delta_log_p < log((2 * gamma_mean * All.MachMin * All.MachMin - (gamma_mean - 1)) / (gamma_mean + 1))) //mach number too low
        {
          return 0;
        }
    }

  double Machnum = calculate_machnumber_cosmic_rays_zone(pth1, pth2, pcr1, pcr2, gamma1, gamma2);

  if(Machnum < All.MachMin)
    {
      return 0;
    }

#endif

#ifdef ZONE_JUMP_V
  double dv = (v_postshock - v_preshock) / All.cf_atime;        //physical peculiar velocity difference

  if(dv > 0)                    //gradient in wrong direction
    {
      return 0;
    }

  if(dv > 2. / (GAMMA + 1) * c_preshock * (1 - All.MachMin * All.MachMin) / All.MachMin)        //mach number too low
    {
      return 0;
    }
#endif

#ifdef ZONE_JUMP_P
  double delta_log_p = log(p_postshock) - log(p_preshock);

  if(delta_log_p < 0)           //pressure jump in wrong direction
    {
      return 0;
    }

  if(delta_log_p < log((2 * GAMMA * All.MachMin * All.MachMin - (GAMMA - 1)) / (GAMMA + 1)))    //mach number too low
    {
      return 0;
    }

#endif

#ifdef ZONE_JUMP_S
  double delta_log_s = log(s_postshock) - log(s_preshock);

  if(delta_log_s < 0)           //entropy jump in wrong direction
    {
      return 0;
    }

  if(delta_log_s < log(((2 * GAMMA * All.MachMin * All.MachMin - GAMMA + 1) / (GAMMA + 1)) * pow((GAMMA - 1 + 2 / (All.MachMin * All.MachMin)) / (GAMMA + 1), GAMMA)))  //mach number too low
    {
      return 0;
    }

#endif

#ifdef ZONE_JUMP_T
  double delta_log_T = log(T_postshock) - log(T_preshock);

  if(delta_log_T < 0)           //temperature jump in wrong direction
    {
      return 0;
    }

  if(delta_log_T < log(((2 * GAMMA * All.MachMin * All.MachMin) - (GAMMA - 1)) * ((GAMMA - 1) * All.MachMin * All.MachMin + 2) / ((GAMMA + 1) * (GAMMA + 1) * All.MachMin * All.MachMin)))      //mach number too low
    {
      return 0;
    }

#endif

  //cell is in the shock zone
  particles_in_shock_zone++;

  //set shock direction for particles in shock zone
  SDATA(k).ShockDir[0] = dir[0];
  SDATA(k).ShockDir[1] = dir[1];
  SDATA(k).ShockDir[2] = dir[2];

  return 1;
}

//calculate the shock direction
static int calculate_shock_direction(int k, double *dir)
{
#ifdef SHOCK_DIR_GRAD_P
  dir[0] = -(SphP[k].Grad.dpress[0]);
  dir[1] = -(SphP[k].Grad.dpress[1]);
  dir[2] = -(SphP[k].Grad.dpress[2]);
#endif

#ifdef SHOCK_DIR_GRAD_S
  dir[0] = -(SphP[k].Grad.dpress[0] / pow(SphP[k].Density, GAMMA) - GAMMA * SphP[k].Pressure * pow(SphP[k].Density, -GAMMA - 1) * SphP[k].Grad.drho[0]);
  dir[1] = -(SphP[k].Grad.dpress[1] / pow(SphP[k].Density, GAMMA) - GAMMA * SphP[k].Pressure * pow(SphP[k].Density, -GAMMA - 1) * SphP[k].Grad.drho[1]);
  dir[2] = -(SphP[k].Grad.dpress[2] / pow(SphP[k].Density, GAMMA) - GAMMA * SphP[k].Pressure * pow(SphP[k].Density, -GAMMA - 1) * SphP[k].Grad.drho[2]);
#endif

#ifdef SHOCK_DIR_GRAD_T
  //a factor not needed, gradient gets normalized
  dir[0] = -SDATA(k).Tgrad[0];
  dir[1] = -SDATA(k).Tgrad[1];
  dir[2] = -SDATA(k).Tgrad[2];
#endif

#ifdef SHOCK_DIR_GRAD_T_CR
  //a factor not needed, gradient gets normalized
  dir[0] = -SDATA(k).CRpseudoTgrad[0];
  dir[1] = -SDATA(k).CRpseudoTgrad[1];
  dir[2] = -SDATA(k).CRpseudoTgrad[2];
#endif

  double norm_square = dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2];

  if(norm_square == 0)          //gradient has zero length
    {
      return 0;
    }

  dir[0] /= sqrt(norm_square);
  dir[1] /= sqrt(norm_square);
  dir[2] /= sqrt(norm_square);

  return 1;
}



static void rays_ini()
{
  //memory for the rays
  MaxRays = (int) ceil((particles_in_shock_zone + 0.03 * NumGas) * All.RayMemFac);
  Rays = (shock_ray *) mymalloc("Rays", MaxRays * sizeof(shock_ray));
  NewRay = Rays;

  //memory for the holes
  MaxNofHoles = MaxRays;
  RayHoles = (shock_ray **) mymalloc("RayHoles", MaxNofHoles * sizeof(shock_ray *));
  NofHoles = 0;



  //memory for export rays
  export_rays.max_elements = MaxRays;
  export_rays.elements = (shock_ray **) mymalloc("export_rays.elements", export_rays.max_elements * sizeof(shock_ray *));
  export_rays.to_task = (unsigned int *) mymalloc("export_rays.to_task", export_rays.max_elements * sizeof(unsigned int));
  export_rays.counter = 0;

  if(Rays == NULL || export_rays.elements == NULL)
    {
      terminate("Error in rays_ini: memory allocation failed!\n");
    }

  //memory for the send and receive counts
  rays_send_count = (unsigned int *) mymalloc("rays_send_count", NTask * sizeof(unsigned int));
  rays_receive_count = (unsigned int *) mymalloc("rays_receive_count", NTask * sizeof(unsigned int));

}

static void rays_finish()
{
  myfree(rays_receive_count);
  myfree(rays_send_count);
  myfree(export_rays.to_task);
  myfree(export_rays.elements);
  myfree(RayHoles);
  myfree(Rays);
}

static inline void insert_ray(double x, double y, double z, double dir_x, double dir_y, double dir_z, double surface_value, int current_cell, char preshock_dir)
{
  if(NewRay - Rays == MaxRays)
    {
      printf("Error: Maximum number of rays: %d\n", MaxRays);
      terminate("Error in insert_ray: Not enough memory for rays. Increase All.RayMemFac!\n");
    }

  NewRay->startpos[0] = x;
  NewRay->startpos[1] = y;
  NewRay->startpos[2] = z;
  NewRay->pos[0] = x;
  NewRay->pos[1] = y;
  NewRay->pos[2] = z;
  NewRay->dir[0] = dir_x;
  NewRay->dir[1] = dir_y;
  NewRay->dir[2] = dir_z;
  NewRay->surface_value = surface_value;
  NewRay->current_cell = current_cell;
  NewRay->preshock_dir = preshock_dir;
  NewRay->previous_cell = current_cell;
  NewRay->original_cell = current_cell;
  NewRay->machnumber = 0;
  NewRay->original_cell_task = ThisTask;
  NewRay->previous_cell_task = ThisTask;
  NewRay->steps = 0;
  NewRay->quantities_set = 0;
  NewRay->rho_pre_shock = 0;
  NewRay->rho_post_shock = 0;
  NewRay->p_post_shock = 0;
  NewRay->v_pre_shock[0] = 0;
  NewRay->v_pre_shock[1] = 0;
  NewRay->v_pre_shock[2] = 0;
  NewRay->v_post_shock[0] = 0;
  NewRay->v_post_shock[1] = 0;
  NewRay->v_post_shock[2] = 0;

#ifdef COSMIC_RAYS
  NewRay->pth_post_shock = 0;
  NewRay->pcr_post_shock = 0;
  NewRay->pth_pre_shock = 0;
  NewRay->pcr_pre_shock = 0;
  NewRay->eth_post_shock = 0;
  NewRay->ecr_post_shock = 0;
  NewRay->gte = 0;
#ifdef COSMIC_RAYS_SHOCK_ACCELERATION
#ifdef COSMIC_RAYS_MAGNETIC_OBLIQUITY
  NewRay->magnetic_obliquity_pre_shock = 0;
#endif
  NewRay->eth_tot_pre_shock = 0;
  NewRay->post_shock_eth_tot_sum = injection_weight(current_cell);
  NewRay->cell_count = 0;
  NewRay->post_shock_cells[0] = -1;
  NewRay->post_shock_cells[1] = -1;
  NewRay->post_shock_cells[2] = -1;
  NewRay->post_shock_cells_tasks[0] = -1;
  NewRay->post_shock_cells_tasks[1] = -1;
  NewRay->post_shock_cells_tasks[2] = -1;
#endif
#else
  NewRay->post_shock_quant = 0;
  NewRay->c_pre_shock = 0;
  NewRay->p_pre_shock = 0;
  NewRay->t_pre_shock = 0;
#endif

  NewRay++;
}


/*!
 * 	\brief Insert a ray to the ray-array of the local task
 *
 *	\return The new ray
 */
static inline shock_ray *insert_ray_empty()
{
  if(NofHoles > 0)
    {
      NofHoles--;
      return RayHoles[NofHoles];
    }

  if(NewRay - Rays == MaxRays)
    {
      printf("Error: Maximum number of rays: %d\n", MaxRays);
      terminate("Error in insert_ray: Not enough memory for rays. Increase All.RayMemFac\n");
    }

  return NewRay++;
}

/*!
 *  \brief Free a ray, indicate that its position in the array can be reused
 */
static inline void free_ray(shock_ray * ray)
{
  RayHoles[NofHoles] = ray;
  NofHoles++;
  assert(NofHoles != MaxNofHoles);
}

/*!
 * 	\brief Insert a ray to the static export_rays struct
 *
 *  \param sr A pointer to the shock ray which has to be exported
 *  \param to_task The task to which the ray has to be exported
 */
static inline void insert_export_ray(shock_ray * sr, unsigned int to_task)
{
  if(export_rays.counter == export_rays.max_elements)
    {
      terminate("Error in insert_export_ray: Not enough memory for rays!\n");
    }

  export_rays.elements[export_rays.counter] = sr;
  export_rays.to_task[export_rays.counter] = to_task;

  export_rays.counter++;
}

//check whether a ray moves in disagreement with the shock direction (helper for move_ray)
static int ray_wrong_direction(MyFloat * shock_dir, shock_ray * ray, int this_task)
{

#ifdef DISABLE_DIRECTION_CHECK
  return 0;
#endif

  if(this_task)
    {
      if(!SDATA(ray->current_cell).ShockZone)   //outside of the shock zone the direction is irrelevant
        {
          return 0;
        }
    }
  else
    {
      if(!SDATAEXCH(ray->current_cell).ShockZone)       //outside of the shock zone the direction is irrelevant
        {
          return 0;
        }
    }

  double projection_square = shock_dir[0] * ray->dir[0] + shock_dir[1] * ray->dir[1] + shock_dir[2] * ray->dir[2];

  if(projection_square > 0 && ray->preshock_dir)        //correct direction
    {
      return 0;
    }
  else if(projection_square < 0 && !ray->preshock_dir)  //correct direction
    {
      return 0;
    }
  else                          //wrong direction
    {
#ifdef SHOCK_FINDER_VERBOSE
      mpi_printf("RAY: WRONG DIRECTION!\n");
#endif
      return 1;
    }
}


//check for high overdensities whether the Mach number jumps are consistent
#ifdef SHOCK_FINDER_FILTER_INCONSISTENT
static int inconsistent_jump(shock_ray * ray)
{

  //check 1d divergence direction
  //

  double div1d =   (ray->dir[0]*ray->v_pre_shock[0] + ray->dir[1]*ray->v_pre_shock[1] + ray->dir[2]*ray->v_pre_shock[2])
                      - (ray->dir[0]*ray->v_post_shock[0] + ray->dir[1]*ray->v_post_shock[1] + ray->dir[2]*ray->v_post_shock[2]);

  if(div1d >= 0)
    {
      return 1; //inconsistent
    }

  //check pressure and density jump directions
  //
 
#ifdef COSMIC_RAYS
  if(ray->pth_pre_shock + ray->pcr_pre_shock > ray->pth_post_shock + ray->pcr_post_shock || ray->rho_pre_shock>ray->rho_post_shock)
    {
      return 1; //inconsistent
    }
#else
  if(ray->p_pre_shock>ray->p_post_shock || ray->rho_pre_shock>ray->rho_post_shock)
    {
      return 1; //inconsistent
    }
#endif

  //for preshock baryonic overdensities above 1000: check consistency of jumps
  //
  double tolerance = 2.0;

  if( All.ComovingIntegrationOn && (All.OmegaBaryon != 0) )
    {
      double overdensity_thresh = 1000;
      double overdensity = ray->rho_pre_shock/(All.OmegaBaryon*3*All.Hubble*All.Hubble/(8*M_PI*All.G))-1; // a factors and h factors cancel

      if( overdensity < overdensity_thresh )
        {
          return 0; // consistent, check only for high density regions
        }
    }

#ifdef COSMIC_RAYS
  double M_e = calculate_machnumber_cosmic_rays_from_energy_equation(ray->pth_pre_shock, ray->pth_post_shock, ray->pcr_pre_shock, ray->pcr_post_shock, ray->rho_pre_shock, ray->rho_post_shock);

  if(ray->machnumber>M_e*tolerance || ray->machnumber<M_e/tolerance)
    {
      return 1; //inconsistent
    }

  return 0; //consistent
#else

  //check for inconsistent pressure jumps
  //

  double M_p = calculate_machnumber_p_jump(ray->p_pre_shock, ray->p_post_shock);

  if(ray->machnumber>M_p*tolerance || ray->machnumber<M_p/tolerance)
    {
      return 1; //inconsistent
    }

  //check for inconsistent density jumps
  //

  double M_rho;

  if(ray->rho_post_shock/ray->rho_pre_shock>=(GAMMA+1.)/(GAMMA-1.)) // unphysical, set Mach number to zero
    {
      M_rho = 0;
    }
  else
    {
      M_rho = calculate_machnumber_rho_jump(ray->rho_pre_shock, ray->rho_post_shock);
    }

  if(!((M_rho>=1 && M_rho<=3) || (ray->machnumber>=1 && ray->machnumber <=3))) //perform the check only if density or temperature Mach number is in range [1,3] (density jump not sensitive for high Mach numbers)
      {
        return 0; //consistent
      }


    if(ray->machnumber>M_rho*tolerance || ray->machnumber<M_rho/tolerance)
    {
      return 1; //inconsistent
    }

  return 0; //consistent
#endif
}
#endif

#ifdef SHOCK_FINDER_FILTER_SFR
//check whether gas is in the effective equation of state range
static int in_effective_eos_range(double rho, double T)
{
  double rho_phys = rho * All.cf_a3inv; //modulo little h, see sfr_eEOS.c: if(dens * All.cf_a3inv >= All.PhysDensThresh) ...


  if(rho_phys<All.PhysDensThresh)
    {
      return 0;
    }

  double T_thresh=1e5*pow(rho_phys/All.PhysDensThresh,0.2);

  if(T>T_thresh)
    {
      return 0;
    }

  return 1;
}
#endif


//set the post shock quantity (helper for move ray)
static void set_post_shock_quant(shock_ray * ray, int task, int original_index)
{
  int this_task = (task == ThisTask);

  if(ray->quantities_set == 1)
    {
      return;
    }
  else
    {
      ray->quantities_set = 1;
    }

  if(this_task)
    {
#ifdef SHOCK_JUMP_CR
#ifdef COSMIC_RAYS_SHOCK_ACCELERATION
      if(ray->steps <= 3)
        {
          assert(ray->post_shock_cells[ray->steps - 1] == -1);
          ray->post_shock_eth_tot_sum += injection_weight(ray->current_cell);
          ray->cell_count++;
          ray->post_shock_cells[ray->steps - 1] = ray->current_cell;
          ray->post_shock_cells_tasks[ray->steps - 1] = ThisTask;
        }
#endif
      ray->pth_post_shock = SphP[ray->current_cell].Pressure;
      ray->pcr_post_shock = SphP[ray->current_cell].CR_Pressure;
      ray->eth_post_shock = SphP[ray->current_cell].Pressure / GAMMA_MINUS1;
      ray->ecr_post_shock = SphP[ray->current_cell].CR_Pressure / (All.GammaCR - 1.);
#endif

      ray->rho_post_shock = SphP[ray->current_cell].Density;
      ray->p_post_shock = SphP[ray->current_cell].Pressure;
      ray->v_post_shock[0] = P[ray->current_cell].Vel[0];
      ray->v_post_shock[1] = P[ray->current_cell].Vel[1];
      ray->v_post_shock[2] = P[ray->current_cell].Vel[2];

#ifdef SHOCK_JUMP_V
      ray->post_shock_quant = P[ray->current_cell].Vel[0] * ray->dir[0] + P[ray->current_cell].Vel[1] * ray->dir[1] + P[ray->current_cell].Vel[2] * ray->dir[2];
#endif

#ifdef SHOCK_JUMP_P
      ray->post_shock_quant = SphP[ray->current_cell].Pressure;
#endif

#ifdef SHOCK_JUMP_T
      ray->post_shock_quant = SDATA(ray->current_cell).Temperature;
#endif

#ifdef SHOCK_JUMP_S
      ray->post_shock_quant = SphP[ray->current_cell].Pressure / pow(SphP[ray->current_cell].Density, GAMMA);
#endif
    }
  else                          //current cell is on an other task
    {
#ifdef SHOCK_JUMP_CR
#ifdef COSMIC_RAYS_SHOCK_ACCELERATION
      if(ray->steps <= 3)
        {
          assert(ray->post_shock_cells[ray->steps - 1] == -1);
          ray->post_shock_eth_tot_sum += injection_weight_exch(ray->current_cell);
          ray->cell_count++;
          ray->post_shock_cells[ray->steps - 1] = original_index;
          ray->post_shock_cells_tasks[ray->steps - 1] = task;
        }
#endif
      ray->pth_post_shock = PrimExch[ray->current_cell].Pressure;
      ray->pcr_post_shock = PrimExch[ray->current_cell].CR_Pressure;
      ray->eth_post_shock = PrimExch[ray->current_cell].Pressure / GAMMA_MINUS1;
      ray->ecr_post_shock = PrimExch[ray->current_cell].CR_Pressure / (All.GammaCR - 1);
#endif

      ray->rho_post_shock = PrimExch[ray->current_cell].Density;
      ray->p_post_shock = PrimExch[ray->current_cell].Pressure;
      ray->v_post_shock[0] = PrimExch[ray->current_cell].VelGas[0];
      ray->v_post_shock[1] = PrimExch[ray->current_cell].VelGas[1];
      ray->v_post_shock[2] = PrimExch[ray->current_cell].VelGas[2];

#ifdef SHOCK_JUMP_V
      ray->post_shock_quant = PrimExch[ray->current_cell].VelGas[0] * ray->dir[0] + PrimExch[ray->current_cell].VelGas[1] * ray->dir[1] + PrimExch[ray->current_cell].VelGas[2] * ray->dir[2];
#endif

#ifdef SHOCK_JUMP_P
      ray->post_shock_quant = PrimExch[ray->current_cell].Pressure;
#endif

#ifdef SHOCK_JUMP_T
      ray->post_shock_quant = SDATAEXCH(ray->current_cell).Temperature;
#endif

#ifdef SHOCK_JUMP_S
      ray->post_shock_quant = PrimExch[ray->current_cell].Pressure / pow(PrimExch[ray->current_cell].Density, GAMMA);
#endif
    }
}

//calculate the mach number (helper for move ray)
static void calculate_machnumber(shock_ray * ray, int this_task)
{
  assert(ray->quantities_set >= 1);     //the post-shock quantity should already be set

  if(ray->quantities_set == 2)
    {
      return;
    }
  else
    {
      ray->quantities_set = 2;
    }

#ifdef SHOCK_JUMP_CR

  //pre shock variables
  double pth1, pcr1, rho1, ecr1, eth1;

  //post shock variables
  double pth2 = ray->pth_post_shock;
  double pcr2 = ray->pcr_post_shock;
  double rho2 = ray->rho_post_shock;
  double ecr2 = ray->ecr_post_shock;
  double eth2 = ray->eth_post_shock;
#endif

  if(this_task)
    {
#ifdef SHOCK_JUMP_CR
      //pre shock quantities
      pth1 = SphP[ray->current_cell].Pressure;
      pcr1 = SphP[ray->current_cell].CR_Pressure;
      rho1 = SphP[ray->current_cell].Density;
      ecr1 = SphP[ray->current_cell].CR_Pressure / (All.GammaCR - 1.);
      eth1 = SphP[ray->current_cell].Pressure / GAMMA_MINUS1;
      ray->pth_pre_shock = SphP[ray->current_cell].Pressure;
      ray->pcr_pre_shock = SphP[ray->current_cell].CR_Pressure;
#ifdef COSMIC_RAYS_SHOCK_ACCELERATION
#ifdef COSMIC_RAYS_MAGNETIC_OBLIQUITY
      ray->magnetic_obliquity_pre_shock = magnetic_shock_obliquity(ray->dir, SphP[ray->current_cell].B);
#endif
      ray->eth_tot_pre_shock = injection_weight(ray->current_cell);
#endif
#else
      //set pre shock qunatities
      ray->c_pre_shock = sqrt(GAMMA * SphP[ray->current_cell].Pressure / SphP[ray->current_cell].Density);      //a factors cancel
      ray->p_pre_shock = SphP[ray->current_cell].Pressure;
      ray->t_pre_shock = SDATA(ray->current_cell).Temperature;
#endif
      ray->rho_pre_shock = SphP[ray->current_cell].Density;     //comoving density
      ray->v_pre_shock[0] = P[ray->current_cell].Vel[0];
      ray->v_pre_shock[1] = P[ray->current_cell].Vel[1];
      ray->v_pre_shock[2] = P[ray->current_cell].Vel[2];

      //set pre shock sound speed in case of temperature floor
#ifdef SHOCK_JUMP_T_FLOOR
      if(SDATA(ray->current_cell).Temperature < SHOCK_JUMP_T_FLOOR)
        {
#ifdef COOLING
          double meanweight = 4. / (3 * HYDROGEN_MASSFRAC + 1 + 4 * HYDROGEN_MASSFRAC * SphP[ray->current_cell].Ne);
#else
          double meanweight = 0.5882352941176471;       //fully ionized
#endif
          ray->c_pre_shock = sqrt(GAMMA * BOLTZMANN * SHOCK_JUMP_T_FLOOR * All.UnitMass_in_g / (All.UnitEnergy_in_cgs * meanweight * PROTONMASS));
        }
#endif

#ifdef SHOCK_JUMP_V

      double v_preshock = P[ray->current_cell].Vel[0] * (-ray->dir[0]) + P[ray->current_cell].Vel[1] * (-ray->dir[1]) + P[ray->current_cell].Vel[2] * (-ray->dir[2]);

      myassert(SphP[ray->current_cell].Density != 0);
      double c_preshock = sqrt(GAMMA * SphP[ray->current_cell].Pressure / SphP[ray->current_cell].Density);

      double dv = (ray->post_shock_quant - v_preshock) / All.cf_atime;

      if(dv > 0)                //gradient in wrong direction
        {
          ray->machnumber = 0;
        }
      else
        {
          ray->machnumber = calculate_machnumber_vdiff(dv, c_preshock);
        }
#endif

#ifdef SHOCK_JUMP_P
      double p_preshock = SphP[ray->current_cell].Pressure;

      if(ray->post_shock_quant - p_preshock < 0)        //jump in wrong direction
        {
          ray->machnumber = 0;
        }
      else
        {
          ray->machnumber = calculate_machnumber_p_jump(p_preshock, ray->post_shock_quant);
        }
#endif

#ifdef SHOCK_JUMP_T

      double T_preshock = SDATA(ray->current_cell).Temperature;

#ifdef SHOCK_JUMP_T_FLOOR
      if(T_preshock < SHOCK_JUMP_T_FLOOR)
        {
          T_preshock = SHOCK_JUMP_T_FLOOR;
          ray->t_pre_shock = SHOCK_JUMP_T_FLOOR;
        }

      if(ray->post_shock_quant < SHOCK_JUMP_T_FLOOR)
        {
          ray->post_shock_quant = SHOCK_JUMP_T_FLOOR;
        }
#endif

      if(ray->post_shock_quant - T_preshock <= 0)        //jump in wrong direction
        {
          ray->machnumber = 0;
        }
#ifdef RESET_WRONG_RHO_JUMPS
      else if(ray->rho_post_shock - ray->rho_pre_shock <= 0)
        {
          ray->machnumber = 0;
        }
#endif
#ifdef RESET_WRONG_P_JUMPS
      else if(ray->p_post_shock - ray->p_pre_shock <= 0)
        {
          ray->machnumber = 0;
        }
#endif
      else
        {
          ray->machnumber = calculate_machnumber_T_jump(T_preshock, ray->post_shock_quant);
        }
#endif

#ifdef SHOCK_JUMP_S

      myassert(SphP[ray->current_cell].Density != 0);
      double s_preshock = SphP[ray->current_cell].Pressure / pow(SphP[ray->current_cell].Density, GAMMA);

      if(ray->post_shock_quant - s_preshock < 0)        //jump in wrong direction
        {
          ray->machnumber = 0;
        }
      else
        {
          ray->machnumber = calculate_machnumber_S_jump(s_preshock, ray->post_shock_quant, All.MachMin);
        }
#endif
    }
  else                          //current cell is on an other task
    {
#ifdef SHOCK_JUMP_CR
      //pre shock quantities
      pth1 = PrimExch[ray->current_cell].Pressure;
      pcr1 = PrimExch[ray->current_cell].CR_Pressure;
      rho1 = PrimExch[ray->current_cell].Density;
      ecr1 = PrimExch[ray->current_cell].CR_Pressure / (All.GammaCR - 1.);
      eth1 = PrimExch[ray->current_cell].Pressure / GAMMA_MINUS1;
      ray->pth_pre_shock = PrimExch[ray->current_cell].Pressure;
      ray->pcr_pre_shock = PrimExch[ray->current_cell].CR_Pressure;
#ifdef COSMIC_RAYS_SHOCK_ACCELERATION
#ifdef COSMIC_RAYS_MAGNETIC_OBLIQUITY
      ray->magnetic_obliquity_pre_shock = magnetic_shock_obliquity(ray->dir, PrimExch[ray->current_cell].B);
#endif
      ray->eth_tot_pre_shock = injection_weight_exch(ray->current_cell);
#endif
#else
      //set pre shock quantities
      ray->c_pre_shock = sqrt(GAMMA * PrimExch[ray->current_cell].Pressure / PrimExch[ray->current_cell].Density);      //a factors cancel
      ray->p_pre_shock = PrimExch[ray->current_cell].Pressure;
      ray->t_pre_shock = SDATAEXCH(ray->current_cell).Temperature;
#endif
      ray->rho_pre_shock = PrimExch[ray->current_cell].Density; //comoving density
      ray->v_pre_shock[0] = PrimExch[ray->current_cell].VelGas[0];
      ray->v_pre_shock[1] = PrimExch[ray->current_cell].VelGas[1];
      ray->v_pre_shock[2] = PrimExch[ray->current_cell].VelGas[2];

      //set pre shock sound speed in case of temperature floor
#ifdef SHOCK_JUMP_T_FLOOR
      if(SDATAEXCH(ray->current_cell).Temperature < SHOCK_JUMP_T_FLOOR)
        {
#ifdef COOLING
          double meanweight = 4. / (3 * HYDROGEN_MASSFRAC + 1 + 4 * HYDROGEN_MASSFRAC * PrimExch[ray->current_cell].Ne);
#else
          double meanweight = 0.5882352941176471;       //fully ionized
#endif
          ray->c_pre_shock = sqrt(GAMMA * BOLTZMANN * SHOCK_JUMP_T_FLOOR * All.UnitMass_in_g / (All.UnitEnergy_in_cgs * meanweight * PROTONMASS));
        }
#endif

#ifdef SHOCK_JUMP_V
      double v_preshock = PrimExch[ray->current_cell].VelGas[0] * (-ray->dir[0]) + PrimExch[ray->current_cell].VelGas[1] * (-ray->dir[1]) + PrimExch[ray->current_cell].VelGas[2] * (-ray->dir[2]);

      myassert(PrimExch[ray->current_cell].Density != 0);
      double c_preshock = sqrt(GAMMA * PrimExch[ray->current_cell].Pressure / PrimExch[ray->current_cell].Density);

      double dv = (ray->post_shock_quant - v_preshock) / All.cf_atime;

      if(dv > 0)                //gradient in wrong direction
        {
          ray->machnumber = 0;
        }
      else
        {
          ray->machnumber = calculate_machnumber_vdiff(dv, c_preshock);
        }
#endif

#ifdef SHOCK_JUMP_P
      double p_preshock = PrimExch[ray->current_cell].Pressure;

      if(ray->post_shock_quant - p_preshock < 0)        //jump in wrong direction
        {
          ray->machnumber = 0;
        }
      else
        {
          ray->machnumber = calculate_machnumber_p_jump(p_preshock, ray->post_shock_quant);
        }
#endif

#ifdef SHOCK_JUMP_T

      double T_preshock = SDATAEXCH(ray->current_cell).Temperature;

#ifdef SHOCK_JUMP_T_FLOOR
      if(T_preshock < SHOCK_JUMP_T_FLOOR)
        {
          T_preshock = SHOCK_JUMP_T_FLOOR;
          ray->t_pre_shock = SHOCK_JUMP_T_FLOOR;
        }

      if(ray->post_shock_quant < SHOCK_JUMP_T_FLOOR)
        {
          ray->post_shock_quant = SHOCK_JUMP_T_FLOOR;
        }
#endif

      if(ray->post_shock_quant - T_preshock <= 0)        //jump in wrong direction
        {
          ray->machnumber = 0;
        }
#ifdef RESET_WRONG_RHO_JUMPS
      else if(ray->rho_post_shock - ray->rho_pre_shock <= 0)
        {
          ray->machnumber = 0;
        }
#endif
#ifdef RESET_WRONG_P_JUMPS
      else if(ray->p_post_shock - ray->p_pre_shock <= 0)
        {
          ray->machnumber = 0;
        }
#endif
      else
        {
          ray->machnumber = calculate_machnumber_T_jump(T_preshock, ray->post_shock_quant);
        }
#endif

#ifdef SHOCK_JUMP_S

      myassert(PrimExch[ray->current_cell].Density != 0);
      double s_preshock = PrimExch[ray->current_cell].Pressure / pow(PrimExch[ray->current_cell].Density, GAMMA);

      if(ray->post_shock_quant - s_preshock < 0)        //jump in wrong direction
        {
          ray->machnumber = 0;
        }
      else
        {
          ray->machnumber = calculate_machnumber_S_jump(s_preshock, ray->post_shock_quant, All.MachMin);
        }
#endif

    }

#ifdef SHOCK_JUMP_CR
  if(rho2 < rho1 || (pth2 + pcr2) < (pth1 + pcr1))      //jump in wrong direction
    {
      ray->machnumber = 0;
      ray->gte = 0;
    }
  else
    {
#ifdef MACHNUMBER_CR_2
      double gamma1 = (pth1 + pcr1) / (eth1 + ecr1) + 1;
      double gamma2 = (pth2 + pcr2) / (eth2 + ecr2) + 1;

      if(2 * fabs(gamma1 - gamma2) / fabs(gamma1 + gamma2) < 0.01)      //fallback to thermal criterion
        {
          myassert(pth1 + pcr1 > 0);

          double gamma_mean = 0.5 * (gamma1 + gamma2);
          ray->machnumber = sqrt((gamma_mean + 1) * (pth2 + pcr2) / (2 * gamma_mean * (pth1 + pcr1)) + (gamma_mean - 1) / (2 * gamma_mean));
        }
      else
        {
          ray->machnumber = calculate_machnumber_cosmic_rays_zone(pth1, pth2, pcr1, pcr2, gamma1, gamma2);
        }
#else
      ray->machnumber = calculate_machnumber_cosmic_rays(pth1, pth2, pcr1, pcr2, rho1, rho2);
#endif
      double c1 = sound_speed_cosmic_rays(pcr1, pth1, rho1);
      ray->gte = generated_thermal_energy_flux_cosmic_rays(eth1, eth2, ecr1, ecr2, rho1, rho2, ray->machnumber, c1);

      assert(ray->machnumber > 0);

      if(ray->gte < 0)          //reset incorrect jumps
        {
          ray->machnumber = 0;
          ray->gte = 0;
        }
    }
#endif

  if(ray->machnumber != 0)
    {
      machnumber_calculated++;
      average_steps_to_preshock += ray->steps;
    }
}

/*!
 *  \brief Moves the ray to the next Voronoi cell
 *
 * 	\param ray The ray to move
 */
static void move_ray(shock_ray * ray)
{
  int old_cell = ray->current_cell;
  int task;
  int original_index;
  int edge;
  int boundary_reached;

  ray->current_cell = find_neighbouring_cell_pos_para(ray->current_cell, ray->pos, ray->dir, ray->previous_cell, ray->previous_cell_task, ray->pos, &task, &original_index, &edge, &boundary_reached);

  int this_task = (task == ThisTask);

#ifdef SHOCK_FINDER_FILTER_SFR
  //exclude cells on the effective equation of state (Springel & Hernquist 2003)
  if(this_task)
    {
      if(SphP[ray->current_cell].Sfr > 0)
        {
          ray->current_cell = -3;       //deactivate ray
          ray->machnumber = 0;

          free_ray(ray);
          return;
        }
    }
  else
    {
      if(PrimExch[ray->current_cell].Sfr > 0)
        {
          ray->current_cell = -3;       //deactivate ray
          ray->machnumber = 0;

          free_ray(ray);
          return;
        }
    }
#endif

  double xtmp;
  double dx = NEAREST_X( ray->startpos[0] - ray->pos[0] );
#if NUMDIMS > 1
  double ytmp;
  double dy = NEAREST_Y( ray->startpos[1] - ray->pos[1] );
#endif
#if NUMDIMS > 2
  double ztmp;
  double dz = NEAREST_Z( ray->startpos[2] - ray->pos[2] );
#endif

#if NUMDIMS ==1
  if(fabs(dx) > 0.25 * boxSize_X)
    {
      ray->current_cell = -3;			//deactivate ray
      ray->machnumber = 0;
      free_ray(ray);
      return;
    }
#else /* NUMDIMS > 1 */
#if NUMDIMS == 2
  if(fabs(dx) > 0.25 * boxSize_X || fabs(dy) > 0.25 * boxSize_Y)
    {
      ray->current_cell = -3;			//deactivate ray
      ray->machnumber = 0;
      free_ray(ray);
      return;
    }
#else /* NUMDIMS > 2 */
  if(fabs(dx) > 0.25 * boxSize_X || fabs(dy) > 0.25 * boxSize_Y || fabs(dz) > 0.25 * boxSize_Z)
    {
      ray->current_cell = -3;			//deactivate ray
      ray->machnumber = 0;
      free_ray(ray);
      return;
    }
#endif /* NUMDIMS > 2 */
#endif

#ifdef SKIP_BORDER
  int reflective_boundary_reached = 0;
  if(boundary_reached)
    {
#ifdef REFLECTIVE_X
      if(boundary_reached & REFL_X_FLAGS)
        {
          reflective_boundary_reached = 1;
        }
#endif
#ifdef REFLECTIVE_Y
      if(boundary_reached & REFL_Y_FLAGS)
        {
          reflective_boundary_reached = 1;
        }
#endif
#ifdef REFLECTIVE_Z
      if(boundary_reached & REFL_Z_FLAGS)
        {
          reflective_boundary_reached = 1;
        }
#endif

      if(reflective_boundary_reached)
        {
          terminate("ray reached the border, increase DistToBorder in param.txt");
        }
    }
#endif

  ray->steps++;

  MOVEPRINTF("task: %d: moved ray to %f %f \tray now in cell: %f %f\n", ThisTask, ray->pos[0], ray->pos[1], P[ray->current_cell].Pos[0], P[ray->current_cell].Pos[1])
    //if the ray went maximum steps: collect quantities
    if(ray->steps == All.RayStepsMax)
    {
      if(!ray->preshock_dir)    //post-shock direction
        {
          set_post_shock_quant(ray, task, original_index);
        }
      else                      //pre-shock direction
        {
          calculate_machnumber(ray, this_task);
        }
    }

  if(this_task)
    {
      if((ray_wrong_direction(SDATA(ray->current_cell).ShockDir, ray, 1) || !SDATA(ray->current_cell).ShockZone) && !ray->preshock_dir) //ray reached the post shock region
        {
          //set the postshock quantity
          set_post_shock_quant(ray, task, original_index);

          //revert direction (send ray to preshock direction)

          ray->dir[0] *= -1;
          ray->dir[1] *= -1;
          ray->dir[2] *= -1;

          ray->preshock_dir = 1;
          ray->previous_cell = ray->current_cell;       //note: dir changed, therefore the current cell has to be ignored in find_next_voronoi_cell
          ray->previous_cell_task = ThisTask;
          ray->current_cell = old_cell;
          ray->steps--;
          ray->steps *= -1;
          return;
        }

      else if((ray_wrong_direction(SDATA(ray->current_cell).ShockDir, ray, 1) || !SDATA(ray->current_cell).ShockZone) && ray->preshock_dir)     //ray reached the pre shock region
        {
          calculate_machnumber(ray, 1);
          ray->current_cell = -1; //deactivate ray
          return;
        }
      else //ray is in the middle of the shock zone (with a direction in agreement with the shock direction and max steps)
        {
#ifdef USE_MAXIMUM_FOR_POST_SHOCK
          ray->rho_post_shock = fmax(ray->rho_post_shock, SphP[ray->current_cell].Density);
          ray->p_post_shock = fmax(ray->p_post_shock, SphP[ray->current_cell].Pressure);
#ifdef SHOCK_JUMP_T
          ray->post_shock_quant = fmax(ray->post_shock_quant, SDATA(ray->current_cell).Temperature);
#else
#error Use USE_MAXIMUM_FOR_POST_SHOCK in combination with SHOCK_JUMP_T.
#endif
#endif

#ifdef COSMIC_RAYS_SHOCK_ACCELERATION
          if(ray->steps <= 3 && !ray->preshock_dir)
            {
              assert(ray->post_shock_cells[ray->steps - 1] == -1);
              ray->post_shock_eth_tot_sum += injection_weight(ray->current_cell);
              ray->cell_count++;
              ray->post_shock_cells[ray->steps - 1] = ray->current_cell;
              ray->post_shock_cells_tasks[ray->steps - 1] = ThisTask;
            }
#endif

          myassert(SDATA(ray->current_cell).ShockZone);

          if(divergence(ray->current_cell) < ray->surface_value)
            {
              ray->current_cell = -3;   //deactivate ray
              ray->machnumber = 0;

              free_ray(ray);
              return;
            }

          ray->previous_cell = old_cell;
          ray->previous_cell_task = ThisTask;
        }
    }
  else                          //this_task == 0, ray crossed the local domain
    {
      MOVEPRINTF("\tRay crossed the domain!\n");

      if((ray_wrong_direction(SDATAEXCH(ray->current_cell).ShockDir, ray, 0) || !SDATAEXCH(ray->current_cell).ShockZone) && !ray->preshock_dir) //ray reached the post shock region
        {
          MOVEPRINTF("\tRay reached post shock region\n");

          //set the postshock quantity
          set_post_shock_quant(ray, task, original_index);

          //revert direction (send ray to preshock direction)

          ray->dir[0] *= -1;
          ray->dir[1] *= -1;
          ray->dir[2] *= -1;

          ray->preshock_dir = 1;
          ray->previous_cell = Mesh.DP[DC[edge].dp_index].originalindex;
          ray->previous_cell_task = DC[edge].task;
          ray->current_cell = old_cell;
          ray->steps--;
          ray->steps *= -1;
          return;

        }
      else if((ray_wrong_direction(SDATAEXCH(ray->current_cell).ShockDir, ray, 0) || !SDATAEXCH(ray->current_cell).ShockZone) && ray->preshock_dir)     //ray reached the pre shock region
        {
          MOVEPRINTF("\tRay reached pre shock region\n");

          calculate_machnumber(ray, 0);
          ray->current_cell = -1;       //deactivate ray
          return;
        }

      else                      //ray is in the middle of the shock zone
        {
#ifdef USE_MAXIMUM_FOR_POST_SHOCK
          ray->rho_post_shock = fmax(ray->rho_post_shock, PrimExch[ray->current_cell].Density);
          ray->p_post_shock = fmax(ray->p_post_shock, PrimExch[ray->current_cell].Pressure);
#ifdef SHOCK_JUMP_T
          ray->post_shock_quant = fmax(ray->post_shock_quant, SDATAEXCH(ray->current_cell).Temperature);
#else
#error Use USE_MAXIMUM_FOR_POST_SHOCK in combination with SHOCK_JUMP_T.
#endif
#endif

          MOVEPRINTF("\tRay is in the middle of the shock zone\n");

#ifdef COSMIC_RAYS_SHOCK_ACCELERATION
          if(ray->steps <= 3 && !ray->preshock_dir)
            {
              assert(ray->post_shock_cells[ray->steps - 1] == -1);
              ray->post_shock_eth_tot_sum += injection_weight_exch(ray->current_cell);
              ray->cell_count++;
              ray->post_shock_cells[ray->steps - 1] = original_index;
              ray->post_shock_cells_tasks[ray->steps - 1] = task;
            }
#endif

          myassert(SDATAEXCH(ray->current_cell).ShockZone);

          if(divergence_exchange(ray->current_cell) < ray->surface_value)
            {
              ray->current_cell = -3;   //deactivate ray
              ray->machnumber = 0;

              free_ray(ray);
              return;
            }

          ray->previous_cell = old_cell;
          ray->previous_cell_task = ThisTask;

          //export this ray. Note: ray gets deactivated (ray->current_cell = -1) in sf_exchange_rays
          ray->current_cell = Mesh.DP[DC[edge].dp_index].originalindex;
          rays_send_count[DC[edge].task]++;
          insert_export_ray(ray, DC[edge].task);
        }
    }
}

static void ray_action()
{
  //create rays
  int idx, k;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      k = TimeBinsHydro.ActiveParticleList[idx];
      if(k < 0)
        continue;

      if(SDATA(k).ShockZone)
        {
          //create ray and send it to postshock region
          insert_ray(SphP[k].Center[0], SphP[k].Center[1], SphP[k].Center[2], (-1 * SDATA(k).ShockDir[0]), (-1 * SDATA(k).ShockDir[1]), (-1 * SDATA(k).ShockDir[2]), divergence(k), k, 0);
        }
    }

#ifdef SHOCK_FINDER_POST_PROCESSING
  mpi_printf_each("\tTask %d: Sent %ld rays!\n", ThisTask, NewRay - Rays);
  unsigned int total_exported = 0;
#endif


  //move rays until they reach the preshockzone

  shock_ray *current_ray = Rays;
  unsigned int active_rays = 0;
  unsigned int total_active_rays = 0;
  int count = 0;


  //evolving the rays
  do
    {
      current_ray = Rays;
      active_rays = 0;

      for(k = 0; k < NTask; k++)
        {
          rays_send_count[k] = 0;
          rays_receive_count[k] = 0;
        }

      export_rays.counter = 0;

      while(current_ray != NewRay)
        {
          if(current_ray->current_cell >= 0)
            {
              active_rays++;
              move_ray(current_ray);
            }

          current_ray++;

        }

      //exchange the rays
      //

      unsigned int imported, exported;
      sf_exchange_rays(&imported, &exported);

      active_rays += imported;
      active_rays -= exported;

      MPI_Allreduce(&active_rays, &total_active_rays, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

#ifdef SHOCK_FINDER_POST_PROCESSING
      MPI_Reduce(&exported, &total_exported, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
      mpi_printf_each("\tTask %d: %d active rays\n", ThisTask, active_rays);
      mpi_printf("\tExported rays: %d\n", total_exported);
      mpi_printf("\t\tTotal active rays: %d\n", total_active_rays);
#endif

      count++;
      if(count >= 100)
        {
          terminate("Too much ray action!\n Hint: Use the flag SKIP_BORDER for reflective boundaries!\n");
        }

    }
  while(total_active_rays != 0);
}

static void set_shock_surface()
{
  shock_ray *current_ray = Rays;

  //communicate the rays to the tasks of the surface cells
  //

  unsigned int j;

  for(j = 0; j < NTask; j++)
    {
      rays_send_count[j] = 0;
      rays_receive_count[j] = 0;
    }

  export_rays.counter = 0;



  while(current_ray != NewRay)
    {
      assert(current_ray->current_cell <= -1);  //all rays should have been processed

      if(current_ray->current_cell == -1 && current_ray->original_cell_task != ThisTask)
        {
          rays_send_count[current_ray->original_cell_task]++;
          insert_export_ray(current_ray, current_ray->original_cell_task);
        }


      current_ray++;

    }

  unsigned int imported, exported;
  sf_exchange_rays(&imported, &exported);

  //mpi_printf_each("\tTask %d: imported: %d, exported: %d\n", ThisTask, imported, exported);


  //set shock surface
  //

  current_ray = Rays;

  while(current_ray != NewRay)
    {
      assert(current_ray->current_cell <= -1);  //all rays should have been processed

      if(current_ray->current_cell == -1 && current_ray->original_cell_task == ThisTask)
        {
          if(SDATA(current_ray->original_cell).ShockSurface != 1)
            {
              SDATA(current_ray->original_cell).ShockSurface = 1;
              particles_in_shock_surface++;

            }
        }

      current_ray++;
    }

#ifdef SHOCK_FINDER_POST_PROCESSING
  mpi_printf_each("\tTask %d: Particles in shock surface: %d\n", ThisTask, particles_in_shock_surface);
#endif
}

/*!
 * Reset wrong jumps and copy fields from the rays to SphP
 */
static void set_machnumber()
{
  shock_ray *current_ray = Rays;

  //set the mach numbers
  //

  current_ray = Rays;

  int wrong_jump = 0;

  while(current_ray != NewRay)
    {
      assert(current_ray->current_cell <= -1);  //all rays should have been processed

      if(current_ray->current_cell == -1 && SDATA(current_ray->original_cell).ShockSurface) //it's a potential shock
        {

          wrong_jump = 0;

          if(current_ray->machnumber == 0 || current_ray->machnumber == 1.0) //jump in wrong direction / no jump
            {
              wrong_jump = 1;
            }

#ifdef SHOCK_FINDER_FILTER_SFR

          if(current_ray == Rays)
            {
              mpi_printf("SHOCKS: Filtering effective equation of state cells (PhysDensThresh=%f) and inconsistent jumps\n", All.PhysDensThresh);
            }

            //for simulations with cooling and star formation: reset the shocks in the effective equation of state range,
            //and check consistency of jumps for high overdensities.
          if(!wrong_jump)
            {
             wrong_jump = in_effective_eos_range(SphP[current_ray->original_cell].Density, SDATA(current_ray->original_cell).Temperature);
            }


          if(!wrong_jump)
            {
              wrong_jump = in_effective_eos_range(current_ray->rho_pre_shock, current_ray->t_pre_shock);
            }

          if(!wrong_jump)
            {
#ifndef SHOCK_JUMP_T
#error current_ray->post_shock_quant != T_post_shock!
#endif
              //post_shock_quant = T_post_shock for the standard configuration
              wrong_jump = in_effective_eos_range(current_ray->rho_post_shock, current_ray->post_shock_quant);
            }
#endif


#ifdef SHOCK_FINDER_FILTER_INCONSISTENT
          if(!wrong_jump)
            {
              wrong_jump = inconsistent_jump(current_ray);
            }
#endif

          if(wrong_jump)
            {
              wrong_jump_direction++;

#ifdef RESET_WRONG_JUMPS
              SDATA(current_ray->original_cell).ShockSurface = 0;
              particles_in_shock_surface--;

              current_ray->machnumber = 0.0;
              current_ray->rho_pre_shock = 0.0;
              current_ray->rho_post_shock = 0.0;
              current_ray->p_post_shock = 0.0;
              current_ray->v_pre_shock[0] = 0.0;
              current_ray->v_pre_shock[1] = 0.0;
              current_ray->v_pre_shock[2] = 0.0;
              current_ray->v_post_shock[0] = 0.0;
              current_ray->v_post_shock[1] = 0.0;
              current_ray->v_post_shock[2] = 0.0;
              current_ray->dir[0] = 0.0;
              current_ray->dir[1] = 0.0;
              current_ray->dir[2] = 0.0;

#ifdef COSMIC_RAYS
              current_ray->gte = 0.0;
              current_ray->pth_post_shock = 0.0;
              current_ray->pcr_post_shock = 0.0;
              current_ray->pth_pre_shock = 0.0;
              current_ray->pcr_pre_shock = 0.0;
              current_ray->eth_post_shock = 0.0;
              current_ray->ecr_post_shock = 0.0;
#else
              current_ray->c_pre_shock = 0.0;
              current_ray->p_pre_shock = 0.0;
              current_ray->t_pre_shock = 0.0;
#endif
#endif
            }

          SphP[current_ray->original_cell].Machnumber = current_ray->machnumber;
          SDATA(current_ray->original_cell).ShockSurfaceArea = shock_surface_area(current_ray->original_cell);

#ifdef COSMIC_RAYS
          SphP[current_ray->original_cell].EnergyDissipation = current_ray->gte * shock_surface_area(current_ray->original_cell);
#else
          SDATA(current_ray->original_cell).CpreShock = current_ray->c_pre_shock;
          SDATA(current_ray->original_cell).PpreShock = current_ray->p_pre_shock;
          SDATA(current_ray->original_cell).PpostShock = current_ray->p_post_shock;
          SDATA(current_ray->original_cell).TpreShock = current_ray->t_pre_shock;

#endif
          SDATA(current_ray->original_cell).RhopreShock = current_ray->rho_pre_shock;
          SDATA(current_ray->original_cell).RhopostShock = current_ray->rho_post_shock;

          SDATA(current_ray->original_cell).VpreShock[0] = current_ray->v_pre_shock[0];
          SDATA(current_ray->original_cell).VpreShock[1] = current_ray->v_pre_shock[1];
          SDATA(current_ray->original_cell).VpreShock[2] = current_ray->v_pre_shock[2];

          SDATA(current_ray->original_cell).VpostShock[0] = current_ray->v_post_shock[0];
          SDATA(current_ray->original_cell).VpostShock[1] = current_ray->v_post_shock[1];
          SDATA(current_ray->original_cell).VpostShock[2] = current_ray->v_post_shock[2];
        }

      current_ray++;

    }
#ifdef COSMIC_RAYS_SHOCK_ACCELERATION
  accelerate_cosmic_rays_at_shock(Rays, NewRay);
#endif
}

static void sf_exchange_rays(unsigned int *imported, unsigned int *exported)
{
  unsigned int nimport = 0;
  unsigned int nexport = 0;
  unsigned int send_offset[NTask];
  unsigned int receive_offset[NTask];
  unsigned int j;

  //get rays_receive_count
  MPI_Alltoall(rays_send_count, 1, MPI_UNSIGNED, rays_receive_count, 1, MPI_UNSIGNED, MPI_COMM_WORLD);

  //calculate nimport/nexport, send_offset/receive_offset
  for(j = 0, send_offset[0] = 0, receive_offset[0] = 0; j < NTask; j++)
    {
      nimport += rays_receive_count[j];
      nexport += rays_send_count[j];

      if(j > 0)
        {
          send_offset[j] = send_offset[j - 1] + rays_send_count[j - 1];
          receive_offset[j] = receive_offset[j - 1] + rays_receive_count[j - 1];
        }
    }

  shock_ray *rays_in = (shock_ray *) mymalloc("rays_in", nimport * sizeof(shock_ray));
  shock_ray *rays_out = (shock_ray *) mymalloc("rays_out", nexport * sizeof(shock_ray));

  //prepare data for export
  for(j = 0; j < NTask; j++)
    {
      rays_send_count[j] = 0;
    }

  unsigned int off;

  for(j = 0; j < export_rays.counter; j++)
    {
      off = send_offset[export_rays.to_task[j]] + rays_send_count[export_rays.to_task[j]]++;
      copy_ray(export_rays.elements[j], &rays_out[off]);

      //deactivate exported ray
      export_rays.elements[j]->current_cell = -2;
      free_ray(export_rays.elements[j]);
    }

  //exchange data
  unsigned int ngrp, receive_task;
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      //send_task = ThisTask;
      receive_task = ThisTask ^ ngrp;

      if(receive_task < NTask)
        {
          if(rays_send_count[receive_task] > 0 || rays_receive_count[receive_task] > 0)
            {
              MPI_Sendrecv(&rays_out[send_offset[receive_task]], rays_send_count[receive_task] * sizeof(shock_ray), MPI_BYTE, receive_task, 0,
                           &rays_in[receive_offset[receive_task]], rays_receive_count[receive_task] * sizeof(shock_ray), MPI_BYTE, receive_task, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  //store the received rays
  for(j = 0; j < nimport; j++)
    {
      copy_ray(&rays_in[j], insert_ray_empty());
    }


  myfree(rays_out);
  myfree(rays_in);

  *imported = nimport;
  *exported = nexport;
}

static void copy_ray(shock_ray * origin, shock_ray * destination)
{
  destination->pos[0] = origin->pos[0];
  destination->pos[1] = origin->pos[1];
  destination->pos[2] = origin->pos[2];

  destination->dir[0] = origin->dir[0];
  destination->dir[1] = origin->dir[1];
  destination->dir[2] = origin->dir[2];


  destination->machnumber = origin->machnumber;
  destination->surface_value = origin->surface_value;

  destination->current_cell = origin->current_cell;
  destination->previous_cell = origin->previous_cell;
  destination->previous_cell_task = origin->previous_cell_task;
  destination->original_cell = origin->original_cell;
  destination->original_cell_task = origin->original_cell_task;
  destination->steps = origin->steps;
  destination->preshock_dir = origin->preshock_dir;
  destination->quantities_set = origin->quantities_set;
  destination->rho_pre_shock = origin->rho_pre_shock;
  destination->rho_post_shock = origin->rho_post_shock;
  destination->p_post_shock = origin->p_post_shock;
  destination->v_post_shock[0] = origin->v_post_shock[0];
  destination->v_post_shock[1] = origin->v_post_shock[1];
  destination->v_post_shock[2] = origin->v_post_shock[2];
  destination->v_pre_shock[0] = origin->v_pre_shock[0];
  destination->v_pre_shock[1] = origin->v_pre_shock[1];
  destination->v_pre_shock[2] = origin->v_pre_shock[2];

#ifdef COSMIC_RAYS
  destination->pth_post_shock = origin->pth_post_shock;
  destination->pcr_post_shock = origin->pcr_post_shock;
  destination->pth_pre_shock = origin->pth_pre_shock;
  destination->pcr_pre_shock = origin->pcr_pre_shock;
  destination->eth_post_shock = origin->eth_post_shock;
  destination->ecr_post_shock = origin->ecr_post_shock;
  destination->gte = origin->gte;
#ifdef COSMIC_RAYS_SHOCK_ACCELERATION
#ifdef COSMIC_RAYS_MAGNETIC_OBLIQUITY
  destination->magnetic_obliquity_pre_shock = origin->magnetic_obliquity_pre_shock;
#endif
  destination->eth_tot_pre_shock = origin->eth_tot_pre_shock;
  destination->post_shock_eth_tot_sum = origin->post_shock_eth_tot_sum;
  destination->cell_count = origin->cell_count;
  destination->post_shock_cells[0] = origin->post_shock_cells[0];
  destination->post_shock_cells[1] = origin->post_shock_cells[1];
  destination->post_shock_cells[2] = origin->post_shock_cells[2];
  destination->post_shock_cells_tasks[0] = origin->post_shock_cells_tasks[0];
  destination->post_shock_cells_tasks[1] = origin->post_shock_cells_tasks[1];
  destination->post_shock_cells_tasks[2] = origin->post_shock_cells_tasks[2];
#endif
#else
  destination->post_shock_quant = origin->post_shock_quant;
  destination->c_pre_shock = origin->c_pre_shock;
  destination->p_pre_shock = origin->p_pre_shock;
  destination->t_pre_shock = origin->t_pre_shock;
#endif
}

#ifdef SHOCK_FINDER_POST_PROCESSING
static void set_additional_info()
{
  long long int total_machnumber_calculated = 0;
  long long int total_particles_in_shock_zone = 0;
  long long int total_particles_in_shock_surface = 0;

  unsigned int total_average_steps_to_preshock = 0;
  double d_total_average_steps_to_preshock = 0;
  unsigned int total_wrong_jump_direction = 0;

  MPI_Allreduce(&particles_in_shock_surface, &total_particles_in_shock_surface, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
  mpi_printf("\tTotal number of particles in shock surface: %d\n", total_particles_in_shock_surface);

  MPI_Allreduce(&particles_in_shock_zone, &total_particles_in_shock_zone, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
  mpi_printf("\tTotal number of particles in shock zone: %d\n", total_particles_in_shock_zone);

  MPI_Reduce(&machnumber_calculated, &total_machnumber_calculated, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&average_steps_to_preshock, &total_average_steps_to_preshock, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  d_total_average_steps_to_preshock = (double) total_average_steps_to_preshock;
  d_total_average_steps_to_preshock /= total_machnumber_calculated;
  mpi_printf("\tAverage steps to pre-shock region: %f.\n", d_total_average_steps_to_preshock);

  MPI_Reduce(&wrong_jump_direction, &total_wrong_jump_direction, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  mpi_printf("\tCells with wrong jump direction: %d\n", total_wrong_jump_direction);

  SfVars.nof_shocksurface_part_total = total_particles_in_shock_surface;
  SfVars.nof_shockzone_part_total = total_particles_in_shock_zone;

  SfVars.nof_shocksurface_part_local = particles_in_shock_surface;
  SfVars.nof_shockzone_part_local = particles_in_shock_zone;
}
#endif

#endif
