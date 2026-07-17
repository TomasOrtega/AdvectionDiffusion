#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

typedef struct {
  double x;
  double y;
  double z;
} linear_advection_3d_coordinates;

typedef struct {
  double f;
} linear_advection_3d_state;

static inline bool
linear_advection_3d_state_valid(const linear_advection_3d_state *U)
{
  double f = U->f;
  
  return true;
}

typedef struct {
  double ax;
  double ay;
  double az;
} linear_advection_3d_parameters;

static inline bool
linear_advection_3d_parameters_valid(const linear_advection_3d_parameters *P)
{
  double ax = P->ax;
  double ay = P->ay;
  double az = P->az;

  return true;
}

typedef struct {
  double flux_f;
} linear_advection_3d_flux;

typedef struct {
  void (*x_flux)(const linear_advection_3d_coordinates *C, const linear_advection_3d_parameters *P,
    const linear_advection_3d_state *U, linear_advection_3d_flux *flux);
  void (*y_flux)(const linear_advection_3d_coordinates *C, const linear_advection_3d_parameters *P,
    const linear_advection_3d_state *U, linear_advection_3d_flux *flux);
  void (*z_flux)(const linear_advection_3d_coordinates *C, const linear_advection_3d_parameters *P,
    const linear_advection_3d_state *U, linear_advection_3d_flux *flux);
} linear_advection_3d_fluxes_3D;

static inline void
linear_advection_3d_x_flux(const linear_advection_3d_coordinates *C, const linear_advection_3d_parameters *P,
  const linear_advection_3d_state *U, linear_advection_3d_flux *flux)
{
  double x = C->x;
  double y = C->y;
  double z = C->z;
  double ax = P->ax;
  double ay = P->ay;
  double az = P->az;
  double f = U->f;

  flux->flux_f = (ax * f);
}

static inline void
linear_advection_3d_y_flux(const linear_advection_3d_coordinates *C, const linear_advection_3d_parameters *P,
  const linear_advection_3d_state *U, linear_advection_3d_flux *flux)
{
  double x = C->x;
  double y = C->y;
  double z = C->z;
  double ax = P->ax;
  double ay = P->ay;
  double az = P->az;
  double f = U->f;

  flux->flux_f = (ay * f);
}

static inline void
linear_advection_3d_z_flux(const linear_advection_3d_coordinates *C, const linear_advection_3d_parameters *P,
  const linear_advection_3d_state *U, linear_advection_3d_flux *flux)
{
  double x = C->x;
  double y = C->y;
  double z = C->z;
  double ax = P->ax;
  double ay = P->ay;
  double az = P->az;
  double f = U->f;

  flux->flux_f = (az * f);
}

static inline linear_advection_3d_fluxes_3D
create_fluxes()
{
  linear_advection_3d_fluxes_3D fluxes;

  fluxes.x_flux = linear_advection_3d_x_flux;
  fluxes.y_flux = linear_advection_3d_y_flux;
  fluxes.z_flux = linear_advection_3d_z_flux;
  
  return fluxes;
}

typedef struct {
  double mu1;
} linear_advection_3d_wavespeed;

typedef struct {
  void (*x_wavespeed)(const linear_advection_3d_coordinates *C, const linear_advection_3d_parameters *P,
    const linear_advection_3d_state *U, linear_advection_3d_wavespeed *wavespeed);
  void (*y_wavespeed)(const linear_advection_3d_coordinates *C, const linear_advection_3d_parameters *P,
    const linear_advection_3d_state *U, linear_advection_3d_wavespeed *wavespeed);
  void (*z_wavespeed)(const linear_advection_3d_coordinates *C, const linear_advection_3d_parameters *P,
    const linear_advection_3d_state *U, linear_advection_3d_wavespeed *wavespeed);
} linear_advection_3d_wavespeed_estimate_3D;

static inline void
linear_advection_3d_x_wavespeed(const linear_advection_3d_coordinates *C, const linear_advection_3d_parameters *P,
  const linear_advection_3d_state *U, linear_advection_3d_wavespeed *wavespeed)
{
  double x = C->x;
  double y = C->y;
  double z = C->z;
  double ax = P->ax;
  double ay = P->ay;
  double az = P->az;
  double f = U->f;

  wavespeed->mu1 = ax;
}

static inline void
linear_advection_3d_y_wavespeed(const linear_advection_3d_coordinates *C, const linear_advection_3d_parameters *P,
  const linear_advection_3d_state *U, linear_advection_3d_wavespeed *wavespeed)
{
  double x = C->x;
  double y = C->y;
  double z = C->z;
  double ax = P->ax;
  double ay = P->ay;
  double az = P->az;
  double f = U->f;

  wavespeed->mu1 = ay;
}

static inline void
linear_advection_3d_z_wavespeed(const linear_advection_3d_coordinates *C, const linear_advection_3d_parameters *P,
  const linear_advection_3d_state *U, linear_advection_3d_wavespeed *wavespeed)
{
  double x = C->x;
  double y = C->y;
  double z = C->z;
  double ax = P->ax;
  double ay = P->ay;
  double az = P->az;
  double f = U->f;

  wavespeed->mu1 = az;
}

static inline linear_advection_3d_wavespeed_estimate_3D
create_wavespeeds()
{
  linear_advection_3d_wavespeed_estimate_3D wavespeed_estimate;

  wavespeed_estimate.x_wavespeed = linear_advection_3d_x_wavespeed;
  wavespeed_estimate.y_wavespeed = linear_advection_3d_y_wavespeed;
  wavespeed_estimate.z_wavespeed = linear_advection_3d_z_wavespeed;
  
  return wavespeed_estimate;
}

typedef struct {
  double f_x;
  double f_y;
  double f_z;
} linear_advection_3d_gradient;

typedef struct {
  double diffusive_flux_f;
} linear_advection_3d_diffusive_flux;

typedef struct {
  void (*x_diffusive_flux)(const linear_advection_3d_coordinates *C, const linear_advection_3d_parameters *P,
    const linear_advection_3d_state *U, linear_advection_3d_gradient *DU, linear_advection_3d_diffusive_flux *flux);
  void (*y_diffusive_flux)(const linear_advection_3d_coordinates *C, const linear_advection_3d_parameters *P,
    const linear_advection_3d_state *U, linear_advection_3d_gradient *DU, linear_advection_3d_diffusive_flux *flux);
  void (*z_diffusive_flux)(const linear_advection_3d_coordinates *C, const linear_advection_3d_parameters *P,
    const linear_advection_3d_state *U, linear_advection_3d_gradient *DU, linear_advection_3d_diffusive_flux *flux);
} linear_advection_3d_diffusive_fluxes_3D;

static inline void
linear_advection_3d_x_diffusive_flux(const linear_advection_3d_coordinates *C, const linear_advection_3d_parameters *P,
  const linear_advection_3d_state *U, linear_advection_3d_gradient *DU, linear_advection_3d_diffusive_flux *diffusive_flux)
{
  double x = C->x;
  double y = C->y;
  double z = C->z;
  double ax = P->ax;
  double ay = P->ay;
  double az = P->az;
  double f = U->f;
  double f_x = DU->f_x;
  double f_y = DU->f_y;
  double f_z = DU->f_z;

  diffusive_flux->diffusive_flux_f = 0.0;
}

static inline void
linear_advection_3d_y_diffusive_flux(const linear_advection_3d_coordinates *C, const linear_advection_3d_parameters *P,
  const linear_advection_3d_state *U, linear_advection_3d_gradient *DU, linear_advection_3d_diffusive_flux *diffusive_flux)
{
  double x = C->x;
  double y = C->y;
  double z = C->z;
  double ax = P->ax;
  double ay = P->ay;
  double az = P->az;
  double f = U->f;
  double f_x = DU->f_x;
  double f_y = DU->f_y;
  double f_z = DU->f_z;

  diffusive_flux->diffusive_flux_f = 0.0;
}

static inline void
linear_advection_3d_z_diffusive_flux(const linear_advection_3d_coordinates *C, const linear_advection_3d_parameters *P,
  const linear_advection_3d_state *U, linear_advection_3d_gradient *DU, linear_advection_3d_diffusive_flux *diffusive_flux)
{
  double x = C->x;
  double y = C->y;
  double z = C->z;
  double ax = P->ax;
  double ay = P->ay;
  double az = P->az;
  double f = U->f;
  double f_x = DU->f_x;
  double f_y = DU->f_y;
  double f_z = DU->f_z;

  diffusive_flux->diffusive_flux_f = 0.0;
}

static inline linear_advection_3d_diffusive_fluxes_3D
create_diffusive_fluxes()
{
  linear_advection_3d_diffusive_fluxes_3D diffusive_fluxes;

  diffusive_fluxes.x_diffusive_flux = linear_advection_3d_x_diffusive_flux;
  diffusive_fluxes.y_diffusive_flux = linear_advection_3d_y_diffusive_flux;
  diffusive_fluxes.z_diffusive_flux = linear_advection_3d_z_diffusive_flux;
  
  return diffusive_fluxes;
}

typedef struct {
  double wave_f;
} linear_advection_3d_lax_friedrichs_wave;

typedef struct {
  linear_advection_3d_lax_friedrichs_wave wave1;
  linear_advection_3d_lax_friedrichs_wave wave2;
} linear_advection_3d_lax_friedrichs_wave_family;

typedef struct {
  void (*x_wave_family)(const linear_advection_3d_coordinates *C_L, const linear_advection_3d_coordinates *C_R,
    const linear_advection_3d_state *U_L, const linear_advection_3d_state *U_R,
    const linear_advection_3d_parameters *P, linear_advection_3d_lax_friedrichs_wave_family *wave_family);
  void (*y_wave_family)(const linear_advection_3d_coordinates *C_L, const linear_advection_3d_coordinates *C_R,
    const linear_advection_3d_state *U_L, const linear_advection_3d_state *U_R,
    const linear_advection_3d_parameters *P, linear_advection_3d_lax_friedrichs_wave_family *wave_family);
  void (*z_wave_family)(const linear_advection_3d_coordinates *C_L, const linear_advection_3d_coordinates *C_R,
    const linear_advection_3d_state *U_L, const linear_advection_3d_state *U_R,
    const linear_advection_3d_parameters *P, linear_advection_3d_lax_friedrichs_wave_family *wave_family);
} linear_advection_3d_lax_friedrichs_waves_3D;

static inline void
linear_advection_3d_lax_friedrichs_x_wave_family(const linear_advection_3d_coordinates *C_L, const linear_advection_3d_coordinates *C_R,
  const linear_advection_3d_state *U_L, const linear_advection_3d_state *U_R,
  const linear_advection_3d_parameters *P, linear_advection_3d_lax_friedrichs_wave_family *wave_family)
{
  double x_L = C_L->x;
  double x_R = C_R->x;
  double y_L = C_L->y;
  double y_R = C_R->y;
  double z_L = C_L->z;
  double z_R = C_R->z;
  double f_L = U_L->f;
  double f_R = U_R->f;
  double ax = P->ax;
  double ay = P->ay;
  double az = P->az;

  wave_family->wave1.wave_f = (0.5 * ((f_R - f_L) - (((ax * f_R) - (ax * f_L)) / fabs(ax))));

  wave_family->wave2.wave_f = (0.5 * ((f_R - f_L) + (((ax * f_R) - (ax * f_L)) / fabs(ax))));
}

static inline void
linear_advection_3d_lax_friedrichs_y_wave_family(const linear_advection_3d_coordinates *C_L, const linear_advection_3d_coordinates *C_R,
  const linear_advection_3d_state *U_L, const linear_advection_3d_state *U_R,
  const linear_advection_3d_parameters *P, linear_advection_3d_lax_friedrichs_wave_family *wave_family)
{
  double x_L = C_L->x;
  double x_R = C_R->x;
  double y_L = C_L->y;
  double y_R = C_R->y;
  double z_L = C_L->z;
  double z_R = C_R->z;
  double f_L = U_L->f;
  double f_R = U_R->f;
  double ax = P->ax;
  double ay = P->ay;
  double az = P->az;

  wave_family->wave1.wave_f = (0.5 * ((f_R - f_L) - (((ay * f_R) - (ay * f_L)) / fabs(ay))));

  wave_family->wave2.wave_f = (0.5 * ((f_R - f_L) + (((ay * f_R) - (ay * f_L)) / fabs(ay))));
}

static inline void
linear_advection_3d_lax_friedrichs_z_wave_family(const linear_advection_3d_coordinates *C_L, const linear_advection_3d_coordinates *C_R,
  const linear_advection_3d_state *U_L, const linear_advection_3d_state *U_R,
  const linear_advection_3d_parameters *P, linear_advection_3d_lax_friedrichs_wave_family *wave_family)
{
  double x_L = C_L->x;
  double x_R = C_R->x;
  double y_L = C_L->y;
  double y_R = C_R->y;
  double z_L = C_L->z;
  double z_R = C_R->z;
  double f_L = U_L->f;
  double f_R = U_R->f;
  double ax = P->ax;
  double ay = P->ay;
  double az = P->az;

  wave_family->wave1.wave_f = (0.5 * ((f_R - f_L) - (((az * f_R) - (az * f_L)) / fabs(az))));

  wave_family->wave2.wave_f = (0.5 * ((f_R - f_L) + (((az * f_R) - (az * f_L)) / fabs(az))));
}

static inline linear_advection_3d_lax_friedrichs_waves_3D
create_waves()
{
  linear_advection_3d_lax_friedrichs_waves_3D waves;

  waves.x_wave_family = linear_advection_3d_lax_friedrichs_x_wave_family;
  waves.y_wave_family = linear_advection_3d_lax_friedrichs_y_wave_family;
  waves.z_wave_family = linear_advection_3d_lax_friedrichs_z_wave_family;
  
  return waves;
}

static inline bool
linear_advection_3d_lax_friedrichs_x_waves_consistent(const linear_advection_3d_coordinates *C, const linear_advection_3d_state *U,
  const linear_advection_3d_parameters *P)
{
  linear_advection_3d_lax_friedrichs_wave_family wave_family;
  linear_advection_3d_lax_friedrichs_x_wave_family(C, C, U, U, P, &wave_family);

  double f_wave1 = wave_family.wave1.wave_f;
  double f_wave2 = wave_family.wave2.wave_f;

  return ((fabs(f_wave1) < 1.0e-8) && (fabs(f_wave2) < 1.0e-8));
}

static inline bool
linear_advection_3d_lax_friedrichs_y_waves_consistent(const linear_advection_3d_coordinates *C, const linear_advection_3d_state *U,
  const linear_advection_3d_parameters *P)
{
  linear_advection_3d_lax_friedrichs_wave_family wave_family;
  linear_advection_3d_lax_friedrichs_y_wave_family(C, C, U, U, P, &wave_family);

  double f_wave1 = wave_family.wave1.wave_f;
  double f_wave2 = wave_family.wave2.wave_f;

  return ((fabs(f_wave1) < 1.0e-8) && (fabs(f_wave2) < 1.0e-8));
}

static inline bool
linear_advection_3d_lax_friedrichs_z_waves_consistent(const linear_advection_3d_coordinates *C, const linear_advection_3d_state *U,
  const linear_advection_3d_parameters *P)
{
  linear_advection_3d_lax_friedrichs_wave_family wave_family;
  linear_advection_3d_lax_friedrichs_z_wave_family(C, C, U, U, P, &wave_family);

  double f_wave1 = wave_family.wave1.wave_f;
  double f_wave2 = wave_family.wave2.wave_f;

  return ((fabs(f_wave1) < 1.0e-8) && (fabs(f_wave2) < 1.0e-8));
}

static inline bool
linear_advection_3d_lax_friedrichs_x_waves_valid(const linear_advection_3d_coordinates *C_L, const linear_advection_3d_coordinates *C_R,
  const linear_advection_3d_state *U_L, const linear_advection_3d_state *U_R, const linear_advection_3d_parameters *P)
{
  double f_L = U_L->f;
  double f_R = U_R->f;

  double f_jump = f_R - f_L;

  linear_advection_3d_lax_friedrichs_wave_family wave_family;
  linear_advection_3d_lax_friedrichs_x_wave_family(C_L, C_R, U_L, U_R, P, &wave_family);

  double f_wave_sum = wave_family.wave1.wave_f + wave_family.wave2.wave_f;

  return (fabs(f_jump - f_wave_sum) < 1.0e-8);
}

static inline bool
linear_advection_3d_lax_friedrichs_y_waves_valid(const linear_advection_3d_coordinates *C_L, const linear_advection_3d_coordinates *C_R,
  const linear_advection_3d_state *U_L, const linear_advection_3d_state *U_R, const linear_advection_3d_parameters *P)
{
  double f_L = U_L->f;
  double f_R = U_R->f;

  double f_jump = f_R - f_L;

  linear_advection_3d_lax_friedrichs_wave_family wave_family;
  linear_advection_3d_lax_friedrichs_y_wave_family(C_L, C_R, U_L, U_R, P, &wave_family);

  double f_wave_sum = wave_family.wave1.wave_f + wave_family.wave2.wave_f;

  return (fabs(f_jump - f_wave_sum) < 1.0e-8);
}

static inline bool
linear_advection_3d_lax_friedrichs_z_waves_valid(const linear_advection_3d_coordinates *C_L, const linear_advection_3d_coordinates *C_R,
  const linear_advection_3d_state *U_L, const linear_advection_3d_state *U_R, const linear_advection_3d_parameters *P)
{
  double f_L = U_L->f;
  double f_R = U_R->f;

  double f_jump = f_R - f_L;

  linear_advection_3d_lax_friedrichs_wave_family wave_family;
  linear_advection_3d_lax_friedrichs_z_wave_family(C_L, C_R, U_L, U_R, P, &wave_family);

  double f_wave_sum = wave_family.wave1.wave_f + wave_family.wave2.wave_f;

  return (fabs(f_jump - f_wave_sum) < 1.0e-8);
}

typedef struct {
  double speed1;
  double speed2;
} linear_advection_3d_lax_friedrichs_speed_family;

typedef struct {
  void (*x_speed_family)(const linear_advection_3d_coordinates *C_L, const linear_advection_3d_coordinates *C_R,
    const linear_advection_3d_state *U_L, const linear_advection_3d_state *U_R,
    const linear_advection_3d_parameters *P, linear_advection_3d_lax_friedrichs_speed_family *speed_family);
  void (*y_speed_family)(const linear_advection_3d_coordinates *C_L, const linear_advection_3d_coordinates *C_R,
    const linear_advection_3d_state *U_L, const linear_advection_3d_state *U_R,
    const linear_advection_3d_parameters *P, linear_advection_3d_lax_friedrichs_speed_family *speed_family);
  void (*z_speed_family)(const linear_advection_3d_coordinates *C_L, const linear_advection_3d_coordinates *C_R,
    const linear_advection_3d_state *U_L, const linear_advection_3d_state *U_R,
    const linear_advection_3d_parameters *P, linear_advection_3d_lax_friedrichs_speed_family *speed_family);
} linear_advection_3d_lax_friedrichs_speeds_3D;

static inline void
linear_advection_3d_lax_friedrichs_x_speed_family(const linear_advection_3d_coordinates *C_L, const linear_advection_3d_coordinates *C_R,
  const linear_advection_3d_state *U_L, const linear_advection_3d_state *U_R,
  const linear_advection_3d_parameters *P, linear_advection_3d_lax_friedrichs_speed_family *speed_family)
{
  double x_L = C_L->x;
  double x_R = C_R->x;
  double y_L = C_L->y;
  double y_R = C_R->y;
  double z_L = C_L->z;
  double z_R = C_R->z;
  double f_L = U_L->f;
  double f_R = U_R->f;
  double ax = P->ax;
  double ay = P->ay;
  double az = P->az;

  speed_family->speed1 = -(fabs(ax));
  speed_family->speed2 = fabs(ax);
}

static inline void
linear_advection_3d_lax_friedrichs_y_speed_family(const linear_advection_3d_coordinates *C_L, const linear_advection_3d_coordinates *C_R,
  const linear_advection_3d_state *U_L, const linear_advection_3d_state *U_R,
  const linear_advection_3d_parameters *P, linear_advection_3d_lax_friedrichs_speed_family *speed_family)
{
  double x_L = C_L->x;
  double x_R = C_R->x;
  double y_L = C_L->y;
  double y_R = C_R->y;
  double z_L = C_L->z;
  double z_R = C_R->z;
  double f_L = U_L->f;
  double f_R = U_R->f;
  double ax = P->ax;
  double ay = P->ay;
  double az = P->az;

  speed_family->speed1 = -(fabs(ay));
  speed_family->speed2 = fabs(ay);
}

static inline void
linear_advection_3d_lax_friedrichs_z_speed_family(const linear_advection_3d_coordinates *C_L, const linear_advection_3d_coordinates *C_R,
  const linear_advection_3d_state *U_L, const linear_advection_3d_state *U_R,
  const linear_advection_3d_parameters *P, linear_advection_3d_lax_friedrichs_speed_family *speed_family)
{
  double x_L = C_L->x;
  double x_R = C_R->x;
  double y_L = C_L->y;
  double y_R = C_R->y;
  double z_L = C_L->z;
  double z_R = C_R->z;
  double f_L = U_L->f;
  double f_R = U_R->f;
  double ax = P->ax;
  double ay = P->ay;
  double az = P->az;

  speed_family->speed1 = -(fabs(az));
  speed_family->speed2 = fabs(az);
}

static inline linear_advection_3d_lax_friedrichs_speeds_3D
create_speeds()
{
  linear_advection_3d_lax_friedrichs_speeds_3D speeds;

  speeds.x_speed_family = linear_advection_3d_lax_friedrichs_x_speed_family;
  speeds.y_speed_family = linear_advection_3d_lax_friedrichs_y_speed_family;
  speeds.z_speed_family = linear_advection_3d_lax_friedrichs_z_speed_family;

  return speeds;
}

typedef struct {
  double fluctuation_f;
} linear_advection_3d_lax_friedrichs_fluctuation;

typedef struct {
  void (*x_fluctuation)(const linear_advection_3d_coordinates *C_L, const linear_advection_3d_coordinates *C_R,
    const linear_advection_3d_state *U_L, const linear_advection_3d_state *U_R,
    const linear_advection_3d_parameters *P, linear_advection_3d_lax_friedrichs_fluctuation *fluctuation);
  void (*y_fluctuation)(const linear_advection_3d_coordinates *C_L, const linear_advection_3d_coordinates *C_R,
    const linear_advection_3d_state *U_L, const linear_advection_3d_state *U_R,
    const linear_advection_3d_parameters *P, linear_advection_3d_lax_friedrichs_fluctuation *fluctuation);
  void (*z_fluctuation)(const linear_advection_3d_coordinates *C_L, const linear_advection_3d_coordinates *C_R,
    const linear_advection_3d_state *U_L, const linear_advection_3d_state *U_R,
    const linear_advection_3d_parameters *P, linear_advection_3d_lax_friedrichs_fluctuation *fluctuation);
} linear_advection_3d_lax_friedrichs_fluctuations_3D;

static inline void
linear_advection_3d_lax_friedrichs_x_left_fluctuation(const linear_advection_3d_coordinates *C_L, const linear_advection_3d_coordinates *C_R,
  const linear_advection_3d_state *U_L, const linear_advection_3d_state *U_R,
  const linear_advection_3d_parameters *P, linear_advection_3d_lax_friedrichs_fluctuation *left_fluctuation)
{
  linear_advection_3d_lax_friedrichs_wave_family wave_family;
  linear_advection_3d_lax_friedrichs_x_wave_family(C_L, C_R, U_L, U_R, P, &wave_family);

  linear_advection_3d_lax_friedrichs_speed_family speed_family;
  linear_advection_3d_lax_friedrichs_x_speed_family(C_L, C_R, U_L, U_R, P, &speed_family);

  left_fluctuation->fluctuation_f = 
    (wave_family.wave1.wave_f * fmin(speed_family.speed1, 0.0)) + 
    (wave_family.wave2.wave_f * fmin(speed_family.speed2, 0.0));
}

static inline void
linear_advection_3d_lax_friedrichs_y_left_fluctuation(const linear_advection_3d_coordinates *C_L, const linear_advection_3d_coordinates *C_R,
  const linear_advection_3d_state *U_L, const linear_advection_3d_state *U_R,
  const linear_advection_3d_parameters *P, linear_advection_3d_lax_friedrichs_fluctuation *left_fluctuation)
{
  linear_advection_3d_lax_friedrichs_wave_family wave_family;
  linear_advection_3d_lax_friedrichs_y_wave_family(C_L, C_R, U_L, U_R, P, &wave_family);

  linear_advection_3d_lax_friedrichs_speed_family speed_family;
  linear_advection_3d_lax_friedrichs_y_speed_family(C_L, C_R, U_L, U_R, P, &speed_family);

  left_fluctuation->fluctuation_f = 
    (wave_family.wave1.wave_f * fmin(speed_family.speed1, 0.0)) + 
    (wave_family.wave2.wave_f * fmin(speed_family.speed2, 0.0));
}

static inline void
linear_advection_3d_lax_friedrichs_z_left_fluctuation(const linear_advection_3d_coordinates *C_L, const linear_advection_3d_coordinates *C_R,
  const linear_advection_3d_state *U_L, const linear_advection_3d_state *U_R,
  const linear_advection_3d_parameters *P, linear_advection_3d_lax_friedrichs_fluctuation *left_fluctuation)
{
  linear_advection_3d_lax_friedrichs_wave_family wave_family;
  linear_advection_3d_lax_friedrichs_z_wave_family(C_L, C_R, U_L, U_R, P, &wave_family);

  linear_advection_3d_lax_friedrichs_speed_family speed_family;
  linear_advection_3d_lax_friedrichs_z_speed_family(C_L, C_R, U_L, U_R, P, &speed_family);

  left_fluctuation->fluctuation_f = 
    (wave_family.wave1.wave_f * fmin(speed_family.speed1, 0.0)) + 
    (wave_family.wave2.wave_f * fmin(speed_family.speed2, 0.0));
}

static inline void
linear_advection_3d_lax_friedrichs_x_right_fluctuation(const linear_advection_3d_coordinates *C_L, const linear_advection_3d_coordinates *C_R,
  const linear_advection_3d_state *U_L, const linear_advection_3d_state *U_R,
  const linear_advection_3d_parameters *P, linear_advection_3d_lax_friedrichs_fluctuation *right_fluctuation)
{
  linear_advection_3d_lax_friedrichs_wave_family wave_family;
  linear_advection_3d_lax_friedrichs_x_wave_family(C_L, C_R, U_L, U_R, P, &wave_family);

  linear_advection_3d_lax_friedrichs_speed_family speed_family;
  linear_advection_3d_lax_friedrichs_x_speed_family(C_L, C_R, U_L, U_R, P, &speed_family);

  right_fluctuation->fluctuation_f = 
    (wave_family.wave1.wave_f * fmax(speed_family.speed1, 0.0)) + 
    (wave_family.wave2.wave_f * fmax(speed_family.speed2, 0.0));
}

static inline void
linear_advection_3d_lax_friedrichs_y_right_fluctuation(const linear_advection_3d_coordinates *C_L, const linear_advection_3d_coordinates *C_R,
  const linear_advection_3d_state *U_L, const linear_advection_3d_state *U_R,
  const linear_advection_3d_parameters *P, linear_advection_3d_lax_friedrichs_fluctuation *right_fluctuation)
{
  linear_advection_3d_lax_friedrichs_wave_family wave_family;
  linear_advection_3d_lax_friedrichs_y_wave_family(C_L, C_R, U_L, U_R, P, &wave_family);

  linear_advection_3d_lax_friedrichs_speed_family speed_family;
  linear_advection_3d_lax_friedrichs_y_speed_family(C_L, C_R, U_L, U_R, P, &speed_family);

  right_fluctuation->fluctuation_f = 
    (wave_family.wave1.wave_f * fmax(speed_family.speed1, 0.0)) + 
    (wave_family.wave2.wave_f * fmax(speed_family.speed2, 0.0));
}

static inline void
linear_advection_3d_lax_friedrichs_z_right_fluctuation(const linear_advection_3d_coordinates *C_L, const linear_advection_3d_coordinates *C_R,
  const linear_advection_3d_state *U_L, const linear_advection_3d_state *U_R,
  const linear_advection_3d_parameters *P, linear_advection_3d_lax_friedrichs_fluctuation *right_fluctuation)
{
  linear_advection_3d_lax_friedrichs_wave_family wave_family;
  linear_advection_3d_lax_friedrichs_z_wave_family(C_L, C_R, U_L, U_R, P, &wave_family);

  linear_advection_3d_lax_friedrichs_speed_family speed_family;
  linear_advection_3d_lax_friedrichs_z_speed_family(C_L, C_R, U_L, U_R, P, &speed_family);

  right_fluctuation->fluctuation_f = 
    (wave_family.wave1.wave_f * fmax(speed_family.speed1, 0.0)) + 
    (wave_family.wave2.wave_f * fmax(speed_family.speed2, 0.0));
}

static inline linear_advection_3d_lax_friedrichs_fluctuations_3D
create_left_fluctuations()
{
  linear_advection_3d_lax_friedrichs_fluctuations_3D left_fluctuations;

  left_fluctuations.x_fluctuation = linear_advection_3d_lax_friedrichs_x_left_fluctuation;
  left_fluctuations.y_fluctuation = linear_advection_3d_lax_friedrichs_y_left_fluctuation;
  left_fluctuations.z_fluctuation = linear_advection_3d_lax_friedrichs_z_left_fluctuation;

  return left_fluctuations;
}

static inline linear_advection_3d_lax_friedrichs_fluctuations_3D
create_right_fluctuations()
{
  linear_advection_3d_lax_friedrichs_fluctuations_3D right_fluctuations;

  right_fluctuations.x_fluctuation = linear_advection_3d_lax_friedrichs_x_right_fluctuation;
  right_fluctuations.y_fluctuation = linear_advection_3d_lax_friedrichs_y_right_fluctuation;
  right_fluctuations.z_fluctuation = linear_advection_3d_lax_friedrichs_z_right_fluctuation;

  return right_fluctuations;
}

static inline bool
linear_advection_3d_lax_friedrichs_x_left_fluctuations_consistent(const linear_advection_3d_coordinates *C, const linear_advection_3d_state *U,
  const linear_advection_3d_parameters *P)
{
  linear_advection_3d_lax_friedrichs_fluctuation left_fluctuation;
  linear_advection_3d_lax_friedrichs_x_left_fluctuation(C, C, U, U, P, &left_fluctuation);

  double f_left_fluctuation = left_fluctuation.fluctuation_f;

  return ((fabs(f_left_fluctuation) < 1.0e-8));
}

static inline bool
linear_advection_3d_lax_friedrichs_y_left_fluctuations_consistent(const linear_advection_3d_coordinates *C, const linear_advection_3d_state *U,
  const linear_advection_3d_parameters *P)
{
  linear_advection_3d_lax_friedrichs_fluctuation left_fluctuation;
  linear_advection_3d_lax_friedrichs_y_left_fluctuation(C, C, U, U, P, &left_fluctuation);

  double f_left_fluctuation = left_fluctuation.fluctuation_f;

  return ((fabs(f_left_fluctuation) < 1.0e-8));
}

static inline bool
linear_advection_3d_lax_friedrichs_z_left_fluctuations_consistent(const linear_advection_3d_coordinates *C, const linear_advection_3d_state *U,
  const linear_advection_3d_parameters *P)
{
  linear_advection_3d_lax_friedrichs_fluctuation left_fluctuation;
  linear_advection_3d_lax_friedrichs_z_left_fluctuation(C, C, U, U, P, &left_fluctuation);

  double f_left_fluctuation = left_fluctuation.fluctuation_f;

  return ((fabs(f_left_fluctuation) < 1.0e-8));
}

static inline bool
linear_advection_3d_lax_friedrichs_x_right_fluctuations_consistent(const linear_advection_3d_coordinates *C, const linear_advection_3d_state *U,
  const linear_advection_3d_parameters *P)
{
  linear_advection_3d_lax_friedrichs_fluctuation right_fluctuation;
  linear_advection_3d_lax_friedrichs_x_right_fluctuation(C, C, U, U, P, &right_fluctuation);

  double f_right_fluctuation = right_fluctuation.fluctuation_f;

  return ((fabs(f_right_fluctuation) < 1.0e-8));
}

static inline bool
linear_advection_3d_lax_friedrichs_y_right_fluctuations_consistent(const linear_advection_3d_coordinates *C, const linear_advection_3d_state *U,
  const linear_advection_3d_parameters *P)
{
  linear_advection_3d_lax_friedrichs_fluctuation right_fluctuation;
  linear_advection_3d_lax_friedrichs_y_right_fluctuation(C, C, U, U, P, &right_fluctuation);

  double f_right_fluctuation = right_fluctuation.fluctuation_f;

  return ((fabs(f_right_fluctuation) < 1.0e-8));
}

static inline bool
linear_advection_3d_lax_friedrichs_z_right_fluctuations_consistent(const linear_advection_3d_coordinates *C, const linear_advection_3d_state *U,
  const linear_advection_3d_parameters *P)
{
  linear_advection_3d_lax_friedrichs_fluctuation right_fluctuation;
  linear_advection_3d_lax_friedrichs_z_right_fluctuation(C, C, U, U, P, &right_fluctuation);

  double f_right_fluctuation = right_fluctuation.fluctuation_f;

  return ((fabs(f_right_fluctuation) < 1.0e-8));
}

static inline bool
linear_advection_3d_lax_friedrichs_x_fluctuations_valid(const linear_advection_3d_coordinates *C_L, const linear_advection_3d_coordinates *C_R,
  const linear_advection_3d_state *U_L, const linear_advection_3d_state *U_R, const linear_advection_3d_parameters *P)
{
  linear_advection_3d_flux left_flux;
  linear_advection_3d_x_flux(C_L, P, U_L, &left_flux);

  linear_advection_3d_flux right_flux;
  linear_advection_3d_x_flux(C_R, P, U_R, &right_flux);

  double f_flux_jump = right_flux.flux_f - left_flux.flux_f;

  linear_advection_3d_lax_friedrichs_fluctuation left_fluctuation;
  linear_advection_3d_lax_friedrichs_x_left_fluctuation(C_L, C_R, U_L, U_R, P, &left_fluctuation);

  linear_advection_3d_lax_friedrichs_fluctuation right_fluctuation;
  linear_advection_3d_lax_friedrichs_x_right_fluctuation(C_L, C_R, U_L, U_R, P, &right_fluctuation);

  double f_fluctuation_sum = left_fluctuation.fluctuation_f + right_fluctuation.fluctuation_f;

  return (fabs(f_flux_jump - f_fluctuation_sum) < 1.0e-8);
}

static inline bool
linear_advection_3d_lax_friedrichs_y_fluctuations_valid(const linear_advection_3d_coordinates *C_L, const linear_advection_3d_coordinates *C_R,
  const linear_advection_3d_state *U_L, const linear_advection_3d_state *U_R, const linear_advection_3d_parameters *P)
{
  linear_advection_3d_flux left_flux;
  linear_advection_3d_y_flux(C_L, P, U_L, &left_flux);

  linear_advection_3d_flux right_flux;
  linear_advection_3d_y_flux(C_R, P, U_R, &right_flux);

  double f_flux_jump = right_flux.flux_f - left_flux.flux_f;

  linear_advection_3d_lax_friedrichs_fluctuation left_fluctuation;
  linear_advection_3d_lax_friedrichs_y_left_fluctuation(C_L, C_R, U_L, U_R, P, &left_fluctuation);

  linear_advection_3d_lax_friedrichs_fluctuation right_fluctuation;
  linear_advection_3d_lax_friedrichs_y_right_fluctuation(C_L, C_R, U_L, U_R, P, &right_fluctuation);

  double f_fluctuation_sum = left_fluctuation.fluctuation_f + right_fluctuation.fluctuation_f;

  return (fabs(f_flux_jump - f_fluctuation_sum) < 1.0e-8);
}

static inline bool
linear_advection_3d_lax_friedrichs_z_fluctuations_valid(const linear_advection_3d_coordinates *C_L, const linear_advection_3d_coordinates *C_R,
  const linear_advection_3d_state *U_L, const linear_advection_3d_state *U_R, const linear_advection_3d_parameters *P)
{
  linear_advection_3d_flux left_flux;
  linear_advection_3d_z_flux(C_L, P, U_L, &left_flux);

  linear_advection_3d_flux right_flux;
  linear_advection_3d_z_flux(C_R, P, U_R, &right_flux);

  double f_flux_jump = right_flux.flux_f - left_flux.flux_f;

  linear_advection_3d_lax_friedrichs_fluctuation left_fluctuation;
  linear_advection_3d_lax_friedrichs_z_left_fluctuation(C_L, C_R, U_L, U_R, P, &left_fluctuation);

  linear_advection_3d_lax_friedrichs_fluctuation right_fluctuation;
  linear_advection_3d_lax_friedrichs_z_right_fluctuation(C_L, C_R, U_L, U_R, P, &right_fluctuation);

  double f_fluctuation_sum = left_fluctuation.fluctuation_f + right_fluctuation.fluctuation_f;

  return (fabs(f_flux_jump - f_fluctuation_sum) < 1.0e-8);
}

typedef struct {
  double reconstruction_f;
} linear_advection_3d_minmod_reconstruction;

typedef struct {
  void (*x_reconstruction)(const linear_advection_3d_coordinates *C_L, const linear_advection_3d_coordinates *C, const linear_advection_3d_coordinates *C_R,
    const linear_advection_3d_state *U_L, const linear_advection_3d_state *U, const linear_advection_3d_state *U_R,
    const linear_advection_3d_parameters *P, linear_advection_3d_minmod_reconstruction *reconstruction);
  void (*y_reconstruction)(const linear_advection_3d_coordinates *C_L, const linear_advection_3d_coordinates *C, const linear_advection_3d_coordinates *C_R,
    const linear_advection_3d_state *U_L, const linear_advection_3d_state *U, const linear_advection_3d_state *U_R,
    const linear_advection_3d_parameters *P, linear_advection_3d_minmod_reconstruction *reconstruction);
  void (*z_reconstruction)(const linear_advection_3d_coordinates *C_L, const linear_advection_3d_coordinates *C, const linear_advection_3d_coordinates *C_R,
    const linear_advection_3d_state *U_L, const linear_advection_3d_state *U, const linear_advection_3d_state *U_R,
    const linear_advection_3d_parameters *P, linear_advection_3d_minmod_reconstruction *reconstruction);
} linear_advection_3d_minmod_reconstruction_3D;

static inline void
linear_advection_3d_minmod_x_left_reconstruction(const linear_advection_3d_coordinates *C_L, const linear_advection_3d_coordinates *C,
  const linear_advection_3d_coordinates *C_R, const linear_advection_3d_state *U_L, const linear_advection_3d_state *U, const linear_advection_3d_state *U_R,
  const linear_advection_3d_parameters *P, linear_advection_3d_minmod_reconstruction *left_reconstruction)
{
  double x_L = C_L->x;
  double x = C->x; 
  double x_R = C_R->x;
  double y_L = C_L->y;
  double y = C->y; 
  double y_R = C_R->y;
  double z_L = C_L->z;
  double z = C->z; 
  double z_R = C_R->z;
  double f_L = U_L->f;
  double f = U->f;
  double f_R = U_R->f;
  double ax = P->ax;
  double ay = P->ay;
  double az = P->az;

  left_reconstruction->reconstruction_f = (f - (0.5 * (fmax(0.0, fmin((f - f_L), (f_R - f))) + fmin(0.0, fmax((f - f_L),
    (f_R - f))))));
}

static inline void
linear_advection_3d_minmod_y_left_reconstruction(const linear_advection_3d_coordinates *C_L, const linear_advection_3d_coordinates *C,
  const linear_advection_3d_coordinates *C_R, const linear_advection_3d_state *U_L, const linear_advection_3d_state *U, const linear_advection_3d_state *U_R,
  const linear_advection_3d_parameters *P, linear_advection_3d_minmod_reconstruction *left_reconstruction)
{
  double x_L = C_L->x;
  double x = C->x; 
  double x_R = C_R->x;
  double y_L = C_L->y;
  double y = C->y; 
  double y_R = C_R->y;
  double z_L = C_L->z;
  double z = C->z; 
  double z_R = C_R->z;
  double f_L = U_L->f;
  double f = U->f;
  double f_R = U_R->f;
  double ax = P->ax;
  double ay = P->ay;
  double az = P->az;

  left_reconstruction->reconstruction_f = (f - (0.5 * (fmax(0.0, fmin((f - f_L), (f_R - f))) + fmin(0.0, fmax((f - f_L),
    (f_R - f))))));
}

static inline void
linear_advection_3d_minmod_z_left_reconstruction(const linear_advection_3d_coordinates *C_L, const linear_advection_3d_coordinates *C,
  const linear_advection_3d_coordinates *C_R, const linear_advection_3d_state *U_L, const linear_advection_3d_state *U, const linear_advection_3d_state *U_R,
  const linear_advection_3d_parameters *P, linear_advection_3d_minmod_reconstruction *left_reconstruction)
{
  double x_L = C_L->x;
  double x = C->x; 
  double x_R = C_R->x;
  double y_L = C_L->y;
  double y = C->y; 
  double y_R = C_R->y;
  double z_L = C_L->z;
  double z = C->z; 
  double z_R = C_R->z;
  double f_L = U_L->f;
  double f = U->f;
  double f_R = U_R->f;
  double ax = P->ax;
  double ay = P->ay;
  double az = P->az;

  left_reconstruction->reconstruction_f = (f - (0.5 * (fmax(0.0, fmin((f - f_L), (f_R - f))) + fmin(0.0, fmax((f - f_L),
    (f_R - f))))));
}

static inline void
linear_advection_3d_minmod_x_right_reconstruction(const linear_advection_3d_coordinates *C_L, const linear_advection_3d_coordinates *C,
  const linear_advection_3d_coordinates *C_R, const linear_advection_3d_state *U_L, const linear_advection_3d_state *U, const linear_advection_3d_state *U_R,
  const linear_advection_3d_parameters *P, linear_advection_3d_minmod_reconstruction *right_reconstruction)
{
  double x_L = C_L->x;
  double x = C->x; 
  double x_R = C_R->x;
  double y_L = C_L->y;
  double y = C->y; 
  double y_R = C_R->y;
  double z_L = C_L->z;
  double z = C->z; 
  double z_R = C_R->z;
  double f_L = U_L->f;
  double f = U->f;
  double f_R = U_R->f;
  double ax = P->ax;
  double ay = P->ay;
  double az = P->az;

  right_reconstruction->reconstruction_f = (f + (0.5 * (fmax(0.0, fmin((f - f_L), (f_R - f))) + fmin(0.0, fmax((f -
    f_L), (f_R - f))))));
}

static inline void
linear_advection_3d_minmod_y_right_reconstruction(const linear_advection_3d_coordinates *C_L, const linear_advection_3d_coordinates *C,
  const linear_advection_3d_coordinates *C_R, const linear_advection_3d_state *U_L, const linear_advection_3d_state *U, const linear_advection_3d_state *U_R,
  const linear_advection_3d_parameters *P, linear_advection_3d_minmod_reconstruction *right_reconstruction)
{
  double x_L = C_L->x;
  double x = C->x; 
  double x_R = C_R->x;
  double y_L = C_L->y;
  double y = C->y; 
  double y_R = C_R->y;
  double z_L = C_L->z;
  double z = C->z; 
  double z_R = C_R->z;
  double f_L = U_L->f;
  double f = U->f;
  double f_R = U_R->f;
  double ax = P->ax;
  double ay = P->ay;
  double az = P->az;

  right_reconstruction->reconstruction_f = (f + (0.5 * (fmax(0.0, fmin((f - f_L), (f_R - f))) + fmin(0.0, fmax((f -
    f_L), (f_R - f))))));
}

static inline void
linear_advection_3d_minmod_z_right_reconstruction(const linear_advection_3d_coordinates *C_L, const linear_advection_3d_coordinates *C,
  const linear_advection_3d_coordinates *C_R, const linear_advection_3d_state *U_L, const linear_advection_3d_state *U, const linear_advection_3d_state *U_R,
  const linear_advection_3d_parameters *P, linear_advection_3d_minmod_reconstruction *right_reconstruction)
{
  double x_L = C_L->x;
  double x = C->x; 
  double x_R = C_R->x;
  double y_L = C_L->y;
  double y = C->y; 
  double y_R = C_R->y;
  double z_L = C_L->z;
  double z = C->z; 
  double z_R = C_R->z;
  double f_L = U_L->f;
  double f = U->f;
  double f_R = U_R->f;
  double ax = P->ax;
  double ay = P->ay;
  double az = P->az;

  right_reconstruction->reconstruction_f = (f + (0.5 * (fmax(0.0, fmin((f - f_L), (f_R - f))) + fmin(0.0, fmax((f -
    f_L), (f_R - f))))));
}

static inline linear_advection_3d_minmod_reconstruction_3D
create_left_reconstruction()
{
  linear_advection_3d_minmod_reconstruction_3D left_reconstruction;

  left_reconstruction.x_reconstruction = linear_advection_3d_minmod_x_left_reconstruction;
  left_reconstruction.y_reconstruction = linear_advection_3d_minmod_y_left_reconstruction;
  left_reconstruction.z_reconstruction = linear_advection_3d_minmod_z_left_reconstruction;

  return left_reconstruction;
}

static inline linear_advection_3d_minmod_reconstruction_3D
create_right_reconstruction()
{
  linear_advection_3d_minmod_reconstruction_3D right_reconstruction;

  right_reconstruction.x_reconstruction = linear_advection_3d_minmod_x_right_reconstruction;
  right_reconstruction.y_reconstruction = linear_advection_3d_minmod_y_right_reconstruction;
  right_reconstruction.z_reconstruction = linear_advection_3d_minmod_z_right_reconstruction;

  return right_reconstruction;
}

static inline bool
linear_advection_3d_minmod_x_left_reconstruction_consistent(const linear_advection_3d_coordinates *C, const linear_advection_3d_state *U,
const linear_advection_3d_parameters *P)
{
  linear_advection_3d_minmod_reconstruction left_reconstruction;
  linear_advection_3d_minmod_x_left_reconstruction(C, C, C, U, U, U, P, &left_reconstruction);

  double f_left_reconstruction = left_reconstruction.reconstruction_f;
  double f = U->f;

  return ((fabs(f_left_reconstruction - f) < 1.0e-8));
}

static inline bool
linear_advection_3d_minmod_y_left_reconstruction_consistent(const linear_advection_3d_coordinates *C, const linear_advection_3d_state *U,
const linear_advection_3d_parameters *P)
{
  linear_advection_3d_minmod_reconstruction left_reconstruction;
  linear_advection_3d_minmod_y_left_reconstruction(C, C, C, U, U, U, P, &left_reconstruction);

  double f_left_reconstruction = left_reconstruction.reconstruction_f;
  double f = U->f;

  return ((fabs(f_left_reconstruction - f) < 1.0e-8));
}

static inline bool
linear_advection_3d_minmod_z_left_reconstruction_consistent(const linear_advection_3d_coordinates *C, const linear_advection_3d_state *U,
const linear_advection_3d_parameters *P)
{
  linear_advection_3d_minmod_reconstruction left_reconstruction;
  linear_advection_3d_minmod_z_left_reconstruction(C, C, C, U, U, U, P, &left_reconstruction);

  double f_left_reconstruction = left_reconstruction.reconstruction_f;
  double f = U->f;

  return ((fabs(f_left_reconstruction - f) < 1.0e-8));
}

static inline bool
linear_advection_3d_minmod_x_right_reconstruction_consistent(const linear_advection_3d_coordinates *C, const linear_advection_3d_state *U,
const linear_advection_3d_parameters *P)
{
  linear_advection_3d_minmod_reconstruction right_reconstruction;
  linear_advection_3d_minmod_x_right_reconstruction(C, C, C, U, U, U, P, &right_reconstruction);

  double f_right_reconstruction = right_reconstruction.reconstruction_f;
  double f = U->f;

  return ((fabs(f_right_reconstruction - f) < 1.0e-8));
}

static inline bool
linear_advection_3d_minmod_y_right_reconstruction_consistent(const linear_advection_3d_coordinates *C, const linear_advection_3d_state *U,
const linear_advection_3d_parameters *P)
{
  linear_advection_3d_minmod_reconstruction right_reconstruction;
  linear_advection_3d_minmod_y_right_reconstruction(C, C, C, U, U, U, P, &right_reconstruction);

  double f_right_reconstruction = right_reconstruction.reconstruction_f;
  double f = U->f;

  return ((fabs(f_right_reconstruction - f) < 1.0e-8));
}

static inline bool
linear_advection_3d_minmod_z_right_reconstruction_consistent(const linear_advection_3d_coordinates *C, const linear_advection_3d_state *U,
const linear_advection_3d_parameters *P)
{
  linear_advection_3d_minmod_reconstruction right_reconstruction;
  linear_advection_3d_minmod_z_right_reconstruction(C, C, C, U, U, U, P, &right_reconstruction);

  double f_right_reconstruction = right_reconstruction.reconstruction_f;
  double f = U->f;

  return ((fabs(f_right_reconstruction - f) < 1.0e-8));
}

int
main(int argc, char *argv[])
{
  // Insert simulation drivers here.
  
  return 0;
}
                