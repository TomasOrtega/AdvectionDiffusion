#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

typedef struct {
  double x;
  double y;
  double z;
} advection_diffusion_full_3d_coordinates;

typedef struct {
  double f;
} advection_diffusion_full_3d_state;

static inline bool
advection_diffusion_full_3d_state_valid(const advection_diffusion_full_3d_state *U)
{
  double f = U->f;
  
  return true;
}

typedef struct {
  double ax;
  double ay;
  double az;
  double Dxx;
  double Dxy;
  double Dxz;
  double Dyx;
  double Dyy;
  double Dyz;
  double Dzx;
  double Dzy;
  double Dzz;
} advection_diffusion_full_3d_parameters;

static inline bool
advection_diffusion_full_3d_parameters_valid(const advection_diffusion_full_3d_parameters *P)
{
  double ax = P->ax;
  double ay = P->ay;
  double az = P->az;
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dxz = P->Dxz;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;
  double Dyz = P->Dyz;
  double Dzx = P->Dzx;
  double Dzy = P->Dzy;
  double Dzz = P->Dzz;

  return (Dxx > 0.0 && Dyy > 0.0 && Dzz > 0.0);
}

typedef struct {
  double flux_f;
} advection_diffusion_full_3d_flux;

typedef struct {
  void (*x_flux)(const advection_diffusion_full_3d_coordinates *C, const advection_diffusion_full_3d_parameters *P,
    const advection_diffusion_full_3d_state *U, advection_diffusion_full_3d_flux *flux);
  void (*y_flux)(const advection_diffusion_full_3d_coordinates *C, const advection_diffusion_full_3d_parameters *P,
    const advection_diffusion_full_3d_state *U, advection_diffusion_full_3d_flux *flux);
  void (*z_flux)(const advection_diffusion_full_3d_coordinates *C, const advection_diffusion_full_3d_parameters *P,
    const advection_diffusion_full_3d_state *U, advection_diffusion_full_3d_flux *flux);
} advection_diffusion_full_3d_fluxes_3D;

static inline void
advection_diffusion_full_3d_x_flux(const advection_diffusion_full_3d_coordinates *C, const advection_diffusion_full_3d_parameters *P,
  const advection_diffusion_full_3d_state *U, advection_diffusion_full_3d_flux *flux)
{
  double x = C->x;
  double y = C->y;
  double z = C->z;
  double ax = P->ax;
  double ay = P->ay;
  double az = P->az;
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dxz = P->Dxz;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;
  double Dyz = P->Dyz;
  double Dzx = P->Dzx;
  double Dzy = P->Dzy;
  double Dzz = P->Dzz;
  double f = U->f;

  flux->flux_f = (ax * f);
}

static inline void
advection_diffusion_full_3d_y_flux(const advection_diffusion_full_3d_coordinates *C, const advection_diffusion_full_3d_parameters *P,
  const advection_diffusion_full_3d_state *U, advection_diffusion_full_3d_flux *flux)
{
  double x = C->x;
  double y = C->y;
  double z = C->z;
  double ax = P->ax;
  double ay = P->ay;
  double az = P->az;
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dxz = P->Dxz;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;
  double Dyz = P->Dyz;
  double Dzx = P->Dzx;
  double Dzy = P->Dzy;
  double Dzz = P->Dzz;
  double f = U->f;

  flux->flux_f = (ay * f);
}

static inline void
advection_diffusion_full_3d_z_flux(const advection_diffusion_full_3d_coordinates *C, const advection_diffusion_full_3d_parameters *P,
  const advection_diffusion_full_3d_state *U, advection_diffusion_full_3d_flux *flux)
{
  double x = C->x;
  double y = C->y;
  double z = C->z;
  double ax = P->ax;
  double ay = P->ay;
  double az = P->az;
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dxz = P->Dxz;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;
  double Dyz = P->Dyz;
  double Dzx = P->Dzx;
  double Dzy = P->Dzy;
  double Dzz = P->Dzz;
  double f = U->f;

  flux->flux_f = (az * f);
}

static inline advection_diffusion_full_3d_fluxes_3D
create_fluxes()
{
  advection_diffusion_full_3d_fluxes_3D fluxes;

  fluxes.x_flux = advection_diffusion_full_3d_x_flux;
  fluxes.y_flux = advection_diffusion_full_3d_y_flux;
  fluxes.z_flux = advection_diffusion_full_3d_z_flux;
  
  return fluxes;
}

typedef struct {
  double mu1;
} advection_diffusion_full_3d_wavespeed;

typedef struct {
  void (*x_wavespeed)(const advection_diffusion_full_3d_coordinates *C, const advection_diffusion_full_3d_parameters *P,
    const advection_diffusion_full_3d_state *U, advection_diffusion_full_3d_wavespeed *wavespeed);
  void (*y_wavespeed)(const advection_diffusion_full_3d_coordinates *C, const advection_diffusion_full_3d_parameters *P,
    const advection_diffusion_full_3d_state *U, advection_diffusion_full_3d_wavespeed *wavespeed);
  void (*z_wavespeed)(const advection_diffusion_full_3d_coordinates *C, const advection_diffusion_full_3d_parameters *P,
    const advection_diffusion_full_3d_state *U, advection_diffusion_full_3d_wavespeed *wavespeed);
} advection_diffusion_full_3d_wavespeed_estimate_3D;

static inline void
advection_diffusion_full_3d_x_wavespeed(const advection_diffusion_full_3d_coordinates *C, const advection_diffusion_full_3d_parameters *P,
  const advection_diffusion_full_3d_state *U, advection_diffusion_full_3d_wavespeed *wavespeed)
{
  double x = C->x;
  double y = C->y;
  double z = C->z;
  double ax = P->ax;
  double ay = P->ay;
  double az = P->az;
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dxz = P->Dxz;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;
  double Dyz = P->Dyz;
  double Dzx = P->Dzx;
  double Dzy = P->Dzy;
  double Dzz = P->Dzz;
  double f = U->f;

  wavespeed->mu1 = ax;
}

static inline void
advection_diffusion_full_3d_y_wavespeed(const advection_diffusion_full_3d_coordinates *C, const advection_diffusion_full_3d_parameters *P,
  const advection_diffusion_full_3d_state *U, advection_diffusion_full_3d_wavespeed *wavespeed)
{
  double x = C->x;
  double y = C->y;
  double z = C->z;
  double ax = P->ax;
  double ay = P->ay;
  double az = P->az;
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dxz = P->Dxz;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;
  double Dyz = P->Dyz;
  double Dzx = P->Dzx;
  double Dzy = P->Dzy;
  double Dzz = P->Dzz;
  double f = U->f;

  wavespeed->mu1 = ay;
}

static inline void
advection_diffusion_full_3d_z_wavespeed(const advection_diffusion_full_3d_coordinates *C, const advection_diffusion_full_3d_parameters *P,
  const advection_diffusion_full_3d_state *U, advection_diffusion_full_3d_wavespeed *wavespeed)
{
  double x = C->x;
  double y = C->y;
  double z = C->z;
  double ax = P->ax;
  double ay = P->ay;
  double az = P->az;
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dxz = P->Dxz;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;
  double Dyz = P->Dyz;
  double Dzx = P->Dzx;
  double Dzy = P->Dzy;
  double Dzz = P->Dzz;
  double f = U->f;

  wavespeed->mu1 = az;
}

static inline advection_diffusion_full_3d_wavespeed_estimate_3D
create_wavespeeds()
{
  advection_diffusion_full_3d_wavespeed_estimate_3D wavespeed_estimate;

  wavespeed_estimate.x_wavespeed = advection_diffusion_full_3d_x_wavespeed;
  wavespeed_estimate.y_wavespeed = advection_diffusion_full_3d_y_wavespeed;
  wavespeed_estimate.z_wavespeed = advection_diffusion_full_3d_z_wavespeed;
  
  return wavespeed_estimate;
}

typedef struct {
  double f_x;
  double f_y;
  double f_z;
} advection_diffusion_full_3d_gradient;

typedef struct {
  double diffusive_flux_f;
} advection_diffusion_full_3d_diffusive_flux;

typedef struct {
  void (*x_diffusive_flux)(const advection_diffusion_full_3d_coordinates *C, const advection_diffusion_full_3d_parameters *P,
    const advection_diffusion_full_3d_state *U, advection_diffusion_full_3d_gradient *DU, advection_diffusion_full_3d_diffusive_flux *flux);
  void (*y_diffusive_flux)(const advection_diffusion_full_3d_coordinates *C, const advection_diffusion_full_3d_parameters *P,
    const advection_diffusion_full_3d_state *U, advection_diffusion_full_3d_gradient *DU, advection_diffusion_full_3d_diffusive_flux *flux);
  void (*z_diffusive_flux)(const advection_diffusion_full_3d_coordinates *C, const advection_diffusion_full_3d_parameters *P,
    const advection_diffusion_full_3d_state *U, advection_diffusion_full_3d_gradient *DU, advection_diffusion_full_3d_diffusive_flux *flux);
} advection_diffusion_full_3d_diffusive_fluxes_3D;

static inline void
advection_diffusion_full_3d_x_diffusive_flux(const advection_diffusion_full_3d_coordinates *C, const advection_diffusion_full_3d_parameters *P,
  const advection_diffusion_full_3d_state *U, advection_diffusion_full_3d_gradient *DU, advection_diffusion_full_3d_diffusive_flux *diffusive_flux)
{
  double x = C->x;
  double y = C->y;
  double z = C->z;
  double ax = P->ax;
  double ay = P->ay;
  double az = P->az;
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dxz = P->Dxz;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;
  double Dyz = P->Dyz;
  double Dzx = P->Dzx;
  double Dzy = P->Dzy;
  double Dzz = P->Dzz;
  double f = U->f;
  double f_x = DU->f_x;
  double f_y = DU->f_y;
  double f_z = DU->f_z;

  diffusive_flux->diffusive_flux_f = ((Dxx * f_x) + (Dxy * f_y) + (Dxz * f_z));
}

static inline void
advection_diffusion_full_3d_y_diffusive_flux(const advection_diffusion_full_3d_coordinates *C, const advection_diffusion_full_3d_parameters *P,
  const advection_diffusion_full_3d_state *U, advection_diffusion_full_3d_gradient *DU, advection_diffusion_full_3d_diffusive_flux *diffusive_flux)
{
  double x = C->x;
  double y = C->y;
  double z = C->z;
  double ax = P->ax;
  double ay = P->ay;
  double az = P->az;
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dxz = P->Dxz;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;
  double Dyz = P->Dyz;
  double Dzx = P->Dzx;
  double Dzy = P->Dzy;
  double Dzz = P->Dzz;
  double f = U->f;
  double f_x = DU->f_x;
  double f_y = DU->f_y;
  double f_z = DU->f_z;

  diffusive_flux->diffusive_flux_f = ((Dyx * f_x) + (Dyy * f_y) + (Dyz * f_z));
}

static inline void
advection_diffusion_full_3d_z_diffusive_flux(const advection_diffusion_full_3d_coordinates *C, const advection_diffusion_full_3d_parameters *P,
  const advection_diffusion_full_3d_state *U, advection_diffusion_full_3d_gradient *DU, advection_diffusion_full_3d_diffusive_flux *diffusive_flux)
{
  double x = C->x;
  double y = C->y;
  double z = C->z;
  double ax = P->ax;
  double ay = P->ay;
  double az = P->az;
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dxz = P->Dxz;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;
  double Dyz = P->Dyz;
  double Dzx = P->Dzx;
  double Dzy = P->Dzy;
  double Dzz = P->Dzz;
  double f = U->f;
  double f_x = DU->f_x;
  double f_y = DU->f_y;
  double f_z = DU->f_z;

  diffusive_flux->diffusive_flux_f = ((Dzx * f_x) + (Dzy * f_y) + (Dzz * f_z));
}

static inline advection_diffusion_full_3d_diffusive_fluxes_3D
create_diffusive_fluxes()
{
  advection_diffusion_full_3d_diffusive_fluxes_3D diffusive_fluxes;

  diffusive_fluxes.x_diffusive_flux = advection_diffusion_full_3d_x_diffusive_flux;
  diffusive_fluxes.y_diffusive_flux = advection_diffusion_full_3d_y_diffusive_flux;
  diffusive_fluxes.z_diffusive_flux = advection_diffusion_full_3d_z_diffusive_flux;
  
  return diffusive_fluxes;
}

typedef struct {
  double wave_f;
} advection_diffusion_full_3d_lax_friedrichs_wave;

typedef struct {
  advection_diffusion_full_3d_lax_friedrichs_wave wave1;
  advection_diffusion_full_3d_lax_friedrichs_wave wave2;
} advection_diffusion_full_3d_lax_friedrichs_wave_family;

typedef struct {
  void (*x_wave_family)(const advection_diffusion_full_3d_coordinates *C_L, const advection_diffusion_full_3d_coordinates *C_R,
    const advection_diffusion_full_3d_state *U_L, const advection_diffusion_full_3d_state *U_R,
    const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_lax_friedrichs_wave_family *wave_family);
  void (*y_wave_family)(const advection_diffusion_full_3d_coordinates *C_L, const advection_diffusion_full_3d_coordinates *C_R,
    const advection_diffusion_full_3d_state *U_L, const advection_diffusion_full_3d_state *U_R,
    const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_lax_friedrichs_wave_family *wave_family);
  void (*z_wave_family)(const advection_diffusion_full_3d_coordinates *C_L, const advection_diffusion_full_3d_coordinates *C_R,
    const advection_diffusion_full_3d_state *U_L, const advection_diffusion_full_3d_state *U_R,
    const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_lax_friedrichs_wave_family *wave_family);
} advection_diffusion_full_3d_lax_friedrichs_waves_3D;

static inline void
advection_diffusion_full_3d_lax_friedrichs_x_wave_family(const advection_diffusion_full_3d_coordinates *C_L, const advection_diffusion_full_3d_coordinates *C_R,
  const advection_diffusion_full_3d_state *U_L, const advection_diffusion_full_3d_state *U_R,
  const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_lax_friedrichs_wave_family *wave_family)
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
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dxz = P->Dxz;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;
  double Dyz = P->Dyz;
  double Dzx = P->Dzx;
  double Dzy = P->Dzy;
  double Dzz = P->Dzz;

  wave_family->wave1.wave_f = (0.5 * ((f_R - f_L) - (((ax * f_R) - (ax * f_L)) / fabs(ax))));

  wave_family->wave2.wave_f = (0.5 * ((f_R - f_L) + (((ax * f_R) - (ax * f_L)) / fabs(ax))));
}

static inline void
advection_diffusion_full_3d_lax_friedrichs_y_wave_family(const advection_diffusion_full_3d_coordinates *C_L, const advection_diffusion_full_3d_coordinates *C_R,
  const advection_diffusion_full_3d_state *U_L, const advection_diffusion_full_3d_state *U_R,
  const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_lax_friedrichs_wave_family *wave_family)
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
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dxz = P->Dxz;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;
  double Dyz = P->Dyz;
  double Dzx = P->Dzx;
  double Dzy = P->Dzy;
  double Dzz = P->Dzz;

  wave_family->wave1.wave_f = (0.5 * ((f_R - f_L) - (((ay * f_R) - (ay * f_L)) / fabs(ay))));

  wave_family->wave2.wave_f = (0.5 * ((f_R - f_L) + (((ay * f_R) - (ay * f_L)) / fabs(ay))));
}

static inline void
advection_diffusion_full_3d_lax_friedrichs_z_wave_family(const advection_diffusion_full_3d_coordinates *C_L, const advection_diffusion_full_3d_coordinates *C_R,
  const advection_diffusion_full_3d_state *U_L, const advection_diffusion_full_3d_state *U_R,
  const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_lax_friedrichs_wave_family *wave_family)
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
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dxz = P->Dxz;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;
  double Dyz = P->Dyz;
  double Dzx = P->Dzx;
  double Dzy = P->Dzy;
  double Dzz = P->Dzz;

  wave_family->wave1.wave_f = (0.5 * ((f_R - f_L) - (((az * f_R) - (az * f_L)) / fabs(az))));

  wave_family->wave2.wave_f = (0.5 * ((f_R - f_L) + (((az * f_R) - (az * f_L)) / fabs(az))));
}

static inline advection_diffusion_full_3d_lax_friedrichs_waves_3D
create_waves()
{
  advection_diffusion_full_3d_lax_friedrichs_waves_3D waves;

  waves.x_wave_family = advection_diffusion_full_3d_lax_friedrichs_x_wave_family;
  waves.y_wave_family = advection_diffusion_full_3d_lax_friedrichs_y_wave_family;
  waves.z_wave_family = advection_diffusion_full_3d_lax_friedrichs_z_wave_family;
  
  return waves;
}

static inline bool
advection_diffusion_full_3d_lax_friedrichs_x_waves_consistent(const advection_diffusion_full_3d_coordinates *C, const advection_diffusion_full_3d_state *U,
  const advection_diffusion_full_3d_parameters *P)
{
  advection_diffusion_full_3d_lax_friedrichs_wave_family wave_family;
  advection_diffusion_full_3d_lax_friedrichs_x_wave_family(C, C, U, U, P, &wave_family);

  double f_wave1 = wave_family.wave1.wave_f;
  double f_wave2 = wave_family.wave2.wave_f;

  return ((fabs(f_wave1) < 1.0e-8) && (fabs(f_wave2) < 1.0e-8));
}

static inline bool
advection_diffusion_full_3d_lax_friedrichs_y_waves_consistent(const advection_diffusion_full_3d_coordinates *C, const advection_diffusion_full_3d_state *U,
  const advection_diffusion_full_3d_parameters *P)
{
  advection_diffusion_full_3d_lax_friedrichs_wave_family wave_family;
  advection_diffusion_full_3d_lax_friedrichs_y_wave_family(C, C, U, U, P, &wave_family);

  double f_wave1 = wave_family.wave1.wave_f;
  double f_wave2 = wave_family.wave2.wave_f;

  return ((fabs(f_wave1) < 1.0e-8) && (fabs(f_wave2) < 1.0e-8));
}

static inline bool
advection_diffusion_full_3d_lax_friedrichs_z_waves_consistent(const advection_diffusion_full_3d_coordinates *C, const advection_diffusion_full_3d_state *U,
  const advection_diffusion_full_3d_parameters *P)
{
  advection_diffusion_full_3d_lax_friedrichs_wave_family wave_family;
  advection_diffusion_full_3d_lax_friedrichs_z_wave_family(C, C, U, U, P, &wave_family);

  double f_wave1 = wave_family.wave1.wave_f;
  double f_wave2 = wave_family.wave2.wave_f;

  return ((fabs(f_wave1) < 1.0e-8) && (fabs(f_wave2) < 1.0e-8));
}

static inline bool
advection_diffusion_full_3d_lax_friedrichs_x_waves_valid(const advection_diffusion_full_3d_coordinates *C_L, const advection_diffusion_full_3d_coordinates *C_R,
  const advection_diffusion_full_3d_state *U_L, const advection_diffusion_full_3d_state *U_R, const advection_diffusion_full_3d_parameters *P)
{
  double f_L = U_L->f;
  double f_R = U_R->f;

  double f_jump = f_R - f_L;

  advection_diffusion_full_3d_lax_friedrichs_wave_family wave_family;
  advection_diffusion_full_3d_lax_friedrichs_x_wave_family(C_L, C_R, U_L, U_R, P, &wave_family);

  double f_wave_sum = wave_family.wave1.wave_f + wave_family.wave2.wave_f;

  return (fabs(f_jump - f_wave_sum) < 1.0e-8);
}

static inline bool
advection_diffusion_full_3d_lax_friedrichs_y_waves_valid(const advection_diffusion_full_3d_coordinates *C_L, const advection_diffusion_full_3d_coordinates *C_R,
  const advection_diffusion_full_3d_state *U_L, const advection_diffusion_full_3d_state *U_R, const advection_diffusion_full_3d_parameters *P)
{
  double f_L = U_L->f;
  double f_R = U_R->f;

  double f_jump = f_R - f_L;

  advection_diffusion_full_3d_lax_friedrichs_wave_family wave_family;
  advection_diffusion_full_3d_lax_friedrichs_y_wave_family(C_L, C_R, U_L, U_R, P, &wave_family);

  double f_wave_sum = wave_family.wave1.wave_f + wave_family.wave2.wave_f;

  return (fabs(f_jump - f_wave_sum) < 1.0e-8);
}

static inline bool
advection_diffusion_full_3d_lax_friedrichs_z_waves_valid(const advection_diffusion_full_3d_coordinates *C_L, const advection_diffusion_full_3d_coordinates *C_R,
  const advection_diffusion_full_3d_state *U_L, const advection_diffusion_full_3d_state *U_R, const advection_diffusion_full_3d_parameters *P)
{
  double f_L = U_L->f;
  double f_R = U_R->f;

  double f_jump = f_R - f_L;

  advection_diffusion_full_3d_lax_friedrichs_wave_family wave_family;
  advection_diffusion_full_3d_lax_friedrichs_z_wave_family(C_L, C_R, U_L, U_R, P, &wave_family);

  double f_wave_sum = wave_family.wave1.wave_f + wave_family.wave2.wave_f;

  return (fabs(f_jump - f_wave_sum) < 1.0e-8);
}

typedef struct {
  double speed1;
  double speed2;
} advection_diffusion_full_3d_lax_friedrichs_speed_family;

typedef struct {
  void (*x_speed_family)(const advection_diffusion_full_3d_coordinates *C_L, const advection_diffusion_full_3d_coordinates *C_R,
    const advection_diffusion_full_3d_state *U_L, const advection_diffusion_full_3d_state *U_R,
    const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_lax_friedrichs_speed_family *speed_family);
  void (*y_speed_family)(const advection_diffusion_full_3d_coordinates *C_L, const advection_diffusion_full_3d_coordinates *C_R,
    const advection_diffusion_full_3d_state *U_L, const advection_diffusion_full_3d_state *U_R,
    const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_lax_friedrichs_speed_family *speed_family);
  void (*z_speed_family)(const advection_diffusion_full_3d_coordinates *C_L, const advection_diffusion_full_3d_coordinates *C_R,
    const advection_diffusion_full_3d_state *U_L, const advection_diffusion_full_3d_state *U_R,
    const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_lax_friedrichs_speed_family *speed_family);
} advection_diffusion_full_3d_lax_friedrichs_speeds_3D;

static inline void
advection_diffusion_full_3d_lax_friedrichs_x_speed_family(const advection_diffusion_full_3d_coordinates *C_L, const advection_diffusion_full_3d_coordinates *C_R,
  const advection_diffusion_full_3d_state *U_L, const advection_diffusion_full_3d_state *U_R,
  const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_lax_friedrichs_speed_family *speed_family)
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
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dxz = P->Dxz;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;
  double Dyz = P->Dyz;
  double Dzx = P->Dzx;
  double Dzy = P->Dzy;
  double Dzz = P->Dzz;

  speed_family->speed1 = -(fabs(ax));
  speed_family->speed2 = fabs(ax);
}

static inline void
advection_diffusion_full_3d_lax_friedrichs_y_speed_family(const advection_diffusion_full_3d_coordinates *C_L, const advection_diffusion_full_3d_coordinates *C_R,
  const advection_diffusion_full_3d_state *U_L, const advection_diffusion_full_3d_state *U_R,
  const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_lax_friedrichs_speed_family *speed_family)
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
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dxz = P->Dxz;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;
  double Dyz = P->Dyz;
  double Dzx = P->Dzx;
  double Dzy = P->Dzy;
  double Dzz = P->Dzz;

  speed_family->speed1 = -(fabs(ay));
  speed_family->speed2 = fabs(ay);
}

static inline void
advection_diffusion_full_3d_lax_friedrichs_z_speed_family(const advection_diffusion_full_3d_coordinates *C_L, const advection_diffusion_full_3d_coordinates *C_R,
  const advection_diffusion_full_3d_state *U_L, const advection_diffusion_full_3d_state *U_R,
  const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_lax_friedrichs_speed_family *speed_family)
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
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dxz = P->Dxz;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;
  double Dyz = P->Dyz;
  double Dzx = P->Dzx;
  double Dzy = P->Dzy;
  double Dzz = P->Dzz;

  speed_family->speed1 = -(fabs(az));
  speed_family->speed2 = fabs(az);
}

static inline advection_diffusion_full_3d_lax_friedrichs_speeds_3D
create_speeds()
{
  advection_diffusion_full_3d_lax_friedrichs_speeds_3D speeds;

  speeds.x_speed_family = advection_diffusion_full_3d_lax_friedrichs_x_speed_family;
  speeds.y_speed_family = advection_diffusion_full_3d_lax_friedrichs_y_speed_family;
  speeds.z_speed_family = advection_diffusion_full_3d_lax_friedrichs_z_speed_family;

  return speeds;
}

typedef struct {
  double fluctuation_f;
} advection_diffusion_full_3d_lax_friedrichs_fluctuation;

typedef struct {
  void (*x_fluctuation)(const advection_diffusion_full_3d_coordinates *C_L, const advection_diffusion_full_3d_coordinates *C_R,
    const advection_diffusion_full_3d_state *U_L, const advection_diffusion_full_3d_state *U_R,
    const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_lax_friedrichs_fluctuation *fluctuation);
  void (*y_fluctuation)(const advection_diffusion_full_3d_coordinates *C_L, const advection_diffusion_full_3d_coordinates *C_R,
    const advection_diffusion_full_3d_state *U_L, const advection_diffusion_full_3d_state *U_R,
    const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_lax_friedrichs_fluctuation *fluctuation);
  void (*z_fluctuation)(const advection_diffusion_full_3d_coordinates *C_L, const advection_diffusion_full_3d_coordinates *C_R,
    const advection_diffusion_full_3d_state *U_L, const advection_diffusion_full_3d_state *U_R,
    const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_lax_friedrichs_fluctuation *fluctuation);
} advection_diffusion_full_3d_lax_friedrichs_fluctuations_3D;

static inline void
advection_diffusion_full_3d_lax_friedrichs_x_left_fluctuation(const advection_diffusion_full_3d_coordinates *C_L, const advection_diffusion_full_3d_coordinates *C_R,
  const advection_diffusion_full_3d_state *U_L, const advection_diffusion_full_3d_state *U_R,
  const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_lax_friedrichs_fluctuation *left_fluctuation)
{
  advection_diffusion_full_3d_lax_friedrichs_wave_family wave_family;
  advection_diffusion_full_3d_lax_friedrichs_x_wave_family(C_L, C_R, U_L, U_R, P, &wave_family);

  advection_diffusion_full_3d_lax_friedrichs_speed_family speed_family;
  advection_diffusion_full_3d_lax_friedrichs_x_speed_family(C_L, C_R, U_L, U_R, P, &speed_family);

  left_fluctuation->fluctuation_f = 
    (wave_family.wave1.wave_f * fmin(speed_family.speed1, 0.0)) + 
    (wave_family.wave2.wave_f * fmin(speed_family.speed2, 0.0));
}

static inline void
advection_diffusion_full_3d_lax_friedrichs_y_left_fluctuation(const advection_diffusion_full_3d_coordinates *C_L, const advection_diffusion_full_3d_coordinates *C_R,
  const advection_diffusion_full_3d_state *U_L, const advection_diffusion_full_3d_state *U_R,
  const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_lax_friedrichs_fluctuation *left_fluctuation)
{
  advection_diffusion_full_3d_lax_friedrichs_wave_family wave_family;
  advection_diffusion_full_3d_lax_friedrichs_y_wave_family(C_L, C_R, U_L, U_R, P, &wave_family);

  advection_diffusion_full_3d_lax_friedrichs_speed_family speed_family;
  advection_diffusion_full_3d_lax_friedrichs_y_speed_family(C_L, C_R, U_L, U_R, P, &speed_family);

  left_fluctuation->fluctuation_f = 
    (wave_family.wave1.wave_f * fmin(speed_family.speed1, 0.0)) + 
    (wave_family.wave2.wave_f * fmin(speed_family.speed2, 0.0));
}

static inline void
advection_diffusion_full_3d_lax_friedrichs_z_left_fluctuation(const advection_diffusion_full_3d_coordinates *C_L, const advection_diffusion_full_3d_coordinates *C_R,
  const advection_diffusion_full_3d_state *U_L, const advection_diffusion_full_3d_state *U_R,
  const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_lax_friedrichs_fluctuation *left_fluctuation)
{
  advection_diffusion_full_3d_lax_friedrichs_wave_family wave_family;
  advection_diffusion_full_3d_lax_friedrichs_z_wave_family(C_L, C_R, U_L, U_R, P, &wave_family);

  advection_diffusion_full_3d_lax_friedrichs_speed_family speed_family;
  advection_diffusion_full_3d_lax_friedrichs_z_speed_family(C_L, C_R, U_L, U_R, P, &speed_family);

  left_fluctuation->fluctuation_f = 
    (wave_family.wave1.wave_f * fmin(speed_family.speed1, 0.0)) + 
    (wave_family.wave2.wave_f * fmin(speed_family.speed2, 0.0));
}

static inline void
advection_diffusion_full_3d_lax_friedrichs_x_right_fluctuation(const advection_diffusion_full_3d_coordinates *C_L, const advection_diffusion_full_3d_coordinates *C_R,
  const advection_diffusion_full_3d_state *U_L, const advection_diffusion_full_3d_state *U_R,
  const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_lax_friedrichs_fluctuation *right_fluctuation)
{
  advection_diffusion_full_3d_lax_friedrichs_wave_family wave_family;
  advection_diffusion_full_3d_lax_friedrichs_x_wave_family(C_L, C_R, U_L, U_R, P, &wave_family);

  advection_diffusion_full_3d_lax_friedrichs_speed_family speed_family;
  advection_diffusion_full_3d_lax_friedrichs_x_speed_family(C_L, C_R, U_L, U_R, P, &speed_family);

  right_fluctuation->fluctuation_f = 
    (wave_family.wave1.wave_f * fmax(speed_family.speed1, 0.0)) + 
    (wave_family.wave2.wave_f * fmax(speed_family.speed2, 0.0));
}

static inline void
advection_diffusion_full_3d_lax_friedrichs_y_right_fluctuation(const advection_diffusion_full_3d_coordinates *C_L, const advection_diffusion_full_3d_coordinates *C_R,
  const advection_diffusion_full_3d_state *U_L, const advection_diffusion_full_3d_state *U_R,
  const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_lax_friedrichs_fluctuation *right_fluctuation)
{
  advection_diffusion_full_3d_lax_friedrichs_wave_family wave_family;
  advection_diffusion_full_3d_lax_friedrichs_y_wave_family(C_L, C_R, U_L, U_R, P, &wave_family);

  advection_diffusion_full_3d_lax_friedrichs_speed_family speed_family;
  advection_diffusion_full_3d_lax_friedrichs_y_speed_family(C_L, C_R, U_L, U_R, P, &speed_family);

  right_fluctuation->fluctuation_f = 
    (wave_family.wave1.wave_f * fmax(speed_family.speed1, 0.0)) + 
    (wave_family.wave2.wave_f * fmax(speed_family.speed2, 0.0));
}

static inline void
advection_diffusion_full_3d_lax_friedrichs_z_right_fluctuation(const advection_diffusion_full_3d_coordinates *C_L, const advection_diffusion_full_3d_coordinates *C_R,
  const advection_diffusion_full_3d_state *U_L, const advection_diffusion_full_3d_state *U_R,
  const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_lax_friedrichs_fluctuation *right_fluctuation)
{
  advection_diffusion_full_3d_lax_friedrichs_wave_family wave_family;
  advection_diffusion_full_3d_lax_friedrichs_z_wave_family(C_L, C_R, U_L, U_R, P, &wave_family);

  advection_diffusion_full_3d_lax_friedrichs_speed_family speed_family;
  advection_diffusion_full_3d_lax_friedrichs_z_speed_family(C_L, C_R, U_L, U_R, P, &speed_family);

  right_fluctuation->fluctuation_f = 
    (wave_family.wave1.wave_f * fmax(speed_family.speed1, 0.0)) + 
    (wave_family.wave2.wave_f * fmax(speed_family.speed2, 0.0));
}

static inline advection_diffusion_full_3d_lax_friedrichs_fluctuations_3D
create_left_fluctuations()
{
  advection_diffusion_full_3d_lax_friedrichs_fluctuations_3D left_fluctuations;

  left_fluctuations.x_fluctuation = advection_diffusion_full_3d_lax_friedrichs_x_left_fluctuation;
  left_fluctuations.y_fluctuation = advection_diffusion_full_3d_lax_friedrichs_y_left_fluctuation;
  left_fluctuations.z_fluctuation = advection_diffusion_full_3d_lax_friedrichs_z_left_fluctuation;

  return left_fluctuations;
}

static inline advection_diffusion_full_3d_lax_friedrichs_fluctuations_3D
create_right_fluctuations()
{
  advection_diffusion_full_3d_lax_friedrichs_fluctuations_3D right_fluctuations;

  right_fluctuations.x_fluctuation = advection_diffusion_full_3d_lax_friedrichs_x_right_fluctuation;
  right_fluctuations.y_fluctuation = advection_diffusion_full_3d_lax_friedrichs_y_right_fluctuation;
  right_fluctuations.z_fluctuation = advection_diffusion_full_3d_lax_friedrichs_z_right_fluctuation;

  return right_fluctuations;
}

static inline bool
advection_diffusion_full_3d_lax_friedrichs_x_left_fluctuations_consistent(const advection_diffusion_full_3d_coordinates *C, const advection_diffusion_full_3d_state *U,
  const advection_diffusion_full_3d_parameters *P)
{
  advection_diffusion_full_3d_lax_friedrichs_fluctuation left_fluctuation;
  advection_diffusion_full_3d_lax_friedrichs_x_left_fluctuation(C, C, U, U, P, &left_fluctuation);

  double f_left_fluctuation = left_fluctuation.fluctuation_f;

  return ((fabs(f_left_fluctuation) < 1.0e-8));
}

static inline bool
advection_diffusion_full_3d_lax_friedrichs_y_left_fluctuations_consistent(const advection_diffusion_full_3d_coordinates *C, const advection_diffusion_full_3d_state *U,
  const advection_diffusion_full_3d_parameters *P)
{
  advection_diffusion_full_3d_lax_friedrichs_fluctuation left_fluctuation;
  advection_diffusion_full_3d_lax_friedrichs_y_left_fluctuation(C, C, U, U, P, &left_fluctuation);

  double f_left_fluctuation = left_fluctuation.fluctuation_f;

  return ((fabs(f_left_fluctuation) < 1.0e-8));
}

static inline bool
advection_diffusion_full_3d_lax_friedrichs_z_left_fluctuations_consistent(const advection_diffusion_full_3d_coordinates *C, const advection_diffusion_full_3d_state *U,
  const advection_diffusion_full_3d_parameters *P)
{
  advection_diffusion_full_3d_lax_friedrichs_fluctuation left_fluctuation;
  advection_diffusion_full_3d_lax_friedrichs_z_left_fluctuation(C, C, U, U, P, &left_fluctuation);

  double f_left_fluctuation = left_fluctuation.fluctuation_f;

  return ((fabs(f_left_fluctuation) < 1.0e-8));
}

static inline bool
advection_diffusion_full_3d_lax_friedrichs_x_right_fluctuations_consistent(const advection_diffusion_full_3d_coordinates *C, const advection_diffusion_full_3d_state *U,
  const advection_diffusion_full_3d_parameters *P)
{
  advection_diffusion_full_3d_lax_friedrichs_fluctuation right_fluctuation;
  advection_diffusion_full_3d_lax_friedrichs_x_right_fluctuation(C, C, U, U, P, &right_fluctuation);

  double f_right_fluctuation = right_fluctuation.fluctuation_f;

  return ((fabs(f_right_fluctuation) < 1.0e-8));
}

static inline bool
advection_diffusion_full_3d_lax_friedrichs_y_right_fluctuations_consistent(const advection_diffusion_full_3d_coordinates *C, const advection_diffusion_full_3d_state *U,
  const advection_diffusion_full_3d_parameters *P)
{
  advection_diffusion_full_3d_lax_friedrichs_fluctuation right_fluctuation;
  advection_diffusion_full_3d_lax_friedrichs_y_right_fluctuation(C, C, U, U, P, &right_fluctuation);

  double f_right_fluctuation = right_fluctuation.fluctuation_f;

  return ((fabs(f_right_fluctuation) < 1.0e-8));
}

static inline bool
advection_diffusion_full_3d_lax_friedrichs_z_right_fluctuations_consistent(const advection_diffusion_full_3d_coordinates *C, const advection_diffusion_full_3d_state *U,
  const advection_diffusion_full_3d_parameters *P)
{
  advection_diffusion_full_3d_lax_friedrichs_fluctuation right_fluctuation;
  advection_diffusion_full_3d_lax_friedrichs_z_right_fluctuation(C, C, U, U, P, &right_fluctuation);

  double f_right_fluctuation = right_fluctuation.fluctuation_f;

  return ((fabs(f_right_fluctuation) < 1.0e-8));
}

static inline bool
advection_diffusion_full_3d_lax_friedrichs_x_fluctuations_valid(const advection_diffusion_full_3d_coordinates *C_L, const advection_diffusion_full_3d_coordinates *C_R,
  const advection_diffusion_full_3d_state *U_L, const advection_diffusion_full_3d_state *U_R, const advection_diffusion_full_3d_parameters *P)
{
  advection_diffusion_full_3d_flux left_flux;
  advection_diffusion_full_3d_x_flux(C_L, P, U_L, &left_flux);

  advection_diffusion_full_3d_flux right_flux;
  advection_diffusion_full_3d_x_flux(C_R, P, U_R, &right_flux);

  double f_flux_jump = right_flux.flux_f - left_flux.flux_f;

  advection_diffusion_full_3d_lax_friedrichs_fluctuation left_fluctuation;
  advection_diffusion_full_3d_lax_friedrichs_x_left_fluctuation(C_L, C_R, U_L, U_R, P, &left_fluctuation);

  advection_diffusion_full_3d_lax_friedrichs_fluctuation right_fluctuation;
  advection_diffusion_full_3d_lax_friedrichs_x_right_fluctuation(C_L, C_R, U_L, U_R, P, &right_fluctuation);

  double f_fluctuation_sum = left_fluctuation.fluctuation_f + right_fluctuation.fluctuation_f;

  return (fabs(f_flux_jump - f_fluctuation_sum) < 1.0e-8);
}

static inline bool
advection_diffusion_full_3d_lax_friedrichs_y_fluctuations_valid(const advection_diffusion_full_3d_coordinates *C_L, const advection_diffusion_full_3d_coordinates *C_R,
  const advection_diffusion_full_3d_state *U_L, const advection_diffusion_full_3d_state *U_R, const advection_diffusion_full_3d_parameters *P)
{
  advection_diffusion_full_3d_flux left_flux;
  advection_diffusion_full_3d_y_flux(C_L, P, U_L, &left_flux);

  advection_diffusion_full_3d_flux right_flux;
  advection_diffusion_full_3d_y_flux(C_R, P, U_R, &right_flux);

  double f_flux_jump = right_flux.flux_f - left_flux.flux_f;

  advection_diffusion_full_3d_lax_friedrichs_fluctuation left_fluctuation;
  advection_diffusion_full_3d_lax_friedrichs_y_left_fluctuation(C_L, C_R, U_L, U_R, P, &left_fluctuation);

  advection_diffusion_full_3d_lax_friedrichs_fluctuation right_fluctuation;
  advection_diffusion_full_3d_lax_friedrichs_y_right_fluctuation(C_L, C_R, U_L, U_R, P, &right_fluctuation);

  double f_fluctuation_sum = left_fluctuation.fluctuation_f + right_fluctuation.fluctuation_f;

  return (fabs(f_flux_jump - f_fluctuation_sum) < 1.0e-8);
}

static inline bool
advection_diffusion_full_3d_lax_friedrichs_z_fluctuations_valid(const advection_diffusion_full_3d_coordinates *C_L, const advection_diffusion_full_3d_coordinates *C_R,
  const advection_diffusion_full_3d_state *U_L, const advection_diffusion_full_3d_state *U_R, const advection_diffusion_full_3d_parameters *P)
{
  advection_diffusion_full_3d_flux left_flux;
  advection_diffusion_full_3d_z_flux(C_L, P, U_L, &left_flux);

  advection_diffusion_full_3d_flux right_flux;
  advection_diffusion_full_3d_z_flux(C_R, P, U_R, &right_flux);

  double f_flux_jump = right_flux.flux_f - left_flux.flux_f;

  advection_diffusion_full_3d_lax_friedrichs_fluctuation left_fluctuation;
  advection_diffusion_full_3d_lax_friedrichs_z_left_fluctuation(C_L, C_R, U_L, U_R, P, &left_fluctuation);

  advection_diffusion_full_3d_lax_friedrichs_fluctuation right_fluctuation;
  advection_diffusion_full_3d_lax_friedrichs_z_right_fluctuation(C_L, C_R, U_L, U_R, P, &right_fluctuation);

  double f_fluctuation_sum = left_fluctuation.fluctuation_f + right_fluctuation.fluctuation_f;

  return (fabs(f_flux_jump - f_fluctuation_sum) < 1.0e-8);
}

typedef struct {
  double reconstruction_f;
} advection_diffusion_full_3d_minmod_reconstruction;

typedef struct {
  void (*x_reconstruction)(const advection_diffusion_full_3d_coordinates *C_L, const advection_diffusion_full_3d_coordinates *C, const advection_diffusion_full_3d_coordinates *C_R,
    const advection_diffusion_full_3d_state *U_L, const advection_diffusion_full_3d_state *U, const advection_diffusion_full_3d_state *U_R,
    const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_minmod_reconstruction *reconstruction);
  void (*y_reconstruction)(const advection_diffusion_full_3d_coordinates *C_L, const advection_diffusion_full_3d_coordinates *C, const advection_diffusion_full_3d_coordinates *C_R,
    const advection_diffusion_full_3d_state *U_L, const advection_diffusion_full_3d_state *U, const advection_diffusion_full_3d_state *U_R,
    const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_minmod_reconstruction *reconstruction);
  void (*z_reconstruction)(const advection_diffusion_full_3d_coordinates *C_L, const advection_diffusion_full_3d_coordinates *C, const advection_diffusion_full_3d_coordinates *C_R,
    const advection_diffusion_full_3d_state *U_L, const advection_diffusion_full_3d_state *U, const advection_diffusion_full_3d_state *U_R,
    const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_minmod_reconstruction *reconstruction);
} advection_diffusion_full_3d_minmod_reconstruction_3D;

static inline void
advection_diffusion_full_3d_minmod_x_left_reconstruction(const advection_diffusion_full_3d_coordinates *C_L, const advection_diffusion_full_3d_coordinates *C,
  const advection_diffusion_full_3d_coordinates *C_R, const advection_diffusion_full_3d_state *U_L, const advection_diffusion_full_3d_state *U, const advection_diffusion_full_3d_state *U_R,
  const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_minmod_reconstruction *left_reconstruction)
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
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dxz = P->Dxz;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;
  double Dyz = P->Dyz;
  double Dzx = P->Dzx;
  double Dzy = P->Dzy;
  double Dzz = P->Dzz;

  left_reconstruction->reconstruction_f = (f - (0.5 * (fmax(0.0, fmin((f - f_L), (f_R - f))) + fmin(0.0, fmax((f - f_L),
    (f_R - f))))));
}

static inline void
advection_diffusion_full_3d_minmod_y_left_reconstruction(const advection_diffusion_full_3d_coordinates *C_L, const advection_diffusion_full_3d_coordinates *C,
  const advection_diffusion_full_3d_coordinates *C_R, const advection_diffusion_full_3d_state *U_L, const advection_diffusion_full_3d_state *U, const advection_diffusion_full_3d_state *U_R,
  const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_minmod_reconstruction *left_reconstruction)
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
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dxz = P->Dxz;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;
  double Dyz = P->Dyz;
  double Dzx = P->Dzx;
  double Dzy = P->Dzy;
  double Dzz = P->Dzz;

  left_reconstruction->reconstruction_f = (f - (0.5 * (fmax(0.0, fmin((f - f_L), (f_R - f))) + fmin(0.0, fmax((f - f_L),
    (f_R - f))))));
}

static inline void
advection_diffusion_full_3d_minmod_z_left_reconstruction(const advection_diffusion_full_3d_coordinates *C_L, const advection_diffusion_full_3d_coordinates *C,
  const advection_diffusion_full_3d_coordinates *C_R, const advection_diffusion_full_3d_state *U_L, const advection_diffusion_full_3d_state *U, const advection_diffusion_full_3d_state *U_R,
  const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_minmod_reconstruction *left_reconstruction)
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
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dxz = P->Dxz;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;
  double Dyz = P->Dyz;
  double Dzx = P->Dzx;
  double Dzy = P->Dzy;
  double Dzz = P->Dzz;

  left_reconstruction->reconstruction_f = (f - (0.5 * (fmax(0.0, fmin((f - f_L), (f_R - f))) + fmin(0.0, fmax((f - f_L),
    (f_R - f))))));
}

static inline void
advection_diffusion_full_3d_minmod_x_right_reconstruction(const advection_diffusion_full_3d_coordinates *C_L, const advection_diffusion_full_3d_coordinates *C,
  const advection_diffusion_full_3d_coordinates *C_R, const advection_diffusion_full_3d_state *U_L, const advection_diffusion_full_3d_state *U, const advection_diffusion_full_3d_state *U_R,
  const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_minmod_reconstruction *right_reconstruction)
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
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dxz = P->Dxz;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;
  double Dyz = P->Dyz;
  double Dzx = P->Dzx;
  double Dzy = P->Dzy;
  double Dzz = P->Dzz;

  right_reconstruction->reconstruction_f = (f + (0.5 * (fmax(0.0, fmin((f - f_L), (f_R - f))) + fmin(0.0, fmax((f -
    f_L), (f_R - f))))));
}

static inline void
advection_diffusion_full_3d_minmod_y_right_reconstruction(const advection_diffusion_full_3d_coordinates *C_L, const advection_diffusion_full_3d_coordinates *C,
  const advection_diffusion_full_3d_coordinates *C_R, const advection_diffusion_full_3d_state *U_L, const advection_diffusion_full_3d_state *U, const advection_diffusion_full_3d_state *U_R,
  const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_minmod_reconstruction *right_reconstruction)
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
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dxz = P->Dxz;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;
  double Dyz = P->Dyz;
  double Dzx = P->Dzx;
  double Dzy = P->Dzy;
  double Dzz = P->Dzz;

  right_reconstruction->reconstruction_f = (f + (0.5 * (fmax(0.0, fmin((f - f_L), (f_R - f))) + fmin(0.0, fmax((f -
    f_L), (f_R - f))))));
}

static inline void
advection_diffusion_full_3d_minmod_z_right_reconstruction(const advection_diffusion_full_3d_coordinates *C_L, const advection_diffusion_full_3d_coordinates *C,
  const advection_diffusion_full_3d_coordinates *C_R, const advection_diffusion_full_3d_state *U_L, const advection_diffusion_full_3d_state *U, const advection_diffusion_full_3d_state *U_R,
  const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_minmod_reconstruction *right_reconstruction)
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
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dxz = P->Dxz;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;
  double Dyz = P->Dyz;
  double Dzx = P->Dzx;
  double Dzy = P->Dzy;
  double Dzz = P->Dzz;

  right_reconstruction->reconstruction_f = (f + (0.5 * (fmax(0.0, fmin((f - f_L), (f_R - f))) + fmin(0.0, fmax((f -
    f_L), (f_R - f))))));
}

static inline advection_diffusion_full_3d_minmod_reconstruction_3D
create_left_reconstruction()
{
  advection_diffusion_full_3d_minmod_reconstruction_3D left_reconstruction;

  left_reconstruction.x_reconstruction = advection_diffusion_full_3d_minmod_x_left_reconstruction;
  left_reconstruction.y_reconstruction = advection_diffusion_full_3d_minmod_y_left_reconstruction;
  left_reconstruction.z_reconstruction = advection_diffusion_full_3d_minmod_z_left_reconstruction;

  return left_reconstruction;
}

static inline advection_diffusion_full_3d_minmod_reconstruction_3D
create_right_reconstruction()
{
  advection_diffusion_full_3d_minmod_reconstruction_3D right_reconstruction;

  right_reconstruction.x_reconstruction = advection_diffusion_full_3d_minmod_x_right_reconstruction;
  right_reconstruction.y_reconstruction = advection_diffusion_full_3d_minmod_y_right_reconstruction;
  right_reconstruction.z_reconstruction = advection_diffusion_full_3d_minmod_z_right_reconstruction;

  return right_reconstruction;
}

static inline bool
advection_diffusion_full_3d_minmod_x_left_reconstruction_consistent(const advection_diffusion_full_3d_coordinates *C, const advection_diffusion_full_3d_state *U,
const advection_diffusion_full_3d_parameters *P)
{
  advection_diffusion_full_3d_minmod_reconstruction left_reconstruction;
  advection_diffusion_full_3d_minmod_x_left_reconstruction(C, C, C, U, U, U, P, &left_reconstruction);

  double f_left_reconstruction = left_reconstruction.reconstruction_f;
  double f = U->f;

  return ((fabs(f_left_reconstruction - f) < 1.0e-8));
}

static inline bool
advection_diffusion_full_3d_minmod_y_left_reconstruction_consistent(const advection_diffusion_full_3d_coordinates *C, const advection_diffusion_full_3d_state *U,
const advection_diffusion_full_3d_parameters *P)
{
  advection_diffusion_full_3d_minmod_reconstruction left_reconstruction;
  advection_diffusion_full_3d_minmod_y_left_reconstruction(C, C, C, U, U, U, P, &left_reconstruction);

  double f_left_reconstruction = left_reconstruction.reconstruction_f;
  double f = U->f;

  return ((fabs(f_left_reconstruction - f) < 1.0e-8));
}

static inline bool
advection_diffusion_full_3d_minmod_z_left_reconstruction_consistent(const advection_diffusion_full_3d_coordinates *C, const advection_diffusion_full_3d_state *U,
const advection_diffusion_full_3d_parameters *P)
{
  advection_diffusion_full_3d_minmod_reconstruction left_reconstruction;
  advection_diffusion_full_3d_minmod_z_left_reconstruction(C, C, C, U, U, U, P, &left_reconstruction);

  double f_left_reconstruction = left_reconstruction.reconstruction_f;
  double f = U->f;

  return ((fabs(f_left_reconstruction - f) < 1.0e-8));
}

static inline bool
advection_diffusion_full_3d_minmod_x_right_reconstruction_consistent(const advection_diffusion_full_3d_coordinates *C, const advection_diffusion_full_3d_state *U,
const advection_diffusion_full_3d_parameters *P)
{
  advection_diffusion_full_3d_minmod_reconstruction right_reconstruction;
  advection_diffusion_full_3d_minmod_x_right_reconstruction(C, C, C, U, U, U, P, &right_reconstruction);

  double f_right_reconstruction = right_reconstruction.reconstruction_f;
  double f = U->f;

  return ((fabs(f_right_reconstruction - f) < 1.0e-8));
}

static inline bool
advection_diffusion_full_3d_minmod_y_right_reconstruction_consistent(const advection_diffusion_full_3d_coordinates *C, const advection_diffusion_full_3d_state *U,
const advection_diffusion_full_3d_parameters *P)
{
  advection_diffusion_full_3d_minmod_reconstruction right_reconstruction;
  advection_diffusion_full_3d_minmod_y_right_reconstruction(C, C, C, U, U, U, P, &right_reconstruction);

  double f_right_reconstruction = right_reconstruction.reconstruction_f;
  double f = U->f;

  return ((fabs(f_right_reconstruction - f) < 1.0e-8));
}

static inline bool
advection_diffusion_full_3d_minmod_z_right_reconstruction_consistent(const advection_diffusion_full_3d_coordinates *C, const advection_diffusion_full_3d_state *U,
const advection_diffusion_full_3d_parameters *P)
{
  advection_diffusion_full_3d_minmod_reconstruction right_reconstruction;
  advection_diffusion_full_3d_minmod_z_right_reconstruction(C, C, C, U, U, U, P, &right_reconstruction);

  double f_right_reconstruction = right_reconstruction.reconstruction_f;
  double f = U->f;

  return ((fabs(f_right_reconstruction - f) < 1.0e-8));
}

typedef struct {
  double interface_f;
} advection_diffusion_full_3d_centered_diffusion_interface_state;

typedef struct {
  void (*x_interface_state)(const advection_diffusion_full_3d_coordinates *C_Lx, const advection_diffusion_full_3d_coordinates *C_Rx,
    const advection_diffusion_full_3d_state *U_Lx, const advection_diffusion_full_3d_state *U_Rx,
    const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_centered_diffusion_interface_state *interface_state);
  void (*y_interface_state)(const advection_diffusion_full_3d_coordinates *C_Ly, const advection_diffusion_full_3d_coordinates *C_Ry,
    const advection_diffusion_full_3d_state *U_Ly, const advection_diffusion_full_3d_state *U_Ry,
    const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_centered_diffusion_interface_state *interface_state);
  void (*z_interface_state)(const advection_diffusion_full_3d_coordinates *C_Lz, const advection_diffusion_full_3d_coordinates *C_Rz,
    const advection_diffusion_full_3d_state *U_Lz, const advection_diffusion_full_3d_state *U_Rz,
    const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_centered_diffusion_interface_state *interface_state);
} advection_diffusion_full_3d_centered_diffusion_interface_states_3D;

static inline void
advection_diffusion_full_3d_centered_diffusion_x_interface_state(const advection_diffusion_full_3d_coordinates *C_Lx, const advection_diffusion_full_3d_coordinates *C_Rx,
  const advection_diffusion_full_3d_state *U_Lx, const advection_diffusion_full_3d_state *U_Rx,
  const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_centered_diffusion_interface_state *interface_state)
{
  double x_Lx = C_Lx->x;
  double x_Rx = C_Rx->x;
  double y_Lx = C_Lx->y;
  double y_Rx = C_Rx->y;
  double z_Lx = C_Lx->z;
  double z_Rx = C_Rx->z;
  double f_Lx = U_Lx->f;
  double f_Rx = U_Rx->f;
  double ax = P->ax;
  double ay = P->ay;
  double az = P->az;
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dxz = P->Dxz;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;
  double Dyz = P->Dyz;
  double Dzx = P->Dzx;
  double Dzy = P->Dzy;
  double Dzz = P->Dzz;

  interface_state->interface_f = (0.5 * (f_Lx + f_Rx));
}

static inline void
advection_diffusion_full_3d_centered_diffusion_y_interface_state(const advection_diffusion_full_3d_coordinates *C_Ly, const advection_diffusion_full_3d_coordinates *C_Ry,
  const advection_diffusion_full_3d_state *U_Ly, const advection_diffusion_full_3d_state *U_Ry,
  const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_centered_diffusion_interface_state *interface_state)
{
  double x_Ly = C_Ly->x;
  double x_Ry = C_Ry->x;
  double y_Ly = C_Ly->y;
  double y_Ry = C_Ry->y;
  double z_Ly = C_Ly->z;
  double z_Ry = C_Ry->z;
  double f_Ly = U_Ly->f;
  double f_Ry = U_Ry->f;
  double ax = P->ax;
  double ay = P->ay;
  double az = P->az;
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dxz = P->Dxz;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;
  double Dyz = P->Dyz;
  double Dzx = P->Dzx;
  double Dzy = P->Dzy;
  double Dzz = P->Dzz;

  interface_state->interface_f = (0.5 * (f_Ly + f_Ry));
}

static inline void
advection_diffusion_full_3d_centered_diffusion_z_interface_state(const advection_diffusion_full_3d_coordinates *C_Lz, const advection_diffusion_full_3d_coordinates *C_Rz,
  const advection_diffusion_full_3d_state *U_Lz, const advection_diffusion_full_3d_state *U_Rz,
  const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_centered_diffusion_interface_state *interface_state)
{
  double x_Lz = C_Lz->x;
  double x_Rz = C_Rz->x;
  double y_Lz = C_Lz->y;
  double y_Rz = C_Rz->y;
  double z_Lz = C_Lz->z;
  double z_Rz = C_Rz->z;
  double f_Lz = U_Lz->f;
  double f_Rz = U_Rz->f;
  double ax = P->ax;
  double ay = P->ay;
  double az = P->az;
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dxz = P->Dxz;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;
  double Dyz = P->Dyz;
  double Dzx = P->Dzx;
  double Dzy = P->Dzy;
  double Dzz = P->Dzz;

  interface_state->interface_f = (0.5 * (f_Lz + f_Rz));
}

static inline advection_diffusion_full_3d_centered_diffusion_interface_states_3D
create_interface_states()
{
  advection_diffusion_full_3d_centered_diffusion_interface_states_3D interface_states;

  interface_states.x_interface_state = advection_diffusion_full_3d_centered_diffusion_x_interface_state;
  interface_states.y_interface_state = advection_diffusion_full_3d_centered_diffusion_y_interface_state;
  interface_states.z_interface_state = advection_diffusion_full_3d_centered_diffusion_z_interface_state;
  
  return interface_states;
}

typedef struct {
  double interface_f_x;
  double interface_f_y;
  double interface_f_z;
} advection_diffusion_full_3d_centered_diffusion_interface_gradient;

typedef struct {
  void (*x_interface_gradient)(const advection_diffusion_full_3d_coordinates *C_Lx, const advection_diffusion_full_3d_coordinates *C_Rx,
    const advection_diffusion_full_3d_state *U_Lx, const advection_diffusion_full_3d_state *U_Rx,
    const advection_diffusion_full_3d_coordinates *C_Lx_Ly, const advection_diffusion_full_3d_coordinates *C_Lx_Ry,
    const advection_diffusion_full_3d_coordinates *C_Rx_Ly, const advection_diffusion_full_3d_coordinates *C_Rx_Ry,
    const advection_diffusion_full_3d_state *U_Lx_Ly, const advection_diffusion_full_3d_state *U_Lx_Ry,
    const advection_diffusion_full_3d_state *U_Rx_Ly, const advection_diffusion_full_3d_state *U_Rx_Ry,
    const advection_diffusion_full_3d_coordinates *C_Lx_Lz, const advection_diffusion_full_3d_coordinates *C_Lx_Rz,
    const advection_diffusion_full_3d_coordinates *C_Rx_Lz, const advection_diffusion_full_3d_coordinates *C_Rx_Rz,
    const advection_diffusion_full_3d_state *U_Lx_Lz, const advection_diffusion_full_3d_state *U_Lx_Rz,
    const advection_diffusion_full_3d_state *U_Rx_Lz, const advection_diffusion_full_3d_state *U_Rx_Rz,
    const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_centered_diffusion_interface_gradient *interface_gradient);
  void (*y_interface_gradient)(const advection_diffusion_full_3d_coordinates *C_Ly_Lx, const advection_diffusion_full_3d_coordinates *C_Ly_Rx,
    const advection_diffusion_full_3d_coordinates *C_Ry_Lx, const advection_diffusion_full_3d_coordinates *C_Ry_Rx,
    const advection_diffusion_full_3d_state *U_Ly_Lx, const advection_diffusion_full_3d_state *U_Ly_Rx,
    const advection_diffusion_full_3d_state *U_Ry_Lx, const advection_diffusion_full_3d_state *U_Ry_Rx,
    const advection_diffusion_full_3d_coordinates *C_Ly, const advection_diffusion_full_3d_coordinates *C_Ry,
    const advection_diffusion_full_3d_state *U_Ly, const advection_diffusion_full_3d_state *U_Ry,
    const advection_diffusion_full_3d_coordinates *C_Ly_Lz, const advection_diffusion_full_3d_coordinates *C_Ly_Rz,
    const advection_diffusion_full_3d_coordinates *C_Ry_Lz, const advection_diffusion_full_3d_coordinates *C_Ry_Rz,
    const advection_diffusion_full_3d_state *U_Ly_Lz, const advection_diffusion_full_3d_state *U_Ly_Rz,
    const advection_diffusion_full_3d_state *U_Ry_Lz, const advection_diffusion_full_3d_state *U_Ry_Rz,
    const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_centered_diffusion_interface_gradient *interface_gradient);
  void (*z_interface_gradient)(const advection_diffusion_full_3d_coordinates *C_Lz_Lx, const advection_diffusion_full_3d_coordinates *C_Lz_Rx,
    const advection_diffusion_full_3d_coordinates *C_Rz_Lx, const advection_diffusion_full_3d_coordinates *C_Rz_Rx,
    const advection_diffusion_full_3d_state *U_Lz_Lx, const advection_diffusion_full_3d_state *U_Lz_Rx,
    const advection_diffusion_full_3d_state *U_Rz_Lx, const advection_diffusion_full_3d_state *U_Rz_Rx,
    const advection_diffusion_full_3d_coordinates *C_Lz_Ly, const advection_diffusion_full_3d_coordinates *C_Lz_Ry,
    const advection_diffusion_full_3d_coordinates *C_Rz_Ly, const advection_diffusion_full_3d_coordinates *C_Rz_Ry,
    const advection_diffusion_full_3d_state *U_Lz_Ly, const advection_diffusion_full_3d_state *U_Lz_Ry,
    const advection_diffusion_full_3d_state *U_Rz_Ly, const advection_diffusion_full_3d_state *U_Rz_Ry,
    const advection_diffusion_full_3d_coordinates *C_Lz, const advection_diffusion_full_3d_coordinates *C_Rz,
    const advection_diffusion_full_3d_state *U_Lz, const advection_diffusion_full_3d_state *U_Rz,
    const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_centered_diffusion_interface_gradient *interface_gradient);
} advection_diffusion_full_3d_centered_diffusion_interface_gradients_3D;

static inline void
advection_diffusion_full_3d_centered_diffusion_x_interface_gradient(const advection_diffusion_full_3d_coordinates *C_Lx, const advection_diffusion_full_3d_coordinates *C_Rx,
  const advection_diffusion_full_3d_state *U_Lx, const advection_diffusion_full_3d_state *U_Rx,
  const advection_diffusion_full_3d_coordinates *C_Lx_Ly, const advection_diffusion_full_3d_coordinates *C_Lx_Ry,
  const advection_diffusion_full_3d_coordinates *C_Rx_Ly, const advection_diffusion_full_3d_coordinates *C_Rx_Ry,
  const advection_diffusion_full_3d_state *U_Lx_Ly, const advection_diffusion_full_3d_state *U_Lx_Ry,
  const advection_diffusion_full_3d_state *U_Rx_Ly, const advection_diffusion_full_3d_state *U_Rx_Ry,
  const advection_diffusion_full_3d_coordinates *C_Lx_Lz, const advection_diffusion_full_3d_coordinates *C_Lx_Rz,
  const advection_diffusion_full_3d_coordinates *C_Rx_Lz, const advection_diffusion_full_3d_coordinates *C_Rx_Rz,
  const advection_diffusion_full_3d_state *U_Lx_Lz, const advection_diffusion_full_3d_state *U_Lx_Rz,
  const advection_diffusion_full_3d_state *U_Rx_Lz, const advection_diffusion_full_3d_state *U_Rx_Rz,
  const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_centered_diffusion_interface_gradient *interface_gradient)
{
  double x_Lx = C_Lx->x;
  double x_Rx = C_Rx->x;
  double y_Lx_Ly = C_Lx_Ly->y;
  double y_Lx_Ry = C_Lx_Ry->y;
  double y_Rx_Ly = C_Rx_Ly->y;
  double y_Rx_Ry = C_Rx_Ry->y;
  double z_Lx_Lz = C_Lx_Lz->z;
  double z_Lx_Rz = C_Lx_Rz->z;
  double z_Rx_Lz = C_Rx_Lz->z;
  double z_Rx_Rz = C_Rx_Rz->z;

  double f_Lx = U_Lx->f;
  double f_Rx = U_Rx->f;
  double f_Lx_Ly = U_Lx_Ly->f;
  double f_Lx_Ry = U_Lx_Ry->f;
  double f_Rx_Ly = U_Rx_Ly->f;
  double f_Rx_Ry = U_Rx_Ry->f;
  double f_Lx_Lz = U_Lx_Lz->f;
  double f_Lx_Rz = U_Lx_Rz->f;
  double f_Rx_Lz = U_Rx_Lz->f;
  double f_Rx_Rz = U_Rx_Rz->f;

  double ax = P->ax;
  double ay = P->ay;
  double az = P->az;
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dxz = P->Dxz;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;
  double Dyz = P->Dyz;
  double Dzx = P->Dzx;
  double Dzy = P->Dzy;
  double Dzz = P->Dzz;

  interface_gradient->interface_f_x = ((f_Rx - f_Lx) / (x_Rx - x_Lx));

  interface_gradient->interface_f_y = (0.5 * (((f_Lx_Ry - f_Lx_Ly) / (y_Lx_Ry - y_Lx_Ly)) + ((f_Rx_Ry - f_Rx_Ly) /
    (y_Rx_Ry - y_Rx_Ly))));

  interface_gradient->interface_f_z = (0.5 * (((f_Lx_Rz - f_Lx_Lz) / (z_Lx_Rz - z_Lx_Lz)) + ((f_Rx_Rz - f_Rx_Lz) /
    (z_Rx_Rz - z_Rx_Lz))));
}

static inline void
advection_diffusion_full_3d_centered_diffusion_y_interface_gradient(const advection_diffusion_full_3d_coordinates *C_Ly_Lx, const advection_diffusion_full_3d_coordinates *C_Ly_Rx,
  const advection_diffusion_full_3d_coordinates *C_Ry_Lx, const advection_diffusion_full_3d_coordinates *C_Ry_Rx,
  const advection_diffusion_full_3d_state *U_Ly_Lx, const advection_diffusion_full_3d_state *U_Ly_Rx,
  const advection_diffusion_full_3d_state *U_Ry_Lx, const advection_diffusion_full_3d_state *U_Ry_Rx,
  const advection_diffusion_full_3d_coordinates *C_Ly, const advection_diffusion_full_3d_coordinates *C_Ry,
  const advection_diffusion_full_3d_state *U_Ly, const advection_diffusion_full_3d_state *U_Ry,
  const advection_diffusion_full_3d_coordinates *C_Ly_Lz, const advection_diffusion_full_3d_coordinates *C_Ly_Rz,
  const advection_diffusion_full_3d_coordinates *C_Ry_Lz, const advection_diffusion_full_3d_coordinates *C_Ry_Rz,
  const advection_diffusion_full_3d_state *U_Ly_Lz, const advection_diffusion_full_3d_state *U_Ly_Rz,
  const advection_diffusion_full_3d_state *U_Ry_Lz, const advection_diffusion_full_3d_state *U_Ry_Rz,
  const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_centered_diffusion_interface_gradient *interface_gradient)
{
  double x_Ly_Lx = C_Ly_Lx->x;
  double x_Ly_Rx = C_Ly_Rx->x;
  double x_Ry_Lx = C_Ry_Lx->x;
  double x_Ry_Rx = C_Ry_Rx->x;
  double y_Ly = C_Ly->y;
  double y_Ry = C_Ry->y;
  double z_Ly_Lz = C_Ly_Lz->z;
  double z_Ly_Rz = C_Ly_Rz->z;
  double z_Ry_Lz = C_Ry_Lz->z;
  double z_Ry_Rz = C_Ry_Rz->z;

  double f_Ly_Lx = U_Ly_Lx->f;
  double f_Ly_Rx = U_Ly_Rx->f;
  double f_Ry_Lx = U_Ry_Lx->f;
  double f_Ry_Rx = U_Ry_Rx->f;
  double f_Ly = U_Ly->f;
  double f_Ry = U_Ry->f;
  double f_Ly_Lz = U_Ly_Lz->f;
  double f_Ly_Rz = U_Ly_Rz->f;
  double f_Ry_Lz = U_Ry_Lz->f;
  double f_Ry_Rz = U_Ry_Rz->f;

  double ax = P->ax;
  double ay = P->ay;
  double az = P->az;
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dxz = P->Dxz;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;
  double Dyz = P->Dyz;
  double Dzx = P->Dzx;
  double Dzy = P->Dzy;
  double Dzz = P->Dzz;

  interface_gradient->interface_f_x = (0.5 * (((f_Ly_Rx - f_Ly_Lx) / (x_Ly_Rx - x_Ly_Lx)) + ((f_Ry_Rx - f_Ry_Lx) /
    (x_Ry_Rx - x_Ry_Lx))));

  interface_gradient->interface_f_y = ((f_Ry - f_Ly) / (y_Ry - y_Ly));

  interface_gradient->interface_f_z = (0.5 * (((f_Ly_Rz - f_Ly_Lz) / (z_Ly_Rz - z_Ly_Lz)) + ((f_Ry_Rz - f_Ry_Lz) /
    (z_Ry_Rz - z_Ry_Lz))));
}

static inline void
advection_diffusion_full_3d_centered_diffusion_z_interface_gradient(const advection_diffusion_full_3d_coordinates *C_Lz_Lx, const advection_diffusion_full_3d_coordinates *C_Lz_Rx,
  const advection_diffusion_full_3d_coordinates *C_Rz_Lx, const advection_diffusion_full_3d_coordinates *C_Rz_Rx,
  const advection_diffusion_full_3d_state *U_Lz_Lx, const advection_diffusion_full_3d_state *U_Lz_Rx,
  const advection_diffusion_full_3d_state *U_Rz_Lx, const advection_diffusion_full_3d_state *U_Rz_Rx,
  const advection_diffusion_full_3d_coordinates *C_Lz_Ly, const advection_diffusion_full_3d_coordinates *C_Lz_Ry,
  const advection_diffusion_full_3d_coordinates *C_Rz_Ly, const advection_diffusion_full_3d_coordinates *C_Rz_Ry,
  const advection_diffusion_full_3d_state *U_Lz_Ly, const advection_diffusion_full_3d_state *U_Lz_Ry,
  const advection_diffusion_full_3d_state *U_Rz_Ly, const advection_diffusion_full_3d_state *U_Rz_Ry,
  const advection_diffusion_full_3d_coordinates *C_Lz, const advection_diffusion_full_3d_coordinates *C_Rz,
  const advection_diffusion_full_3d_state *U_Lz, const advection_diffusion_full_3d_state *U_Rz,
  const advection_diffusion_full_3d_parameters *P, advection_diffusion_full_3d_centered_diffusion_interface_gradient *interface_gradient)
{
  double x_Lz_Lx = C_Lz_Lx->x;
  double x_Lz_Rx = C_Lz_Rx->x;
  double x_Rz_Lx = C_Rz_Lx->x;
  double x_Rz_Rx = C_Rz_Rx->x;
  double y_Lz_Ly = C_Lz_Ly->y;
  double y_Lz_Ry = C_Lz_Ry->y;
  double y_Rz_Ly = C_Rz_Ly->y;
  double y_Rz_Ry = C_Rz_Ry->y;
  double z_Lz = C_Lz->z;
  double z_Rz = C_Rz->z;

  double f_Lz_Lx = U_Lz_Lx->f;
  double f_Lz_Rx = U_Lz_Rx->f;
  double f_Rz_Lx = U_Rz_Lx->f;
  double f_Rz_Rx = U_Rz_Rx->f;
  double f_Lz_Ly = U_Lz_Ly->f;
  double f_Lz_Ry = U_Lz_Ry->f;
  double f_Rz_Ly = U_Rz_Ly->f;
  double f_Rz_Ry = U_Rz_Ry->f;
  double f_Lz = U_Lz->f;
  double f_Rz = U_Rz->f;

  double ax = P->ax;
  double ay = P->ay;
  double az = P->az;
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dxz = P->Dxz;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;
  double Dyz = P->Dyz;
  double Dzx = P->Dzx;
  double Dzy = P->Dzy;
  double Dzz = P->Dzz;

  interface_gradient->interface_f_x = (0.5 * (((f_Lz_Rx - f_Lz_Lx) / (x_Lz_Rx - x_Lz_Lx)) + ((f_Rz_Rx - f_Rz_Lx) /
    (x_Rz_Rx - x_Rz_Lx))));

  interface_gradient->interface_f_y = (0.5 * (((f_Lz_Ry - f_Lz_Ly) / (y_Lz_Ry - y_Lz_Ly)) + ((f_Rz_Ry - f_Rz_Ly) /
    (y_Rz_Ry - y_Rz_Ly))));

  interface_gradient->interface_f_z = ((f_Rz - f_Lz) / (z_Rz - z_Lz));
}

static inline advection_diffusion_full_3d_centered_diffusion_interface_gradients_3D
create_interface_gradients()
{
  advection_diffusion_full_3d_centered_diffusion_interface_gradients_3D interface_gradients;

  interface_gradients.x_interface_gradient = advection_diffusion_full_3d_centered_diffusion_x_interface_gradient;
  interface_gradients.y_interface_gradient = advection_diffusion_full_3d_centered_diffusion_y_interface_gradient;
  interface_gradients.z_interface_gradient = advection_diffusion_full_3d_centered_diffusion_z_interface_gradient;
  
  return interface_gradients;
}

int
main(int argc, char *argv[])
{
  // Insert simulation drivers here.
  
  return 0;
}
                