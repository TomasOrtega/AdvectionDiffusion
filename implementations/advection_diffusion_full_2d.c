#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

typedef struct {
  double x;
  double y;
} advection_diffusion_full_2d_coordinates;

typedef struct {
  double f;
} advection_diffusion_full_2d_state;

static inline bool
advection_diffusion_full_2d_state_valid(const advection_diffusion_full_2d_state *U)
{
  double f = U->f;
  
  return true;
}

typedef struct {
  double ax;
  double ay;
  double Dxx;
  double Dxy;
  double Dyx;
  double Dyy;
} advection_diffusion_full_2d_parameters;

static inline bool
advection_diffusion_full_2d_parameters_valid(const advection_diffusion_full_2d_parameters *P)
{
  double ax = P->ax;
  double ay = P->ay;
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;

  return (Dxx > 0.0 && Dyy > 0.0);
}

typedef struct {
  double flux_f;
} advection_diffusion_full_2d_flux;

typedef struct {
  void (*x_flux)(const advection_diffusion_full_2d_coordinates *C, const advection_diffusion_full_2d_parameters *P,
    const advection_diffusion_full_2d_state *U, advection_diffusion_full_2d_flux *flux);
  void (*y_flux)(const advection_diffusion_full_2d_coordinates *C, const advection_diffusion_full_2d_parameters *P,
    const advection_diffusion_full_2d_state *U, advection_diffusion_full_2d_flux *flux);
} advection_diffusion_full_2d_fluxes_2D;

static inline void
advection_diffusion_full_2d_x_flux(const advection_diffusion_full_2d_coordinates *C, const advection_diffusion_full_2d_parameters *P,
  const advection_diffusion_full_2d_state *U, advection_diffusion_full_2d_flux *flux)
{
  double x = C->x;
  double y = C->y;
  double ax = P->ax;
  double ay = P->ay;
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;
  double f = U->f;

  flux->flux_f = (ax * f);
}

static inline void
advection_diffusion_full_2d_y_flux(const advection_diffusion_full_2d_coordinates *C, const advection_diffusion_full_2d_parameters *P,
  const advection_diffusion_full_2d_state *U, advection_diffusion_full_2d_flux *flux)
{
  double x = C->x;
  double y = C->y;
  double ax = P->ax;
  double ay = P->ay;
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;
  double f = U->f;

  flux->flux_f = (ay * f);
}

static inline advection_diffusion_full_2d_fluxes_2D
create_fluxes()
{
  advection_diffusion_full_2d_fluxes_2D fluxes;

  fluxes.x_flux = advection_diffusion_full_2d_x_flux;
  fluxes.y_flux = advection_diffusion_full_2d_y_flux;
  
  return fluxes;
}

typedef struct {
  double mu1;
} advection_diffusion_full_2d_wavespeed;

typedef struct {
  void (*x_wavespeed)(const advection_diffusion_full_2d_coordinates *C, const advection_diffusion_full_2d_parameters *P,
    const advection_diffusion_full_2d_state *U, advection_diffusion_full_2d_wavespeed *wavespeed);
  void (*y_wavespeed)(const advection_diffusion_full_2d_coordinates *C, const advection_diffusion_full_2d_parameters *P,
    const advection_diffusion_full_2d_state *U, advection_diffusion_full_2d_wavespeed *wavespeed);
} advection_diffusion_full_2d_wavespeed_estimate_2D;

static inline void
advection_diffusion_full_2d_x_wavespeed(const advection_diffusion_full_2d_coordinates *C, const advection_diffusion_full_2d_parameters *P,
  const advection_diffusion_full_2d_state *U, advection_diffusion_full_2d_wavespeed *wavespeed)
{
  double x = C->x;
  double y = C->y;
  double ax = P->ax;
  double ay = P->ay;
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;
  double f = U->f;

  wavespeed->mu1 = ax;
}

static inline void
advection_diffusion_full_2d_y_wavespeed(const advection_diffusion_full_2d_coordinates *C, const advection_diffusion_full_2d_parameters *P,
  const advection_diffusion_full_2d_state *U, advection_diffusion_full_2d_wavespeed *wavespeed)
{
  double x = C->x;
  double y = C->y;
  double ax = P->ax;
  double ay = P->ay;
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;
  double f = U->f;

  wavespeed->mu1 = ay;
}

static inline advection_diffusion_full_2d_wavespeed_estimate_2D
create_wavespeeds()
{
  advection_diffusion_full_2d_wavespeed_estimate_2D wavespeed_estimate;

  wavespeed_estimate.x_wavespeed = advection_diffusion_full_2d_x_wavespeed;
  wavespeed_estimate.y_wavespeed = advection_diffusion_full_2d_y_wavespeed;
  
  return wavespeed_estimate;
}

typedef struct {
  double f_x;
  double f_y;
} advection_diffusion_full_2d_gradient;

typedef struct {
  double diffusive_flux_f;
} advection_diffusion_full_2d_diffusive_flux;

typedef struct {
  void (*x_diffusive_flux)(const advection_diffusion_full_2d_coordinates *C, const advection_diffusion_full_2d_parameters *P,
    const advection_diffusion_full_2d_state *U, advection_diffusion_full_2d_gradient *DU, advection_diffusion_full_2d_diffusive_flux *flux);
  void (*y_diffusive_flux)(const advection_diffusion_full_2d_coordinates *C, const advection_diffusion_full_2d_parameters *P,
    const advection_diffusion_full_2d_state *U, advection_diffusion_full_2d_gradient *DU, advection_diffusion_full_2d_diffusive_flux *flux);
} advection_diffusion_full_2d_diffusive_fluxes_2D;

static inline void
advection_diffusion_full_2d_x_diffusive_flux(const advection_diffusion_full_2d_coordinates *C, const advection_diffusion_full_2d_parameters *P,
  const advection_diffusion_full_2d_state *U, advection_diffusion_full_2d_gradient *DU, advection_diffusion_full_2d_diffusive_flux *diffusive_flux)
{
  double x = C->x;
  double y = C->y;
  double ax = P->ax;
  double ay = P->ay;
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;
  double f = U->f;
  double f_x = DU->f_x;
  double f_y = DU->f_y;

  diffusive_flux->diffusive_flux_f = ((Dxx * f_x) + (Dxy * f_y));
}

static inline void
advection_diffusion_full_2d_y_diffusive_flux(const advection_diffusion_full_2d_coordinates *C, const advection_diffusion_full_2d_parameters *P,
  const advection_diffusion_full_2d_state *U, advection_diffusion_full_2d_gradient *DU, advection_diffusion_full_2d_diffusive_flux *diffusive_flux)
{
  double x = C->x;
  double y = C->y;
  double ax = P->ax;
  double ay = P->ay;
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;
  double f = U->f;
  double f_x = DU->f_x;
  double f_y = DU->f_y;

  diffusive_flux->diffusive_flux_f = ((Dyx * f_x) + (Dyy * f_y));
}

static inline advection_diffusion_full_2d_diffusive_fluxes_2D
create_diffusive_fluxes()
{
  advection_diffusion_full_2d_diffusive_fluxes_2D diffusive_fluxes;

  diffusive_fluxes.x_diffusive_flux = advection_diffusion_full_2d_x_diffusive_flux;
  diffusive_fluxes.y_diffusive_flux = advection_diffusion_full_2d_y_diffusive_flux;
  
  return diffusive_fluxes;
}

typedef struct {
  double wave_f;
} advection_diffusion_full_2d_lax_friedrichs_wave;

typedef struct {
  advection_diffusion_full_2d_lax_friedrichs_wave wave1;
  advection_diffusion_full_2d_lax_friedrichs_wave wave2;
} advection_diffusion_full_2d_lax_friedrichs_wave_family;

typedef struct {
  void (*x_wave_family)(const advection_diffusion_full_2d_coordinates *C_L, const advection_diffusion_full_2d_coordinates *C_R,
    const advection_diffusion_full_2d_state *U_L, const advection_diffusion_full_2d_state *U_R,
    const advection_diffusion_full_2d_parameters *P, advection_diffusion_full_2d_lax_friedrichs_wave_family *wave_family);
  void (*y_wave_family)(const advection_diffusion_full_2d_coordinates *C_L, const advection_diffusion_full_2d_coordinates *C_R,
    const advection_diffusion_full_2d_state *U_L, const advection_diffusion_full_2d_state *U_R,
    const advection_diffusion_full_2d_parameters *P, advection_diffusion_full_2d_lax_friedrichs_wave_family *wave_family);
} advection_diffusion_full_2d_lax_friedrichs_waves_2D;

static inline void
advection_diffusion_full_2d_lax_friedrichs_x_wave_family(const advection_diffusion_full_2d_coordinates *C_L, const advection_diffusion_full_2d_coordinates *C_R,
  const advection_diffusion_full_2d_state *U_L, const advection_diffusion_full_2d_state *U_R,
  const advection_diffusion_full_2d_parameters *P, advection_diffusion_full_2d_lax_friedrichs_wave_family *wave_family)
{
  double x_L = C_L->x;
  double x_R = C_R->x;
  double y_L = C_L->y;
  double y_R = C_R->y;
  double f_L = U_L->f;
  double f_R = U_R->f;
  double ax = P->ax;
  double ay = P->ay;
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;

  wave_family->wave1.wave_f = (0.5 * ((f_R - f_L) - (((ax * f_R) - (ax * f_L)) / fabs(ax))));

  wave_family->wave2.wave_f = (0.5 * ((f_R - f_L) + (((ax * f_R) - (ax * f_L)) / fabs(ax))));
}

static inline void
advection_diffusion_full_2d_lax_friedrichs_y_wave_family(const advection_diffusion_full_2d_coordinates *C_L, const advection_diffusion_full_2d_coordinates *C_R,
  const advection_diffusion_full_2d_state *U_L, const advection_diffusion_full_2d_state *U_R,
  const advection_diffusion_full_2d_parameters *P, advection_diffusion_full_2d_lax_friedrichs_wave_family *wave_family)
{
  double x_L = C_L->x;
  double x_R = C_R->x;
  double y_L = C_L->y;
  double y_R = C_R->y;
  double f_L = U_L->f;
  double f_R = U_R->f;
  double ax = P->ax;
  double ay = P->ay;
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;

  wave_family->wave1.wave_f = (0.5 * ((f_R - f_L) - (((ay * f_R) - (ay * f_L)) / fabs(ay))));

  wave_family->wave2.wave_f = (0.5 * ((f_R - f_L) + (((ay * f_R) - (ay * f_L)) / fabs(ay))));
}

static inline advection_diffusion_full_2d_lax_friedrichs_waves_2D
create_waves()
{
  advection_diffusion_full_2d_lax_friedrichs_waves_2D waves;

  waves.x_wave_family = advection_diffusion_full_2d_lax_friedrichs_x_wave_family;
  waves.y_wave_family = advection_diffusion_full_2d_lax_friedrichs_y_wave_family;
  
  return waves;
}

static inline bool
advection_diffusion_full_2d_lax_friedrichs_x_waves_consistent(const advection_diffusion_full_2d_coordinates *C, const advection_diffusion_full_2d_state *U,
  const advection_diffusion_full_2d_parameters *P)
{
  advection_diffusion_full_2d_lax_friedrichs_wave_family wave_family;
  advection_diffusion_full_2d_lax_friedrichs_x_wave_family(C, C, U, U, P, &wave_family);

  double f_wave1 = wave_family.wave1.wave_f;
  double f_wave2 = wave_family.wave2.wave_f;

  return ((fabs(f_wave1) < 1.0e-8) && (fabs(f_wave2) < 1.0e-8));
}

static inline bool
advection_diffusion_full_2d_lax_friedrichs_y_waves_consistent(const advection_diffusion_full_2d_coordinates *C, const advection_diffusion_full_2d_state *U,
  const advection_diffusion_full_2d_parameters *P)
{
  advection_diffusion_full_2d_lax_friedrichs_wave_family wave_family;
  advection_diffusion_full_2d_lax_friedrichs_y_wave_family(C, C, U, U, P, &wave_family);

  double f_wave1 = wave_family.wave1.wave_f;
  double f_wave2 = wave_family.wave2.wave_f;

  return ((fabs(f_wave1) < 1.0e-8) && (fabs(f_wave2) < 1.0e-8));
}

static inline bool
advection_diffusion_full_2d_lax_friedrichs_x_waves_valid(const advection_diffusion_full_2d_coordinates *C_L, const advection_diffusion_full_2d_coordinates *C_R,
  const advection_diffusion_full_2d_state *U_L, const advection_diffusion_full_2d_state *U_R, const advection_diffusion_full_2d_parameters *P)
{
  double f_L = U_L->f;
  double f_R = U_R->f;

  double f_jump = f_R - f_L;

  advection_diffusion_full_2d_lax_friedrichs_wave_family wave_family;
  advection_diffusion_full_2d_lax_friedrichs_x_wave_family(C_L, C_R, U_L, U_R, P, &wave_family);

  double f_wave_sum = wave_family.wave1.wave_f + wave_family.wave2.wave_f;

  return (fabs(f_jump - f_wave_sum) < 1.0e-8);
}

static inline bool
advection_diffusion_full_2d_lax_friedrichs_y_waves_valid(const advection_diffusion_full_2d_coordinates *C_L, const advection_diffusion_full_2d_coordinates *C_R,
  const advection_diffusion_full_2d_state *U_L, const advection_diffusion_full_2d_state *U_R, const advection_diffusion_full_2d_parameters *P)
{
  double f_L = U_L->f;
  double f_R = U_R->f;

  double f_jump = f_R - f_L;

  advection_diffusion_full_2d_lax_friedrichs_wave_family wave_family;
  advection_diffusion_full_2d_lax_friedrichs_y_wave_family(C_L, C_R, U_L, U_R, P, &wave_family);

  double f_wave_sum = wave_family.wave1.wave_f + wave_family.wave2.wave_f;

  return (fabs(f_jump - f_wave_sum) < 1.0e-8);
}

typedef struct {
  double speed1;
  double speed2;
} advection_diffusion_full_2d_lax_friedrichs_speed_family;

typedef struct {
  void (*x_speed_family)(const advection_diffusion_full_2d_coordinates *C_L, const advection_diffusion_full_2d_coordinates *C_R,
    const advection_diffusion_full_2d_state *U_L, const advection_diffusion_full_2d_state *U_R,
    const advection_diffusion_full_2d_parameters *P, advection_diffusion_full_2d_lax_friedrichs_speed_family *speed_family);
  void (*y_speed_family)(const advection_diffusion_full_2d_coordinates *C_L, const advection_diffusion_full_2d_coordinates *C_R,
    const advection_diffusion_full_2d_state *U_L, const advection_diffusion_full_2d_state *U_R,
    const advection_diffusion_full_2d_parameters *P, advection_diffusion_full_2d_lax_friedrichs_speed_family *speed_family);
} advection_diffusion_full_2d_lax_friedrichs_speeds_2D;

static inline void
advection_diffusion_full_2d_lax_friedrichs_x_speed_family(const advection_diffusion_full_2d_coordinates *C_L, const advection_diffusion_full_2d_coordinates *C_R,
  const advection_diffusion_full_2d_state *U_L, const advection_diffusion_full_2d_state *U_R,
  const advection_diffusion_full_2d_parameters *P, advection_diffusion_full_2d_lax_friedrichs_speed_family *speed_family)
{
  double x_L = C_L->x;
  double x_R = C_R->x;
  double y_L = C_L->y;
  double y_R = C_R->y;
  double f_L = U_L->f;
  double f_R = U_R->f;
  double ax = P->ax;
  double ay = P->ay;
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;

  speed_family->speed1 = -(fabs(ax));
  speed_family->speed2 = fabs(ax);
}

static inline void
advection_diffusion_full_2d_lax_friedrichs_y_speed_family(const advection_diffusion_full_2d_coordinates *C_L, const advection_diffusion_full_2d_coordinates *C_R,
  const advection_diffusion_full_2d_state *U_L, const advection_diffusion_full_2d_state *U_R,
  const advection_diffusion_full_2d_parameters *P, advection_diffusion_full_2d_lax_friedrichs_speed_family *speed_family)
{
  double x_L = C_L->x;
  double x_R = C_R->x;
  double y_L = C_L->y;
  double y_R = C_R->y;
  double f_L = U_L->f;
  double f_R = U_R->f;
  double ax = P->ax;
  double ay = P->ay;
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;

  speed_family->speed1 = -(fabs(ay));
  speed_family->speed2 = fabs(ay);
}

static inline advection_diffusion_full_2d_lax_friedrichs_speeds_2D
create_speeds()
{
  advection_diffusion_full_2d_lax_friedrichs_speeds_2D speeds;

  speeds.x_speed_family = advection_diffusion_full_2d_lax_friedrichs_x_speed_family;
  speeds.y_speed_family = advection_diffusion_full_2d_lax_friedrichs_y_speed_family;

  return speeds;
}

typedef struct {
  double fluctuation_f;
} advection_diffusion_full_2d_lax_friedrichs_fluctuation;

typedef struct {
  void (*x_fluctuation)(const advection_diffusion_full_2d_coordinates *C_L, const advection_diffusion_full_2d_coordinates *C_R,
    const advection_diffusion_full_2d_state *U_L, const advection_diffusion_full_2d_state *U_R,
    const advection_diffusion_full_2d_parameters *P, advection_diffusion_full_2d_lax_friedrichs_fluctuation *fluctuation);
  void (*y_fluctuation)(const advection_diffusion_full_2d_coordinates *C_L, const advection_diffusion_full_2d_coordinates *C_R,
    const advection_diffusion_full_2d_state *U_L, const advection_diffusion_full_2d_state *U_R,
    const advection_diffusion_full_2d_parameters *P, advection_diffusion_full_2d_lax_friedrichs_fluctuation *fluctuation);
} advection_diffusion_full_2d_lax_friedrichs_fluctuations_2D;

static inline void
advection_diffusion_full_2d_lax_friedrichs_x_left_fluctuation(const advection_diffusion_full_2d_coordinates *C_L, const advection_diffusion_full_2d_coordinates *C_R,
  const advection_diffusion_full_2d_state *U_L, const advection_diffusion_full_2d_state *U_R,
  const advection_diffusion_full_2d_parameters *P, advection_diffusion_full_2d_lax_friedrichs_fluctuation *left_fluctuation)
{
  advection_diffusion_full_2d_lax_friedrichs_wave_family wave_family;
  advection_diffusion_full_2d_lax_friedrichs_x_wave_family(C_L, C_R, U_L, U_R, P, &wave_family);

  advection_diffusion_full_2d_lax_friedrichs_speed_family speed_family;
  advection_diffusion_full_2d_lax_friedrichs_x_speed_family(C_L, C_R, U_L, U_R, P, &speed_family);

  left_fluctuation->fluctuation_f = 
    (wave_family.wave1.wave_f * fmin(speed_family.speed1, 0.0)) + 
    (wave_family.wave2.wave_f * fmin(speed_family.speed2, 0.0));
}

static inline void
advection_diffusion_full_2d_lax_friedrichs_y_left_fluctuation(const advection_diffusion_full_2d_coordinates *C_L, const advection_diffusion_full_2d_coordinates *C_R,
  const advection_diffusion_full_2d_state *U_L, const advection_diffusion_full_2d_state *U_R,
  const advection_diffusion_full_2d_parameters *P, advection_diffusion_full_2d_lax_friedrichs_fluctuation *left_fluctuation)
{
  advection_diffusion_full_2d_lax_friedrichs_wave_family wave_family;
  advection_diffusion_full_2d_lax_friedrichs_y_wave_family(C_L, C_R, U_L, U_R, P, &wave_family);

  advection_diffusion_full_2d_lax_friedrichs_speed_family speed_family;
  advection_diffusion_full_2d_lax_friedrichs_y_speed_family(C_L, C_R, U_L, U_R, P, &speed_family);

  left_fluctuation->fluctuation_f = 
    (wave_family.wave1.wave_f * fmin(speed_family.speed1, 0.0)) + 
    (wave_family.wave2.wave_f * fmin(speed_family.speed2, 0.0));
}

static inline void
advection_diffusion_full_2d_lax_friedrichs_x_right_fluctuation(const advection_diffusion_full_2d_coordinates *C_L, const advection_diffusion_full_2d_coordinates *C_R,
  const advection_diffusion_full_2d_state *U_L, const advection_diffusion_full_2d_state *U_R,
  const advection_diffusion_full_2d_parameters *P, advection_diffusion_full_2d_lax_friedrichs_fluctuation *right_fluctuation)
{
  advection_diffusion_full_2d_lax_friedrichs_wave_family wave_family;
  advection_diffusion_full_2d_lax_friedrichs_x_wave_family(C_L, C_R, U_L, U_R, P, &wave_family);

  advection_diffusion_full_2d_lax_friedrichs_speed_family speed_family;
  advection_diffusion_full_2d_lax_friedrichs_x_speed_family(C_L, C_R, U_L, U_R, P, &speed_family);

  right_fluctuation->fluctuation_f = 
    (wave_family.wave1.wave_f * fmax(speed_family.speed1, 0.0)) + 
    (wave_family.wave2.wave_f * fmax(speed_family.speed2, 0.0));
}

static inline void
advection_diffusion_full_2d_lax_friedrichs_y_right_fluctuation(const advection_diffusion_full_2d_coordinates *C_L, const advection_diffusion_full_2d_coordinates *C_R,
  const advection_diffusion_full_2d_state *U_L, const advection_diffusion_full_2d_state *U_R,
  const advection_diffusion_full_2d_parameters *P, advection_diffusion_full_2d_lax_friedrichs_fluctuation *right_fluctuation)
{
  advection_diffusion_full_2d_lax_friedrichs_wave_family wave_family;
  advection_diffusion_full_2d_lax_friedrichs_y_wave_family(C_L, C_R, U_L, U_R, P, &wave_family);

  advection_diffusion_full_2d_lax_friedrichs_speed_family speed_family;
  advection_diffusion_full_2d_lax_friedrichs_y_speed_family(C_L, C_R, U_L, U_R, P, &speed_family);

  right_fluctuation->fluctuation_f = 
    (wave_family.wave1.wave_f * fmax(speed_family.speed1, 0.0)) + 
    (wave_family.wave2.wave_f * fmax(speed_family.speed2, 0.0));
}

static inline advection_diffusion_full_2d_lax_friedrichs_fluctuations_2D
create_left_fluctuations()
{
  advection_diffusion_full_2d_lax_friedrichs_fluctuations_2D left_fluctuations;

  left_fluctuations.x_fluctuation = advection_diffusion_full_2d_lax_friedrichs_x_left_fluctuation;
  left_fluctuations.y_fluctuation = advection_diffusion_full_2d_lax_friedrichs_y_left_fluctuation;

  return left_fluctuations;
}

static inline advection_diffusion_full_2d_lax_friedrichs_fluctuations_2D
create_right_fluctuations()
{
  advection_diffusion_full_2d_lax_friedrichs_fluctuations_2D right_fluctuations;

  right_fluctuations.x_fluctuation = advection_diffusion_full_2d_lax_friedrichs_x_right_fluctuation;
  right_fluctuations.y_fluctuation = advection_diffusion_full_2d_lax_friedrichs_y_right_fluctuation;

  return right_fluctuations;
}

static inline bool
advection_diffusion_full_2d_lax_friedrichs_x_left_fluctuations_consistent(const advection_diffusion_full_2d_coordinates *C, const advection_diffusion_full_2d_state *U,
  const advection_diffusion_full_2d_parameters *P)
{
  advection_diffusion_full_2d_lax_friedrichs_fluctuation left_fluctuation;
  advection_diffusion_full_2d_lax_friedrichs_x_left_fluctuation(C, C, U, U, P, &left_fluctuation);

  double f_left_fluctuation = left_fluctuation.fluctuation_f;

  return ((fabs(f_left_fluctuation) < 1.0e-8));
}

static inline bool
advection_diffusion_full_2d_lax_friedrichs_y_left_fluctuations_consistent(const advection_diffusion_full_2d_coordinates *C, const advection_diffusion_full_2d_state *U,
  const advection_diffusion_full_2d_parameters *P)
{
  advection_diffusion_full_2d_lax_friedrichs_fluctuation left_fluctuation;
  advection_diffusion_full_2d_lax_friedrichs_y_left_fluctuation(C, C, U, U, P, &left_fluctuation);

  double f_left_fluctuation = left_fluctuation.fluctuation_f;

  return ((fabs(f_left_fluctuation) < 1.0e-8));
}

static inline bool
advection_diffusion_full_2d_lax_friedrichs_x_right_fluctuations_consistent(const advection_diffusion_full_2d_coordinates *C, const advection_diffusion_full_2d_state *U,
  const advection_diffusion_full_2d_parameters *P)
{
  advection_diffusion_full_2d_lax_friedrichs_fluctuation right_fluctuation;
  advection_diffusion_full_2d_lax_friedrichs_x_right_fluctuation(C, C, U, U, P, &right_fluctuation);

  double f_right_fluctuation = right_fluctuation.fluctuation_f;

  return ((fabs(f_right_fluctuation) < 1.0e-8));
}

static inline bool
advection_diffusion_full_2d_lax_friedrichs_y_right_fluctuations_consistent(const advection_diffusion_full_2d_coordinates *C, const advection_diffusion_full_2d_state *U,
  const advection_diffusion_full_2d_parameters *P)
{
  advection_diffusion_full_2d_lax_friedrichs_fluctuation right_fluctuation;
  advection_diffusion_full_2d_lax_friedrichs_y_right_fluctuation(C, C, U, U, P, &right_fluctuation);

  double f_right_fluctuation = right_fluctuation.fluctuation_f;

  return ((fabs(f_right_fluctuation) < 1.0e-8));
}

static inline bool
advection_diffusion_full_2d_lax_friedrichs_x_fluctuations_valid(const advection_diffusion_full_2d_coordinates *C_L, const advection_diffusion_full_2d_coordinates *C_R,
  const advection_diffusion_full_2d_state *U_L, const advection_diffusion_full_2d_state *U_R, const advection_diffusion_full_2d_parameters *P)
{
  advection_diffusion_full_2d_flux left_flux;
  advection_diffusion_full_2d_x_flux(C_L, P, U_L, &left_flux);

  advection_diffusion_full_2d_flux right_flux;
  advection_diffusion_full_2d_x_flux(C_R, P, U_R, &right_flux);

  double f_flux_jump = right_flux.flux_f - left_flux.flux_f;

  advection_diffusion_full_2d_lax_friedrichs_fluctuation left_fluctuation;
  advection_diffusion_full_2d_lax_friedrichs_x_left_fluctuation(C_L, C_R, U_L, U_R, P, &left_fluctuation);

  advection_diffusion_full_2d_lax_friedrichs_fluctuation right_fluctuation;
  advection_diffusion_full_2d_lax_friedrichs_x_right_fluctuation(C_L, C_R, U_L, U_R, P, &right_fluctuation);

  double f_fluctuation_sum = left_fluctuation.fluctuation_f + right_fluctuation.fluctuation_f;

  return (fabs(f_flux_jump - f_fluctuation_sum) < 1.0e-8);
}

static inline bool
advection_diffusion_full_2d_lax_friedrichs_y_fluctuations_valid(const advection_diffusion_full_2d_coordinates *C_L, const advection_diffusion_full_2d_coordinates *C_R,
  const advection_diffusion_full_2d_state *U_L, const advection_diffusion_full_2d_state *U_R, const advection_diffusion_full_2d_parameters *P)
{
  advection_diffusion_full_2d_flux left_flux;
  advection_diffusion_full_2d_y_flux(C_L, P, U_L, &left_flux);

  advection_diffusion_full_2d_flux right_flux;
  advection_diffusion_full_2d_y_flux(C_R, P, U_R, &right_flux);

  double f_flux_jump = right_flux.flux_f - left_flux.flux_f;

  advection_diffusion_full_2d_lax_friedrichs_fluctuation left_fluctuation;
  advection_diffusion_full_2d_lax_friedrichs_y_left_fluctuation(C_L, C_R, U_L, U_R, P, &left_fluctuation);

  advection_diffusion_full_2d_lax_friedrichs_fluctuation right_fluctuation;
  advection_diffusion_full_2d_lax_friedrichs_y_right_fluctuation(C_L, C_R, U_L, U_R, P, &right_fluctuation);

  double f_fluctuation_sum = left_fluctuation.fluctuation_f + right_fluctuation.fluctuation_f;

  return (fabs(f_flux_jump - f_fluctuation_sum) < 1.0e-8);
}

typedef struct {
  double reconstruction_f;
} advection_diffusion_full_2d_minmod_reconstruction;

typedef struct {
  void (*x_reconstruction)(const advection_diffusion_full_2d_coordinates *C_L, const advection_diffusion_full_2d_coordinates *C, const advection_diffusion_full_2d_coordinates *C_R,
    const advection_diffusion_full_2d_state *U_L, const advection_diffusion_full_2d_state *U, const advection_diffusion_full_2d_state *U_R,
    const advection_diffusion_full_2d_parameters *P, advection_diffusion_full_2d_minmod_reconstruction *reconstruction);
  void (*y_reconstruction)(const advection_diffusion_full_2d_coordinates *C_L, const advection_diffusion_full_2d_coordinates *C, const advection_diffusion_full_2d_coordinates *C_R,
    const advection_diffusion_full_2d_state *U_L, const advection_diffusion_full_2d_state *U, const advection_diffusion_full_2d_state *U_R,
    const advection_diffusion_full_2d_parameters *P, advection_diffusion_full_2d_minmod_reconstruction *reconstruction);
} advection_diffusion_full_2d_minmod_reconstruction_2D;

static inline void
advection_diffusion_full_2d_minmod_x_left_reconstruction(const advection_diffusion_full_2d_coordinates *C_L, const advection_diffusion_full_2d_coordinates *C,
  const advection_diffusion_full_2d_coordinates *C_R, const advection_diffusion_full_2d_state *U_L, const advection_diffusion_full_2d_state *U, const advection_diffusion_full_2d_state *U_R,
  const advection_diffusion_full_2d_parameters *P, advection_diffusion_full_2d_minmod_reconstruction *left_reconstruction)
{
  double x_L = C_L->x;
  double x = C->x; 
  double x_R = C_R->x;
  double y_L = C_L->y;
  double y = C->y; 
  double y_R = C_R->y;
  double f_L = U_L->f;
  double f = U->f;
  double f_R = U_R->f;
  double ax = P->ax;
  double ay = P->ay;
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;

  left_reconstruction->reconstruction_f = (f - (0.5 * (fmax(0.0, fmin((f - f_L), (f_R - f))) + fmin(0.0, fmax((f - f_L),
    (f_R - f))))));
}

static inline void
advection_diffusion_full_2d_minmod_y_left_reconstruction(const advection_diffusion_full_2d_coordinates *C_L, const advection_diffusion_full_2d_coordinates *C,
  const advection_diffusion_full_2d_coordinates *C_R, const advection_diffusion_full_2d_state *U_L, const advection_diffusion_full_2d_state *U, const advection_diffusion_full_2d_state *U_R,
  const advection_diffusion_full_2d_parameters *P, advection_diffusion_full_2d_minmod_reconstruction *left_reconstruction)
{
  double x_L = C_L->x;
  double x = C->x; 
  double x_R = C_R->x;
  double y_L = C_L->y;
  double y = C->y; 
  double y_R = C_R->y;
  double f_L = U_L->f;
  double f = U->f;
  double f_R = U_R->f;
  double ax = P->ax;
  double ay = P->ay;
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;

  left_reconstruction->reconstruction_f = (f - (0.5 * (fmax(0.0, fmin((f - f_L), (f_R - f))) + fmin(0.0, fmax((f - f_L),
    (f_R - f))))));
}

static inline void
advection_diffusion_full_2d_minmod_x_right_reconstruction(const advection_diffusion_full_2d_coordinates *C_L, const advection_diffusion_full_2d_coordinates *C,
  const advection_diffusion_full_2d_coordinates *C_R, const advection_diffusion_full_2d_state *U_L, const advection_diffusion_full_2d_state *U, const advection_diffusion_full_2d_state *U_R,
  const advection_diffusion_full_2d_parameters *P, advection_diffusion_full_2d_minmod_reconstruction *right_reconstruction)
{
  double x_L = C_L->x;
  double x = C->x; 
  double x_R = C_R->x;
  double y_L = C_L->y;
  double y = C->y; 
  double y_R = C_R->y;
  double f_L = U_L->f;
  double f = U->f;
  double f_R = U_R->f;
  double ax = P->ax;
  double ay = P->ay;
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;

  right_reconstruction->reconstruction_f = (f + (0.5 * (fmax(0.0, fmin((f - f_L), (f_R - f))) + fmin(0.0, fmax((f -
    f_L), (f_R - f))))));
}

static inline void
advection_diffusion_full_2d_minmod_y_right_reconstruction(const advection_diffusion_full_2d_coordinates *C_L, const advection_diffusion_full_2d_coordinates *C,
  const advection_diffusion_full_2d_coordinates *C_R, const advection_diffusion_full_2d_state *U_L, const advection_diffusion_full_2d_state *U, const advection_diffusion_full_2d_state *U_R,
  const advection_diffusion_full_2d_parameters *P, advection_diffusion_full_2d_minmod_reconstruction *right_reconstruction)
{
  double x_L = C_L->x;
  double x = C->x; 
  double x_R = C_R->x;
  double y_L = C_L->y;
  double y = C->y; 
  double y_R = C_R->y;
  double f_L = U_L->f;
  double f = U->f;
  double f_R = U_R->f;
  double ax = P->ax;
  double ay = P->ay;
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;

  right_reconstruction->reconstruction_f = (f + (0.5 * (fmax(0.0, fmin((f - f_L), (f_R - f))) + fmin(0.0, fmax((f -
    f_L), (f_R - f))))));
}

static inline advection_diffusion_full_2d_minmod_reconstruction_2D
create_left_reconstruction()
{
  advection_diffusion_full_2d_minmod_reconstruction_2D left_reconstruction;

  left_reconstruction.x_reconstruction = advection_diffusion_full_2d_minmod_x_left_reconstruction;
  left_reconstruction.y_reconstruction = advection_diffusion_full_2d_minmod_y_left_reconstruction;

  return left_reconstruction;
}

static inline advection_diffusion_full_2d_minmod_reconstruction_2D
create_right_reconstruction()
{
  advection_diffusion_full_2d_minmod_reconstruction_2D right_reconstruction;

  right_reconstruction.x_reconstruction = advection_diffusion_full_2d_minmod_x_right_reconstruction;
  right_reconstruction.y_reconstruction = advection_diffusion_full_2d_minmod_y_right_reconstruction;

  return right_reconstruction;
}

static inline bool
advection_diffusion_full_2d_minmod_x_left_reconstruction_consistent(const advection_diffusion_full_2d_coordinates *C, const advection_diffusion_full_2d_state *U,
const advection_diffusion_full_2d_parameters *P)
{
  advection_diffusion_full_2d_minmod_reconstruction left_reconstruction;
  advection_diffusion_full_2d_minmod_x_left_reconstruction(C, C, C, U, U, U, P, &left_reconstruction);

  double f_left_reconstruction = left_reconstruction.reconstruction_f;
  double f = U->f;

  return ((fabs(f_left_reconstruction - f) < 1.0e-8));
}

static inline bool
advection_diffusion_full_2d_minmod_y_left_reconstruction_consistent(const advection_diffusion_full_2d_coordinates *C, const advection_diffusion_full_2d_state *U,
const advection_diffusion_full_2d_parameters *P)
{
  advection_diffusion_full_2d_minmod_reconstruction left_reconstruction;
  advection_diffusion_full_2d_minmod_y_left_reconstruction(C, C, C, U, U, U, P, &left_reconstruction);

  double f_left_reconstruction = left_reconstruction.reconstruction_f;
  double f = U->f;

  return ((fabs(f_left_reconstruction - f) < 1.0e-8));
}

static inline bool
advection_diffusion_full_2d_minmod_x_right_reconstruction_consistent(const advection_diffusion_full_2d_coordinates *C, const advection_diffusion_full_2d_state *U,
const advection_diffusion_full_2d_parameters *P)
{
  advection_diffusion_full_2d_minmod_reconstruction right_reconstruction;
  advection_diffusion_full_2d_minmod_x_right_reconstruction(C, C, C, U, U, U, P, &right_reconstruction);

  double f_right_reconstruction = right_reconstruction.reconstruction_f;
  double f = U->f;

  return ((fabs(f_right_reconstruction - f) < 1.0e-8));
}

static inline bool
advection_diffusion_full_2d_minmod_y_right_reconstruction_consistent(const advection_diffusion_full_2d_coordinates *C, const advection_diffusion_full_2d_state *U,
const advection_diffusion_full_2d_parameters *P)
{
  advection_diffusion_full_2d_minmod_reconstruction right_reconstruction;
  advection_diffusion_full_2d_minmod_y_right_reconstruction(C, C, C, U, U, U, P, &right_reconstruction);

  double f_right_reconstruction = right_reconstruction.reconstruction_f;
  double f = U->f;

  return ((fabs(f_right_reconstruction - f) < 1.0e-8));
}

typedef struct {
  double interface_f;
} advection_diffusion_full_2d_centered_diffusion_interface_state;

typedef struct {
  void (*x_interface_state)(const advection_diffusion_full_2d_coordinates *C_Lx, const advection_diffusion_full_2d_coordinates *C_Rx,
    const advection_diffusion_full_2d_state *U_Lx, const advection_diffusion_full_2d_state *U_Rx,
    const advection_diffusion_full_2d_parameters *P, advection_diffusion_full_2d_centered_diffusion_interface_state *interface_state);
  void (*y_interface_state)(const advection_diffusion_full_2d_coordinates *C_Ly, const advection_diffusion_full_2d_coordinates *C_Ry,
    const advection_diffusion_full_2d_state *U_Ly, const advection_diffusion_full_2d_state *U_Ry,
    const advection_diffusion_full_2d_parameters *P, advection_diffusion_full_2d_centered_diffusion_interface_state *interface_state);
} advection_diffusion_full_2d_centered_diffusion_interface_states_2D;

static inline void
advection_diffusion_full_2d_centered_diffusion_x_interface_state(const advection_diffusion_full_2d_coordinates *C_Lx, const advection_diffusion_full_2d_coordinates *C_Rx,
  const advection_diffusion_full_2d_state *U_Lx, const advection_diffusion_full_2d_state *U_Rx,
  const advection_diffusion_full_2d_parameters *P, advection_diffusion_full_2d_centered_diffusion_interface_state *interface_state)
{
  double x_Lx = C_Lx->x;
  double x_Rx = C_Rx->x;
  double y_Lx = C_Lx->y;
  double y_Rx = C_Rx->y;
  double f_Lx = U_Lx->f;
  double f_Rx = U_Rx->f;
  double ax = P->ax;
  double ay = P->ay;
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;

  interface_state->interface_f = (0.5 * (f_Lx + f_Rx));
}

static inline void
advection_diffusion_full_2d_centered_diffusion_y_interface_state(const advection_diffusion_full_2d_coordinates *C_Ly, const advection_diffusion_full_2d_coordinates *C_Ry,
  const advection_diffusion_full_2d_state *U_Ly, const advection_diffusion_full_2d_state *U_Ry,
  const advection_diffusion_full_2d_parameters *P, advection_diffusion_full_2d_centered_diffusion_interface_state *interface_state)
{
  double x_Ly = C_Ly->x;
  double x_Ry = C_Ry->x;
  double y_Ly = C_Ly->y;
  double y_Ry = C_Ry->y;
  double f_Ly = U_Ly->f;
  double f_Ry = U_Ry->f;
  double ax = P->ax;
  double ay = P->ay;
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;

  interface_state->interface_f = (0.5 * (f_Ly + f_Ry));
}

static inline advection_diffusion_full_2d_centered_diffusion_interface_states_2D
create_interface_states()
{
  advection_diffusion_full_2d_centered_diffusion_interface_states_2D interface_states;

  interface_states.x_interface_state = advection_diffusion_full_2d_centered_diffusion_x_interface_state;
  interface_states.y_interface_state = advection_diffusion_full_2d_centered_diffusion_y_interface_state;
  
  return interface_states;
}

typedef struct {
  double interface_f_x;
  double interface_f_y;
} advection_diffusion_full_2d_centered_diffusion_interface_gradient;

typedef struct {
  void (*x_interface_gradient)(const advection_diffusion_full_2d_coordinates *C_Lx, const advection_diffusion_full_2d_coordinates *C_Rx,
    const advection_diffusion_full_2d_state *U_Lx, const advection_diffusion_full_2d_state *U_Rx,
    const advection_diffusion_full_2d_coordinates *C_Lx_Ly, const advection_diffusion_full_2d_coordinates *C_Lx_Ry,
    const advection_diffusion_full_2d_coordinates *C_Rx_Ly, const advection_diffusion_full_2d_coordinates *C_Rx_Ry,
    const advection_diffusion_full_2d_state *U_Lx_Ly, const advection_diffusion_full_2d_state *U_Lx_Ry,
    const advection_diffusion_full_2d_state *U_Rx_Ly, const advection_diffusion_full_2d_state *U_Rx_Ry,
    const advection_diffusion_full_2d_parameters *P, advection_diffusion_full_2d_centered_diffusion_interface_gradient *interface_gradient);
  void (*y_interface_gradient)(const advection_diffusion_full_2d_coordinates *C_Ly_Lx, const advection_diffusion_full_2d_coordinates *C_Ly_Rx,
    const advection_diffusion_full_2d_coordinates *C_Ry_Lx, const advection_diffusion_full_2d_coordinates *C_Ry_Rx,
    const advection_diffusion_full_2d_state *U_Ly_Lx, const advection_diffusion_full_2d_state *U_Ly_Rx,
    const advection_diffusion_full_2d_state *U_Ry_Lx, const advection_diffusion_full_2d_state *U_Ry_Rx,
    const advection_diffusion_full_2d_coordinates *C_Ly, const advection_diffusion_full_2d_coordinates *C_Ry,
    const advection_diffusion_full_2d_state *U_Ly, const advection_diffusion_full_2d_state *U_Ry,
    const advection_diffusion_full_2d_parameters *P, advection_diffusion_full_2d_centered_diffusion_interface_gradient *interface_gradient);
} advection_diffusion_full_2d_centered_diffusion_interface_gradients_2D;

static inline void
advection_diffusion_full_2d_centered_diffusion_x_interface_gradient(const advection_diffusion_full_2d_coordinates *C_Lx, const advection_diffusion_full_2d_coordinates *C_Rx,
  const advection_diffusion_full_2d_state *U_Lx, const advection_diffusion_full_2d_state *U_Rx,
  const advection_diffusion_full_2d_coordinates *C_Lx_Ly, const advection_diffusion_full_2d_coordinates *C_Lx_Ry,
  const advection_diffusion_full_2d_coordinates *C_Rx_Ly, const advection_diffusion_full_2d_coordinates *C_Rx_Ry,
  const advection_diffusion_full_2d_state *U_Lx_Ly, const advection_diffusion_full_2d_state *U_Lx_Ry,
  const advection_diffusion_full_2d_state *U_Rx_Ly, const advection_diffusion_full_2d_state *U_Rx_Ry,
  const advection_diffusion_full_2d_parameters *P, advection_diffusion_full_2d_centered_diffusion_interface_gradient *interface_gradient)
{
  double x_Lx = C_Lx->x;
  double x_Rx = C_Rx->x;
  double y_Lx_Ly = C_Lx_Ly->y;
  double y_Lx_Ry = C_Lx_Ry->y;
  double y_Rx_Ly = C_Rx_Ly->y;
  double y_Rx_Ry = C_Rx_Ry->y;

  double f_Lx = U_Lx->f;
  double f_Rx = U_Rx->f;
  double f_Lx_Ly = U_Lx_Ly->f;
  double f_Lx_Ry = U_Lx_Ry->f;
  double f_Rx_Ly = U_Rx_Ly->f;
  double f_Rx_Ry = U_Rx_Ry->f;

  double ax = P->ax;
  double ay = P->ay;
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;

  interface_gradient->interface_f_x = ((f_Rx - f_Lx) / (x_Rx - x_Lx));

  interface_gradient->interface_f_y = (0.5 * (((f_Lx_Ry - f_Lx_Ly) / (y_Lx_Ry - y_Lx_Ly)) + ((f_Rx_Ry - f_Rx_Ly) /
    (y_Rx_Ry - y_Rx_Ly))));
}

static inline void
advection_diffusion_full_2d_centered_diffusion_y_interface_gradient(const advection_diffusion_full_2d_coordinates *C_Ly_Lx, const advection_diffusion_full_2d_coordinates *C_Ly_Rx,
  const advection_diffusion_full_2d_coordinates *C_Ry_Lx, const advection_diffusion_full_2d_coordinates *C_Ry_Rx,
  const advection_diffusion_full_2d_state *U_Ly_Lx, const advection_diffusion_full_2d_state *U_Ly_Rx,
  const advection_diffusion_full_2d_state *U_Ry_Lx, const advection_diffusion_full_2d_state *U_Ry_Rx,
  const advection_diffusion_full_2d_coordinates *C_Ly, const advection_diffusion_full_2d_coordinates *C_Ry,
  const advection_diffusion_full_2d_state *U_Ly, const advection_diffusion_full_2d_state *U_Ry,
  const advection_diffusion_full_2d_parameters *P, advection_diffusion_full_2d_centered_diffusion_interface_gradient *interface_gradient)
{
  double x_Ly_Lx = C_Ly_Lx->x;
  double x_Ly_Rx = C_Ly_Rx->x;
  double x_Ry_Lx = C_Ry_Lx->x;
  double x_Ry_Rx = C_Ry_Rx->x;
  double y_Ly = C_Ly->y;
  double y_Ry = C_Ry->y;

  double f_Ly_Lx = U_Ly_Lx->f;
  double f_Ly_Rx = U_Ly_Rx->f;
  double f_Ry_Lx = U_Ry_Lx->f;
  double f_Ry_Rx = U_Ry_Rx->f;
  double f_Ly = U_Ly->f;
  double f_Ry = U_Ry->f;

  double ax = P->ax;
  double ay = P->ay;
  double Dxx = P->Dxx;
  double Dxy = P->Dxy;
  double Dyx = P->Dyx;
  double Dyy = P->Dyy;

  interface_gradient->interface_f_x = (0.5 * (((f_Ly_Rx - f_Ly_Lx) / (x_Ly_Rx - x_Ly_Lx)) + ((f_Ry_Rx - f_Ry_Lx) /
    (x_Ry_Rx - x_Ry_Lx))));

  interface_gradient->interface_f_y = ((f_Ry - f_Ly) / (y_Ry - y_Ly));
}

static inline advection_diffusion_full_2d_centered_diffusion_interface_gradients_2D
create_interface_gradients()
{
  advection_diffusion_full_2d_centered_diffusion_interface_gradients_2D interface_gradients;

  interface_gradients.x_interface_gradient = advection_diffusion_full_2d_centered_diffusion_x_interface_gradient;
  interface_gradients.y_interface_gradient = advection_diffusion_full_2d_centered_diffusion_y_interface_gradient;
  
  return interface_gradients;
}

int
main(int argc, char *argv[])
{
  // Insert simulation drivers here.
  
  return 0;
}
                