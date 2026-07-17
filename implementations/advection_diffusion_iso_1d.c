#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

typedef struct {
  double x;
} advection_diffusion_iso_1d_coordinates;

typedef struct {
  double f;
} advection_diffusion_iso_1d_state;

static inline bool
advection_diffusion_iso_1d_state_valid(const advection_diffusion_iso_1d_state *U)
{
  double f = U->f;
  
  return true;
}

typedef struct {
  double a;
  double D;
} advection_diffusion_iso_1d_parameters;

static inline bool
advection_diffusion_iso_1d_parameters_valid(const advection_diffusion_iso_1d_parameters *P)
{
  double a = P->a;
  double D = P->D;

  return (D > 0.0);
}

typedef struct {
  double flux_f;
} advection_diffusion_iso_1d_flux;

typedef struct {
  void (*x_flux)(const advection_diffusion_iso_1d_coordinates *C, const advection_diffusion_iso_1d_parameters *P,
    const advection_diffusion_iso_1d_state *U, advection_diffusion_iso_1d_flux *flux);
} advection_diffusion_iso_1d_fluxes_1D;

static inline void
advection_diffusion_iso_1d_x_flux(const advection_diffusion_iso_1d_coordinates *C, const advection_diffusion_iso_1d_parameters *P,
  const advection_diffusion_iso_1d_state *U, advection_diffusion_iso_1d_flux *flux)
{
  double x = C->x;
  double a = P->a;
  double D = P->D;
  double f = U->f;

  flux->flux_f = (a * f);
}

static inline advection_diffusion_iso_1d_fluxes_1D
create_fluxes()
{
  advection_diffusion_iso_1d_fluxes_1D fluxes;

  fluxes.x_flux = advection_diffusion_iso_1d_x_flux;
  
  return fluxes;
}

typedef struct {
  double mu1;
} advection_diffusion_iso_1d_wavespeed;

typedef struct {
  void (*x_wavespeed)(const advection_diffusion_iso_1d_coordinates *C, const advection_diffusion_iso_1d_parameters *P,
    const advection_diffusion_iso_1d_state *U, advection_diffusion_iso_1d_wavespeed *wavespeed);
} advection_diffusion_iso_1d_wavespeed_estimate_1D;

static inline void
advection_diffusion_iso_1d_x_wavespeed(const advection_diffusion_iso_1d_coordinates *C, const advection_diffusion_iso_1d_parameters *P,
  const advection_diffusion_iso_1d_state *U, advection_diffusion_iso_1d_wavespeed *wavespeed)
{
  double x = C->x;
  double a = P->a;
  double D = P->D;
  double f = U->f;

  wavespeed->mu1 = a;
}

static inline advection_diffusion_iso_1d_wavespeed_estimate_1D
create_wavespeeds()
{
  advection_diffusion_iso_1d_wavespeed_estimate_1D wavespeed_estimate;

  wavespeed_estimate.x_wavespeed = advection_diffusion_iso_1d_x_wavespeed;
  
  return wavespeed_estimate;
}

typedef struct {
  double f_x;
} advection_diffusion_iso_1d_gradient;

typedef struct {
  double diffusive_flux_f;
} advection_diffusion_iso_1d_diffusive_flux;

typedef struct {
  void (*x_diffusive_flux)(const advection_diffusion_iso_1d_coordinates *C, const advection_diffusion_iso_1d_parameters *P,
    const advection_diffusion_iso_1d_state *U, advection_diffusion_iso_1d_gradient *DU, advection_diffusion_iso_1d_diffusive_flux *flux);
} advection_diffusion_iso_1d_diffusive_fluxes_1D;

static inline void
advection_diffusion_iso_1d_x_diffusive_flux(const advection_diffusion_iso_1d_coordinates *C, const advection_diffusion_iso_1d_parameters *P,
  const advection_diffusion_iso_1d_state *U, advection_diffusion_iso_1d_gradient *DU, advection_diffusion_iso_1d_diffusive_flux *diffusive_flux)
{
  double x = C->x;
  double a = P->a;
  double D = P->D;
  double f = U->f;
  double f_x = DU->f_x;

  diffusive_flux->diffusive_flux_f = (D * f_x);
}

static inline advection_diffusion_iso_1d_diffusive_fluxes_1D
create_diffusive_fluxes()
{
  advection_diffusion_iso_1d_diffusive_fluxes_1D diffusive_fluxes;

  diffusive_fluxes.x_diffusive_flux = advection_diffusion_iso_1d_x_diffusive_flux;
  
  return diffusive_fluxes;
}

typedef struct {
  double wave_f;
} advection_diffusion_iso_1d_lax_friedrichs_wave;

typedef struct {
  advection_diffusion_iso_1d_lax_friedrichs_wave wave1;
  advection_diffusion_iso_1d_lax_friedrichs_wave wave2;
} advection_diffusion_iso_1d_lax_friedrichs_wave_family;

typedef struct {
  void (*x_wave_family)(const advection_diffusion_iso_1d_coordinates *C_L, const advection_diffusion_iso_1d_coordinates *C_R,
    const advection_diffusion_iso_1d_state *U_L, const advection_diffusion_iso_1d_state *U_R,
    const advection_diffusion_iso_1d_parameters *P, advection_diffusion_iso_1d_lax_friedrichs_wave_family *wave_family);
} advection_diffusion_iso_1d_lax_friedrichs_waves_1D;

static inline void
advection_diffusion_iso_1d_lax_friedrichs_x_wave_family(const advection_diffusion_iso_1d_coordinates *C_L, const advection_diffusion_iso_1d_coordinates *C_R,
  const advection_diffusion_iso_1d_state *U_L, const advection_diffusion_iso_1d_state *U_R,
  const advection_diffusion_iso_1d_parameters *P, advection_diffusion_iso_1d_lax_friedrichs_wave_family *wave_family)
{
  double x_L = C_L->x;
  double x_R = C_R->x;
  double f_L = U_L->f;
  double f_R = U_R->f;
  double a = P->a;
  double D = P->D;

  wave_family->wave1.wave_f = (0.5 * ((f_R - f_L) - (((a * f_R) - (a * f_L)) / fabs(a))));

  wave_family->wave2.wave_f = (0.5 * ((f_R - f_L) + (((a * f_R) - (a * f_L)) / fabs(a))));
}

static inline advection_diffusion_iso_1d_lax_friedrichs_waves_1D
create_waves()
{
  advection_diffusion_iso_1d_lax_friedrichs_waves_1D waves;

  waves.x_wave_family = advection_diffusion_iso_1d_lax_friedrichs_x_wave_family;
  
  return waves;
}

static inline bool
advection_diffusion_iso_1d_lax_friedrichs_x_waves_consistent(const advection_diffusion_iso_1d_coordinates *C, const advection_diffusion_iso_1d_state *U,
  const advection_diffusion_iso_1d_parameters *P)
{
  advection_diffusion_iso_1d_lax_friedrichs_wave_family wave_family;
  advection_diffusion_iso_1d_lax_friedrichs_x_wave_family(C, C, U, U, P, &wave_family);

  double f_wave1 = wave_family.wave1.wave_f;
  double f_wave2 = wave_family.wave2.wave_f;

  return ((fabs(f_wave1) < 1.0e-8) && (fabs(f_wave2) < 1.0e-8));
}

static inline bool
advection_diffusion_iso_1d_lax_friedrichs_x_waves_valid(const advection_diffusion_iso_1d_coordinates *C_L, const advection_diffusion_iso_1d_coordinates *C_R,
  const advection_diffusion_iso_1d_state *U_L, const advection_diffusion_iso_1d_state *U_R, const advection_diffusion_iso_1d_parameters *P)
{
  double f_L = U_L->f;
  double f_R = U_R->f;

  double f_jump = f_R - f_L;

  advection_diffusion_iso_1d_lax_friedrichs_wave_family wave_family;
  advection_diffusion_iso_1d_lax_friedrichs_x_wave_family(C_L, C_R, U_L, U_R, P, &wave_family);

  double f_wave_sum = wave_family.wave1.wave_f + wave_family.wave2.wave_f;

  return (fabs(f_jump - f_wave_sum) < 1.0e-8);
}

typedef struct {
  double speed1;
  double speed2;
} advection_diffusion_iso_1d_lax_friedrichs_speed_family;

typedef struct {
  void (*x_speed_family)(const advection_diffusion_iso_1d_coordinates *C_L, const advection_diffusion_iso_1d_coordinates *C_R,
    const advection_diffusion_iso_1d_state *U_L, const advection_diffusion_iso_1d_state *U_R,
    const advection_diffusion_iso_1d_parameters *P, advection_diffusion_iso_1d_lax_friedrichs_speed_family *speed_family);
} advection_diffusion_iso_1d_lax_friedrichs_speeds_1D;

static inline void
advection_diffusion_iso_1d_lax_friedrichs_x_speed_family(const advection_diffusion_iso_1d_coordinates *C_L, const advection_diffusion_iso_1d_coordinates *C_R,
  const advection_diffusion_iso_1d_state *U_L, const advection_diffusion_iso_1d_state *U_R,
  const advection_diffusion_iso_1d_parameters *P, advection_diffusion_iso_1d_lax_friedrichs_speed_family *speed_family)
{
  double x_L = C_L->x;
  double x_R = C_R->x;
  double f_L = U_L->f;
  double f_R = U_R->f;
  double a = P->a;
  double D = P->D;

  speed_family->speed1 = -(fabs(a));
  speed_family->speed2 = fabs(a);
}

static inline advection_diffusion_iso_1d_lax_friedrichs_speeds_1D
create_speeds()
{
  advection_diffusion_iso_1d_lax_friedrichs_speeds_1D speeds;

  speeds.x_speed_family = advection_diffusion_iso_1d_lax_friedrichs_x_speed_family;

  return speeds;
}

typedef struct {
  double fluctuation_f;
} advection_diffusion_iso_1d_lax_friedrichs_fluctuation;

typedef struct {
  void (*x_fluctuation)(const advection_diffusion_iso_1d_coordinates *C_L, const advection_diffusion_iso_1d_coordinates *C_R,
    const advection_diffusion_iso_1d_state *U_L, const advection_diffusion_iso_1d_state *U_R,
    const advection_diffusion_iso_1d_parameters *P, advection_diffusion_iso_1d_lax_friedrichs_fluctuation *fluctuation);
} advection_diffusion_iso_1d_lax_friedrichs_fluctuations_1D;

static inline void
advection_diffusion_iso_1d_lax_friedrichs_x_left_fluctuation(const advection_diffusion_iso_1d_coordinates *C_L, const advection_diffusion_iso_1d_coordinates *C_R,
  const advection_diffusion_iso_1d_state *U_L, const advection_diffusion_iso_1d_state *U_R,
  const advection_diffusion_iso_1d_parameters *P, advection_diffusion_iso_1d_lax_friedrichs_fluctuation *left_fluctuation)
{
  advection_diffusion_iso_1d_lax_friedrichs_wave_family wave_family;
  advection_diffusion_iso_1d_lax_friedrichs_x_wave_family(C_L, C_R, U_L, U_R, P, &wave_family);

  advection_diffusion_iso_1d_lax_friedrichs_speed_family speed_family;
  advection_diffusion_iso_1d_lax_friedrichs_x_speed_family(C_L, C_R, U_L, U_R, P, &speed_family);

  left_fluctuation->fluctuation_f = 
    (wave_family.wave1.wave_f * fmin(speed_family.speed1, 0.0)) + 
    (wave_family.wave2.wave_f * fmin(speed_family.speed2, 0.0));
}

static inline void
advection_diffusion_iso_1d_lax_friedrichs_x_right_fluctuation(const advection_diffusion_iso_1d_coordinates *C_L, const advection_diffusion_iso_1d_coordinates *C_R,
  const advection_diffusion_iso_1d_state *U_L, const advection_diffusion_iso_1d_state *U_R,
  const advection_diffusion_iso_1d_parameters *P, advection_diffusion_iso_1d_lax_friedrichs_fluctuation *right_fluctuation)
{
  advection_diffusion_iso_1d_lax_friedrichs_wave_family wave_family;
  advection_diffusion_iso_1d_lax_friedrichs_x_wave_family(C_L, C_R, U_L, U_R, P, &wave_family);

  advection_diffusion_iso_1d_lax_friedrichs_speed_family speed_family;
  advection_diffusion_iso_1d_lax_friedrichs_x_speed_family(C_L, C_R, U_L, U_R, P, &speed_family);

  right_fluctuation->fluctuation_f = 
    (wave_family.wave1.wave_f * fmax(speed_family.speed1, 0.0)) + 
    (wave_family.wave2.wave_f * fmax(speed_family.speed2, 0.0));
}

static inline advection_diffusion_iso_1d_lax_friedrichs_fluctuations_1D
create_left_fluctuations()
{
  advection_diffusion_iso_1d_lax_friedrichs_fluctuations_1D left_fluctuations;

  left_fluctuations.x_fluctuation = advection_diffusion_iso_1d_lax_friedrichs_x_left_fluctuation;

  return left_fluctuations;
}

static inline advection_diffusion_iso_1d_lax_friedrichs_fluctuations_1D
create_right_fluctuations()
{
  advection_diffusion_iso_1d_lax_friedrichs_fluctuations_1D right_fluctuations;

  right_fluctuations.x_fluctuation = advection_diffusion_iso_1d_lax_friedrichs_x_right_fluctuation;

  return right_fluctuations;
}

static inline bool
advection_diffusion_iso_1d_lax_friedrichs_x_left_fluctuations_consistent(const advection_diffusion_iso_1d_coordinates *C, const advection_diffusion_iso_1d_state *U,
  const advection_diffusion_iso_1d_parameters *P)
{
  advection_diffusion_iso_1d_lax_friedrichs_fluctuation left_fluctuation;
  advection_diffusion_iso_1d_lax_friedrichs_x_left_fluctuation(C, C, U, U, P, &left_fluctuation);

  double f_left_fluctuation = left_fluctuation.fluctuation_f;

  return ((fabs(f_left_fluctuation) < 1.0e-8));
}

static inline bool
advection_diffusion_iso_1d_lax_friedrichs_x_right_fluctuations_consistent(const advection_diffusion_iso_1d_coordinates *C, const advection_diffusion_iso_1d_state *U,
  const advection_diffusion_iso_1d_parameters *P)
{
  advection_diffusion_iso_1d_lax_friedrichs_fluctuation right_fluctuation;
  advection_diffusion_iso_1d_lax_friedrichs_x_right_fluctuation(C, C, U, U, P, &right_fluctuation);

  double f_right_fluctuation = right_fluctuation.fluctuation_f;

  return ((fabs(f_right_fluctuation) < 1.0e-8));
}

static inline bool
advection_diffusion_iso_1d_lax_friedrichs_x_fluctuations_valid(const advection_diffusion_iso_1d_coordinates *C_L, const advection_diffusion_iso_1d_coordinates *C_R,
  const advection_diffusion_iso_1d_state *U_L, const advection_diffusion_iso_1d_state *U_R, const advection_diffusion_iso_1d_parameters *P)
{
  advection_diffusion_iso_1d_flux left_flux;
  advection_diffusion_iso_1d_x_flux(C_L, P, U_L, &left_flux);

  advection_diffusion_iso_1d_flux right_flux;
  advection_diffusion_iso_1d_x_flux(C_R, P, U_R, &right_flux);

  double f_flux_jump = right_flux.flux_f - left_flux.flux_f;

  advection_diffusion_iso_1d_lax_friedrichs_fluctuation left_fluctuation;
  advection_diffusion_iso_1d_lax_friedrichs_x_left_fluctuation(C_L, C_R, U_L, U_R, P, &left_fluctuation);

  advection_diffusion_iso_1d_lax_friedrichs_fluctuation right_fluctuation;
  advection_diffusion_iso_1d_lax_friedrichs_x_right_fluctuation(C_L, C_R, U_L, U_R, P, &right_fluctuation);

  double f_fluctuation_sum = left_fluctuation.fluctuation_f + right_fluctuation.fluctuation_f;

  return (fabs(f_flux_jump - f_fluctuation_sum) < 1.0e-8);
}

typedef struct {
  double reconstruction_f;
} advection_diffusion_iso_1d_minmod_reconstruction;

typedef struct {
  void (*x_reconstruction)(const advection_diffusion_iso_1d_coordinates *C_L, const advection_diffusion_iso_1d_coordinates *C, const advection_diffusion_iso_1d_coordinates *C_R,
    const advection_diffusion_iso_1d_state *U_L, const advection_diffusion_iso_1d_state *U, const advection_diffusion_iso_1d_state *U_R,
    const advection_diffusion_iso_1d_parameters *P, advection_diffusion_iso_1d_minmod_reconstruction *reconstruction);
} advection_diffusion_iso_1d_minmod_reconstruction_1D;

static inline void
advection_diffusion_iso_1d_minmod_x_left_reconstruction(const advection_diffusion_iso_1d_coordinates *C_L, const advection_diffusion_iso_1d_coordinates *C,
  const advection_diffusion_iso_1d_coordinates *C_R, const advection_diffusion_iso_1d_state *U_L, const advection_diffusion_iso_1d_state *U, const advection_diffusion_iso_1d_state *U_R,
  const advection_diffusion_iso_1d_parameters *P, advection_diffusion_iso_1d_minmod_reconstruction *left_reconstruction)
{
  double x_L = C_L->x;
  double x = C->x; 
  double x_R = C_R->x;
  double f_L = U_L->f;
  double f = U->f;
  double f_R = U_R->f;
  double a = P->a;
  double D = P->D;

  left_reconstruction->reconstruction_f = (f - (0.5 * (fmax(0.0, fmin((f - f_L), (f_R - f))) + fmin(0.0, fmax((f - f_L),
    (f_R - f))))));
}

static inline void
advection_diffusion_iso_1d_minmod_x_right_reconstruction(const advection_diffusion_iso_1d_coordinates *C_L, const advection_diffusion_iso_1d_coordinates *C,
  const advection_diffusion_iso_1d_coordinates *C_R, const advection_diffusion_iso_1d_state *U_L, const advection_diffusion_iso_1d_state *U, const advection_diffusion_iso_1d_state *U_R,
  const advection_diffusion_iso_1d_parameters *P, advection_diffusion_iso_1d_minmod_reconstruction *right_reconstruction)
{
  double x_L = C_L->x;
  double x = C->x; 
  double x_R = C_R->x;
  double f_L = U_L->f;
  double f = U->f;
  double f_R = U_R->f;
  double a = P->a;
  double D = P->D;

  right_reconstruction->reconstruction_f = (f + (0.5 * (fmax(0.0, fmin((f - f_L), (f_R - f))) + fmin(0.0, fmax((f -
    f_L), (f_R - f))))));
}

static inline advection_diffusion_iso_1d_minmod_reconstruction_1D
create_left_reconstruction()
{
  advection_diffusion_iso_1d_minmod_reconstruction_1D left_reconstruction;

  left_reconstruction.x_reconstruction = advection_diffusion_iso_1d_minmod_x_left_reconstruction;

  return left_reconstruction;
}

static inline advection_diffusion_iso_1d_minmod_reconstruction_1D
create_right_reconstruction()
{
  advection_diffusion_iso_1d_minmod_reconstruction_1D right_reconstruction;

  right_reconstruction.x_reconstruction = advection_diffusion_iso_1d_minmod_x_right_reconstruction;

  return right_reconstruction;
}

static inline bool
advection_diffusion_iso_1d_minmod_x_left_reconstruction_consistent(const advection_diffusion_iso_1d_coordinates *C, const advection_diffusion_iso_1d_state *U,
const advection_diffusion_iso_1d_parameters *P)
{
  advection_diffusion_iso_1d_minmod_reconstruction left_reconstruction;
  advection_diffusion_iso_1d_minmod_x_left_reconstruction(C, C, C, U, U, U, P, &left_reconstruction);

  double f_left_reconstruction = left_reconstruction.reconstruction_f;
  double f = U->f;

  return ((fabs(f_left_reconstruction - f) < 1.0e-8));
}

static inline bool
advection_diffusion_iso_1d_minmod_x_right_reconstruction_consistent(const advection_diffusion_iso_1d_coordinates *C, const advection_diffusion_iso_1d_state *U,
const advection_diffusion_iso_1d_parameters *P)
{
  advection_diffusion_iso_1d_minmod_reconstruction right_reconstruction;
  advection_diffusion_iso_1d_minmod_x_right_reconstruction(C, C, C, U, U, U, P, &right_reconstruction);

  double f_right_reconstruction = right_reconstruction.reconstruction_f;
  double f = U->f;

  return ((fabs(f_right_reconstruction - f) < 1.0e-8));
}

typedef struct {
  double interface_f;
} advection_diffusion_iso_1d_centered_diffusion_interface_state;

typedef struct {
  void (*x_interface_state)(const advection_diffusion_iso_1d_coordinates *C_Lx, const advection_diffusion_iso_1d_coordinates *C_Rx,
    const advection_diffusion_iso_1d_state *U_Lx, const advection_diffusion_iso_1d_state *U_Rx,
    const advection_diffusion_iso_1d_parameters *P, advection_diffusion_iso_1d_centered_diffusion_interface_state *interface_state);
} advection_diffusion_iso_1d_centered_diffusion_interface_states_1D;

static inline void
advection_diffusion_iso_1d_centered_diffusion_x_interface_state(const advection_diffusion_iso_1d_coordinates *C_Lx, const advection_diffusion_iso_1d_coordinates *C_Rx,
  const advection_diffusion_iso_1d_state *U_Lx, const advection_diffusion_iso_1d_state *U_Rx,
  const advection_diffusion_iso_1d_parameters *P, advection_diffusion_iso_1d_centered_diffusion_interface_state *interface_state)
{
  double x_Lx = C_Lx->x;
  double x_Rx = C_Rx->x;
  double f_Lx = U_Lx->f;
  double f_Rx = U_Rx->f;
  double a = P->a;
  double D = P->D;

  interface_state->interface_f = (0.5 * (f_Lx + f_Rx));
}

static inline advection_diffusion_iso_1d_centered_diffusion_interface_states_1D
create_interface_states()
{
  advection_diffusion_iso_1d_centered_diffusion_interface_states_1D interface_states;

  interface_states.x_interface_state = advection_diffusion_iso_1d_centered_diffusion_x_interface_state;
  
  return interface_states;
}

typedef struct {
  double interface_f_x;
} advection_diffusion_iso_1d_centered_diffusion_interface_gradient;

typedef struct {
  void (*x_interface_gradient)(const advection_diffusion_iso_1d_coordinates *C_Lx, const advection_diffusion_iso_1d_coordinates *C_Rx,
    const advection_diffusion_iso_1d_state *U_Lx, const advection_diffusion_iso_1d_state *U_Rx,
    const advection_diffusion_iso_1d_parameters *P, advection_diffusion_iso_1d_centered_diffusion_interface_gradient *interface_gradient);
} advection_diffusion_iso_1d_centered_diffusion_interface_gradients_1D;

static inline void
advection_diffusion_iso_1d_centered_diffusion_x_interface_gradient(const advection_diffusion_iso_1d_coordinates *C_Lx, const advection_diffusion_iso_1d_coordinates *C_Rx,
  const advection_diffusion_iso_1d_state *U_Lx, const advection_diffusion_iso_1d_state *U_Rx,
  const advection_diffusion_iso_1d_parameters *P, advection_diffusion_iso_1d_centered_diffusion_interface_gradient *interface_gradient)
{
  double x_Lx = C_Lx->x;
  double x_Rx = C_Rx->x;

  double f_Lx = U_Lx->f;
  double f_Rx = U_Rx->f;

  double a = P->a;
  double D = P->D;

  interface_gradient->interface_f_x = ((f_Rx - f_Lx) / (x_Rx - x_Lx));
}

static inline advection_diffusion_iso_1d_centered_diffusion_interface_gradients_1D
create_interface_gradients()
{
  advection_diffusion_iso_1d_centered_diffusion_interface_gradients_1D interface_gradients;

  interface_gradients.x_interface_gradient = advection_diffusion_iso_1d_centered_diffusion_x_interface_gradient;
  
  return interface_gradients;
}

int
main(int argc, char *argv[])
{
  // Insert simulation drivers here.
  
  return 0;
}
                