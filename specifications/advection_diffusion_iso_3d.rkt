
#lang racket

(define advection-diffusion-iso-3d
  (hash
   'name "advection-diffusion-iso-3d"
   'coordinates (list `x `y `z)
   'state (list `f)
   'state-assumptions (list)
   'parameters (list `ax `ay `az `D)
   'parameters-assumptions (list `(> D 0.0))
   'fluxes (list
            (list `(* ax f))
            (list `(* ay f))
            (list `(* az f)))
   'wavespeeds (list
                (list `ax)
                (list `ay)
                (list `az))
   'diffusive-fluxes (list
                      (list `(* D f_x))
                      (list `(* D f_y))
                      (list `(* D f_z)))
   )
  )

(define lax-friedrichs-3d
  (hash
   'name "lax-friedrichs"
   'waves (list
           (let ([flux_f_L `(* ax f_L)]
                 [flux_f_R `(* ax f_R)]
                 [amax `(abs ax)])
             (list
              (list `(* 0.5 (- (- f_R f_L) (/ (- ,flux_f_R ,flux_f_L) ,amax))))
              (list `(* 0.5 (+ (- f_R f_L) (/ (- ,flux_f_R ,flux_f_L) ,amax))))
              )
             )
           (let ([flux_f_L `(* ay f_L)]
                 [flux_f_R `(* ay f_R)]
                 [amax `(abs ay)])
             (list
              (list `(* 0.5 (- (- f_R f_L) (/ (- ,flux_f_R ,flux_f_L) ,amax))))
              (list `(* 0.5 (+ (- f_R f_L) (/ (- ,flux_f_R ,flux_f_L) ,amax))))
              )
             )
           (let ([flux_f_L `(* az f_L)]
                 [flux_f_R `(* az f_R)]
                 [amax `(abs az)])
             (list
              (list `(* 0.5 (- (- f_R f_L) (/ (- ,flux_f_R ,flux_f_L) ,amax))))
              (list `(* 0.5 (+ (- f_R f_L) (/ (- ,flux_f_R ,flux_f_L) ,amax))))
              )
             )
           )
   'speeds (list
            (let ([amax `(abs ax)])
              (list `(- ,amax) amax))
            (let ([amax `(abs ay)])
              (list `(- ,amax) amax))
            (let ([amax `(abs az)])
              (list `(- ,amax) amax)))
   )
  )

(define minmod-3d
  (hash
   'name "minmod"
   'order 2
   'left-reconstruction (list
                         (let ([minmod_f `(+ (max 0.0 (min (- f f_L) (- f_R f))) (min 0.0 (max (- f f_L) (- f_R f))))])
                           (list `(- f (* 0.5 ,minmod_f))))
                         (let ([minmod_f `(+ (max 0.0 (min (- f f_L) (- f_R f))) (min 0.0 (max (- f f_L) (- f_R f))))])
                           (list `(- f (* 0.5 ,minmod_f))))
                         (let ([minmod_f `(+ (max 0.0 (min (- f f_L) (- f_R f))) (min 0.0 (max (- f f_L) (- f_R f))))])
                           (list `(- f (* 0.5 ,minmod_f)))))
   'right-reconstruction (list
                          (let ([minmod_f `(+ (max 0.0 (min (- f f_L) (- f_R f))) (min 0.0 (max (- f f_L) (- f_R f))))])
                            (list `(+ f (* 0.5 ,minmod_f))))
                          (let ([minmod_f `(+ (max 0.0 (min (- f f_L) (- f_R f))) (min 0.0 (max (- f f_L) (- f_R f))))])
                            (list `(+ f (* 0.5 ,minmod_f))))
                          (let ([minmod_f `(+ (max 0.0 (min (- f f_L) (- f_R f))) (min 0.0 (max (- f f_L) (- f_R f))))])
                            (list `(+ f (* 0.5 ,minmod_f)))))
   )
  )

(define centered-diffusion-3d
  (hash
   'name "centered-diffusion"
   'interface-state (list
                     (list `(* 0.5 (+ f_Lx f_Rx)))
                     (list `(* 0.5 (+ f_Ly f_Ry)))
                     (list `(* 0.5 (+ f_Lz f_Rz))))
   'interface-gradients (list
                         (list
                          (list `(/ (- f_Rx f_Lx) (- x_Rx x_Lx)))
                          (list `(* 0.5 (+ (/ (- f_Lx_Ry f_Lx_Ly) (- y_Lx_Ry y_Lx_Ly))
                                           (/ (- f_Rx_Ry f_Rx_Ly) (- y_Rx_Ry y_Rx_Ly)))))
                          (list `(* 0.5 (+ (/ (- f_Lx_Rz f_Lx_Lz) (- z_Lx_Rz z_Lx_Lz))
                                           (/ (- f_Rx_Rz f_Rx_Lz) (- z_Rx_Rz z_Rx_Lz))))))
                         (list
                          (list `(* 0.5 (+ (/ (- f_Ly_Rx f_Ly_Lx) (- x_Ly_Rx x_Ly_Lx))
                                           (/ (- f_Ry_Rx f_Ry_Lx) (- x_Ry_Rx x_Ry_Lx)))))
                          (list `(/ (- f_Ry f_Ly) (- y_Ry y_Ly)))
                          (list `(* 0.5 (+ (/ (- f_Ly_Rz f_Ly_Lz) (- z_Ly_Rz z_Ly_Lz))
                                           (/ (- f_Ry_Rz f_Ry_Lz) (- z_Ry_Rz z_Ry_Lz))))))
                         (list
                          (list `(* 0.5 (+ (/ (- f_Lz_Rx f_Lz_Lx) (- x_Lz_Rx x_Lz_Lx))
                                           (/ (- f_Rz_Rx f_Rz_Lx) (- x_Rz_Rx x_Rz_Lx)))))
                          (list `(* 0.5 (+ (/ (- f_Lz_Ry f_Lz_Ly) (- y_Lz_Ry y_Lz_Ly))
                                           (/ (- f_Rz_Ry f_Rz_Ly) (- y_Rz_Ry y_Rz_Ly)))))
                          (list `(/ (- f_Rz f_Lz) (- z_Rz z_Lz)))))
   )
  )
            