
#lang racket

(define advection-diffusion-iso-2d
  (hash
   'name "advection-diffusion-iso-2d"
   'coordinates (list `x `y)
   'state (list `f)
   'state-assumptions (list)
   'parameters (list `ax `ay `D)
   'parameters-assumptions (list `(> D 0.0))
   'fluxes (list
            (list `(* ax f))
            (list `(* ay f)))
   'wavespeeds (list
                (list `ax)
                (list `ay))
   'diffusive-fluxes (list
                      (list `(* D f_x))
                      (list `(* D f_y)))
   )
  )

(define lax-friedrichs-2d
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
           )
   'speeds (list
            (let ([amax `(abs ax)])
              (list `(- ,amax) amax))
            (let ([amax `(abs ay)])
              (list `(- ,amax) amax)))
   )
  )

(define minmod-2d
  (hash
   'name "minmod"
   'order 2
   'left-reconstruction (list
                         (let ([minmod_f `(+ (max 0.0 (min (- f f_L) (- f_R f))) (min 0.0 (max (- f f_L) (- f_R f))))])
                           (list `(- f (* 0.5 ,minmod_f))))
                         (let ([minmod_f `(+ (max 0.0 (min (- f f_L) (- f_R f))) (min 0.0 (max (- f f_L) (- f_R f))))])
                           (list `(- f (* 0.5 ,minmod_f)))))
   'right-reconstruction (list
                          (let ([minmod_f `(+ (max 0.0 (min (- f f_L) (- f_R f))) (min 0.0 (max (- f f_L) (- f_R f))))])
                            (list `(+ f (* 0.5 ,minmod_f))))
                          (let ([minmod_f `(+ (max 0.0 (min (- f f_L) (- f_R f))) (min 0.0 (max (- f f_L) (- f_R f))))])
                            (list `(+ f (* 0.5 ,minmod_f)))))
   )
  )

(define centered-diffusion-2d
  (hash
   'name "centered-diffusion"
   'interface-state (list
                     (list `(* 0.5 (+ f_Lx f_Rx)))
                     (list `(* 0.5 (+ f_Ly f_Ry))))
   'interface-gradients (list
                         (list
                          (list `(/ (- f_Rx f_Lx) (- x_Rx x_Lx)))
                          (list `(* 0.5 (+ (/ (- f_Lx_Ry f_Lx_Ly) (- y_Lx_Ry y_Lx_Ly))
                                           (/ (- f_Rx_Ry f_Rx_Ly) (- y_Rx_Ry y_Rx_Ly))))))
                         (list
                          (list `(* 0.5 (+ (/ (- f_Ly_Rx f_Ly_Lx) (- x_Ly_Rx x_Ly_Lx))
                                           (/ (- f_Ry_Rx f_Ry_Lx) (- x_Ry_Rx x_Ry_Lx)))))
                          (list `(/ (- f_Ry f_Ly) (- y_Ry y_Ly)))))
   )
  )
            