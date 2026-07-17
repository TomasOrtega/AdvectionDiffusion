
#lang racket

(define linear-advection-3d
  (hash
   'name "linear-advection-3d"
   'coordinates (list `x `y `z)
   'state (list `f)
   'state-assumptions (list)
   'parameters (list `ax `ay `az)
   'parameters-assumptions (list)
   'fluxes (list
            (list `(* ax f))
            (list `(* ay f))
            (list `(* az f)))
   'wavespeeds (list
                (list `ax)
                (list `ay)
                (list `az))
   'diffusive-fluxes (list
                      (list 0.0)
                      (list 0.0)
                      (list 0.0))
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
            