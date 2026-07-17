
#lang racket

(define linear-advection-2d
  (hash
   'name "linear-advection-2d"
   'coordinates (list `x `y)
   'state (list `f)
   'state-assumptions (list)
   'parameters (list `ax `ay)
   'parameters-assumptions (list)
   'fluxes (list
            (list `(* ax f))
            (list `(* ay f)))
   'wavespeeds (list
                (list `ax)
                (list `ay))
   'diffusive-fluxes (list
                      (list 0.0)
                      (list 0.0))
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
            