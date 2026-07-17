
#lang racket

(define linear-advection-1d
  (hash
   'name "linear-advection-1d"
   'coordinates (list `x)
   'state (list `f)
   'state-assumptions (list)
   'parameters (list `a)
   'parameters-assumptions (list)
   'fluxes (list
            (list `(* a f)))
   'wavespeeds (list
                (list `a))
   'diffusive-fluxes (list
                      (list 0.0))
   )
  )

(define lax-friedrichs-1d
  (hash
   'name "lax-friedrichs"
   'waves (list
           (let ([flux_f_L `(* a f_L)]
                 [flux_f_R `(* a f_R)]
                 [amax `(abs a)])
             (list
              (list `(* 0.5 (- (- f_R f_L) (/ (- ,flux_f_R ,flux_f_L) ,amax))))
              (list `(* 0.5 (+ (- f_R f_L) (/ (- ,flux_f_R ,flux_f_L) ,amax))))
              )
             )
           )
   'speeds (list
            (let ([amax `(abs a)])
              (list `(- ,amax) amax)))
   )
  )

(define minmod-1d
  (hash
   'name "minmod"
   'order 2
   'left-reconstruction (list
                         (let ([minmod_f `(+ (max 0.0 (min (- f f_L) (- f_R f))) (min 0.0 (max (- f f_L) (- f_R f))))])
                           (list `(- f (* 0.5 ,minmod_f)))))
   'right-reconstruction (list
                          (let ([minmod_f `(+ (max 0.0 (min (- f f_L) (- f_R f))) (min 0.0 (max (- f f_L) (- f_R f))))])
                            (list `(+ f (* 0.5 ,minmod_f)))))
   )
  )
            