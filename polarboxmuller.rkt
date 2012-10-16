#lang racket

;; polar-box-muller: number number -> (  -> number)
;; this isn't finished.  Should take a mean and stdv and generate random numbers
;; sampled from the corresponding normal distributions.
;; This returns the function that generates one such number each time it is called.
(define (polar-box-muller mean stdv)
  (local ((define y1 0)
          (define y2 0)
          (define counter 1)
          
          (define (do-2-more x1 x2 w)
            (cond [(>= w 1.0) 
                   (let ((fw (sqrt (/ (* -2.0 (log w)) w))))
                     (set! y1 (* x1 fw))
                     (set! y2 (* x2 fw)))]
                  [else 
                   (let ((lx1 (- (* 2.0 (random)) 1.0))
                         (lx2 (- (* 2.0 (random)) 1.0)))
                     (do-2-more lx1 lx2 (+ (sqr lx1) (sqr lx2))))]))
          
          (define (get-one)
            (set! counter (modulo (add1 counter) 2))
            (cond [(zero? counter)
                   (do-2-more 0 0 0)
                   y1]
                  [else y2]))
          )
    get-one))


           
                  