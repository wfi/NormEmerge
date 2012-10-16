#lang racket

(require games/cards)

(provide (all-defined-out))

;;--- Basic Utilities

;; average : (listof number) -> number
(define (average l)
  (/ (foldr + 0 l) (length l)))

;; stdv : (listof number) -> number
(define (stdv l)
  (sqrt (variance l)))

;; variance: (listof number) -> number
(define (variance l)
  (let ((mu (average l)))
    (/ (foldr + 0 (map (lambda (a) (sqr (- mu a))) l))
       (sub1 (length l)))))

;;--- Population Utilities

;; make-rand-pop: number number number -> population
;; create a population having n agents centered around a mean with a spread scale * [-0.05,0.05]
;; thus, acceptable ranges for spread are 0 (no variation at all) up to 10 with the proviso that the mean is 0.5.
(define (make-rand-pop n mean spread)
  (build-list n 
              (lambda (_) (max 0 (min 1 (+ mean (* spread (- (/ (random 1000) 10000.0) 0.05))))))))


;; pair-up : population -> (listof (list agent agent))
;; create list of random pairs of the agents in the population
(define (pair-up pop)
  (local ((define rand-indices (shuffle-list (build-list (length pop) (lambda (x) x)) 7))
          ;; build-em : (listof N) -> (listof (pair agent agent))
          (define (build-em loi)
            (cond [(empty? loi) empty]
                  [else (cons (list (list-ref pop (first loi))
                                    (list-ref pop (second loi)))
                              (build-em (drop loi 2)))]))
            )
    (build-em rand-indices)
    ))

;; pibar : agent number number -> number
;; given an agent, a, the average bias of a population including agent a, pbar, and the population size, n, 
;; compute the average bias not counting the given agent
(define (pibar a pbar n)
  (/ (- (* n pbar) a)
     (sub1 n)))


