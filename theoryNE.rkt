#lang racket

(require "neutils.rkt")

(provide (all-defined-out))

;; Theory code for Norm Emergence Analysis
;; Two-choice task

(define P 0.05) ;; proportional update factor


;;--- compute expected Delta-i's and Delta-bars based either pibar or on the estimate, pbar

;; delta-i-p{i}bar-propA : agent number number -> number
;; for the propA method (pbar - pi) and (pibar - pi) respectively
(define (delta-i-pbar-propA a pbar psize)
  (- pbar a))
(define (delta-i-pibar-propA a pbar psize)
  (- (pibar a pbar psize) a))

;; delta-i-pibar-propB : agent number number -> number
;; compute the increment/decrement to the bias of this agent based on the pibar
(define (delta-i-pibar-propB a pbar psize)
  (* 2 P (- a (sqr a)) (- (* 2 (pibar a pbar psize)) 1.0)))

;; delta-i-pbar-propB : agent number number -> number
;; compute the increment/decrement to the bias of this agent based on the ESTIMATE pbar for actual pibar
(define (delta-i-pbar-propB a pbar psize)
  (* 2 P (- a (sqr a)) (- (* 2 pbar) 1.0)))


;; avg-pbar-delta-propA: population -> number
;; using the PROPORTIONAL-A method, 1/n * (P * summate( pbar - a )), compute the average delta based on the given population values.
;; NOTE, this IS making the pi_bar = pbar estimation
(define (avg-pbar-delta-propA pop)
  (delta-bar pop delta-i-pbar-propA))

;; avg-pibar-delta-propA: population -> number
;; this is not making the pi_bar = pbar estimation
(define (avg-pibar-delta-propA pop)
  (delta-bar pop delta-i-pibar-propA))

;; avg-pbar-delta-propB: population -> (listof number)
;; using the theory derived from the other proportional update rule, but return only the update value, Delta_i, for this agent
;; based on the update rule:
;; a_i =  + P*(1-a_i) if chose act and coordinate
;;     =  - P*(1-a_i) if chose act and conflict
;;     =  + P*a_i if chose NOT-act and conflict
;;     =  - P*a_i if chose NOT-act and coordinate
;; new-pop-proportional-B: population -> population
;; compute the new population based on the theory using pbar as the estimate for pibar
(define (avg-pbar-delta-propB pop)
  (delta-bar pop delta-i-pbar-propB))

;; avg-pibar-delta-propB: population -> (listof number)
;; compute the Delta-bar for the given population based on the actual pibar value
(define (avg-pibar-delta-propB pop)
  (delta-bar pop delta-i-pibar-propB))

;; delta-bar: population (agent number -> number) -> population
;; abstract version of the average functions for propA/probB with pbar/pibar
(define (delta-bar pop f-delta-i)
  (let ((pbar (average pop))
        (psize (length pop)))
    (map (lambda (a) (f-delta-i a pbar psize)) pop)))

;; variance-t-plus-one : population 
(define (variance-t-plus-one pop adj-f)
  (local ([define pbar-t-plus-one (/ (foldl + 0 (map adj-f pop)) (length pop))])
    (/ (foldl + 0 (map (lambda (pi) (sqr (- (adj-f pi) pbar-t-plus-one))) pop)) (length pop))))
