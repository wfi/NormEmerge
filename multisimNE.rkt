#lang racket

(require "neutils.rkt")

(define P 0.05)
(define epsilon 0.00001)

;; Norm Emergence
;; Simulation of emergence of norms among populations of agents in N-choice task

;; an agent is a (vectorof number)
;; where the vector has as many cells as there are possible actions, the values are probabilities and must sum to 1.0

;; a population (pop) is a (listof agent)


;;--- local utilities ---------------

;; vfold: (Y X -> X) X (vectorof Y) -> (vectorof X)
;; fold for vectors
(define (vfold f b v)
  (local ((define (help i sofar)
            (cond [(= i (vector-length v)) sofar]
                  [else (help (add1 i) (f (vector-ref v i) sofar))])))
    (help 0 b)))

;; normalize-agent: agent -> agent
;; ensure that the agent's probabilities sum to 1.0 and therefore may represent an agent
(define (normalize-agent a)
  (let* ((tsum (vfold + 0 a)))
    (list->vector (map (lambda (p) (/ p tsum)) tsum))))
;; normalize-agent!: agent -> agent
;; destructive version of same
(define (normalize-agent! a)
  (let* ((tsum (vfold + 0 a)))
    (vector-map! (lambda (p) (/ p tsum)) a)))

;; make-rand-agent: number -> agent
;; make an agent of the given size with random probabilities for each of the asize possible actions
(define (make-rand-agent asize)
  (normalize-agent! (build-vector asize (lambda (_) (random)))))

;; agent-act: agent -> number
;; probabilistically select a random action based on the agent's current biases
(define (agent-act a)
  (local ((define (pick-one i val)
            (cond [(<= val (vector-ref a i)) i]
                  [else (pick-one (add1 i) (- val (vector-ref a i)))])))
    (pick-one 0 (random))))

;; interact-X: agent agent (agent number (1 OR -1) -> (listof number)) -> (list agent agent)
;; simulate the interaction of two agents and return their respective updates for the given delta-operator
(define (interact-X a1 a2 delta-op)
  (let ((a1-act (agent-act a1))
        (a2-act (agent-act a2)))
    (cond [(= a1-act a2-act)
           (list (delta-op a1 a1-act 1) (delta-op a2 a2-act 1))]
          [else (list (delta-op a1 a1-act -1) (delta-op a2 a2-act -1))])))
                      
          

;; deltajsA: agent number (1 OR -1) -> (listof number)
;; compute delta for each action of agent a having taken action act based on a uniform redstribution
(define (deltajsA a act outcome-sign)
  (local ((define (updates i)
            (cond [(= i (vector-length a)) empty]
                  [(= i act) (cons (* outcome-sign P (- 1.0 (vector-ref a i)))
                                   (updates (add1 i)))]
                  [else (cons (/ (* outcome-sign -1 P (- 1.0 (vector-ref a act)))
                                 (sub1 (vector-length a)))
                              (updates (add1 i)))])))
    (updates 0)))

;; deltajsB: agent number (1 OR -1) -> (listof number)
;; compute delta for each action of agent a having taken action act based on a probability weighted redistribution
;; The outcome-sign should be 1 when agent coordinated and -1 when conflict.
(define (deltajsB a act outcome-sign)
  (local ((define (updates i)
            (cond [(= i (vector-length a)) empty]
                  [(= i act) (cons (* outcome-sign P (+ (- 1.0 (vector-ref a i)) epsilon)) ;; this and following epsilon treatment is not quite right
                                   (updates (add1 i)))]
                  [else (cons (* outcome-sign -1 P (+ (vector-ref a i) epsilon) (+ (vector-ref a act) epsilon))
                              (updates (add1 i)))])))
    (updates 0)))

;; agent-update!: agent (listof number) -> agent
;; given the list of updates from interact, destructively modify the agent accordingly
(define (agent-update! a adj)
  (normalize-agent! 
   (vector-map! (lambda (p delta-j)
                  (if (< delta-j 0)
                      (max (+ p delta-j) 0)
                      (min (+ p delta-j) 1))) 
                a
                (list->vector adj))))


;;--- do simple tests -------------------

(define 4pop (build-list 10 (lambda (_) (make-rand-agent 4))))
(define 8pop (build-list 10 (lambda (_) (make-rand-agent 8))))
(define 16pop (build-list 10 (lambda (_) (make-rand-agent 16))))

;; one-pop-step: population (agent number (1 OR -1) -> (listof number) ) -> (list number)
;; destructively update the population after pairing them up but then return a sum of each agent's probabilities as check on consistency
(define (one-pop-step pop deltajs-op)
  (map (lambda (pair)
         (let* ((res (interact-X (first pair) (second pair) deltajs-op))
                )
           (agent-update! (first pair) (first res))
           (agent-update! (second pair) (second res))))
       (pair-up pop))
  (map (lambda (a) (vfold + 0 a)) pop))

;; n-pop-steps: population number (agent number (1 OR -1) -> population
;; run one-pop-step n times returning the population when done [everything is destructive]
(define (n-pop-steps pop n deltajs-op)
  (do ((i n (sub1 i)))
    ((zero? i) pop)
    (one-pop-step pop deltajs-op)))
