#lang racket

(require "neutils.rkt")
(require "theoryNE.rkt")
(require plot)
(require plot/utils)

;; Norm Emergence
;; Simulation of emergence of norms among populations of agents in 2-choice task

;; an agent is a number

;; a population (pop) is a (listof agent)

(define X 0.01)
(define CONST-INC X)
(define Pp 0.01);; proportion for coordination
(define Pn 0.02);; proportion for conflict


;;--- Agent Update Methods ---------------------

;; generic-adjust: agent boolean boolean -> agent
;; current agent, the action taken (true=a, false=(1-a)) and the outcome of interaction (true=confilct, false=coordinate)
;; yields an adjusted agent according to some method


;; adjust-bias : agent boolean boolean -> agent
;; given the selected act and the outcome (true means conflict), return the constant-incremented or decremented agent
(define (adjust-bias a act outcome-conflict)
  (cond [(or (and act outcome-conflict)
             (and (not act) (not outcome-conflict)))
         (max (- a CONST-INC) 0)]
        [else (min (+ a CONST-INC) 1)]))

;; adj-bias-no-clip: agent boolean boolean -> agent
;; given the selected act and outcome (true means conflict), return the constant-incremented or decremented agent (bias) WITHOUT clipping at 0.0 or 1.0
(define (adj-bias-no-clip a act oc)
  (cond [(or (and act oc) (and (not act) (not oc)))
         (- a CONST-INC)]
        [else (+ a CONST-INC)]))

;; adj-proportional-A: agent boolean boolean -> agent
;; given the act chosen by the given agent and the outcome (true means conflict) of that action for the recent pairing,
;; compute the new bias according to the PROPORTIONAL-A update rule:
;; a_i = a_i + P*(1-a_i) if chose act and coordinate,
;; = a_i - P*(a_i) if chose act and conflict,
;; = a_i + P*(1-a_i) if chose NOT-act and conflict
;; = a_i - P*(a_i) if chose NOT-act and coordinate
;; as of 10/9/2011, this version is consistant with the analysis that simplifies to (p-bar_i-hat - p_i)
(define (adj-proportional-A a act outcome)
  (cond [act (cond [outcome (- a (* P a))]
                   [else (+ a (* P (- 1.0 a)))])]
        [else (cond [outcome (+ a (* P (- 1.0 a)))]
                    [else (- a (* P a))])]))

;; adj-proportional-B: agent boolean boolean -> agent
;; using the other proportional update rule
(define (adj-proportional-B a act outcome)
  (+ a (adj-prop-update-B a act outcome)))

;; adj-prop-update-B: agent boolean boolean -> agent
;; using the other proportional update rule, but return only the update, Delta_i, for this agent
;; a_i =  + P*(1-a_i) if chose act and coordinate
;;     =  - P*(1-a_i) if chose act and conflict
;;     =  + P*a_i if chose NOT-act and conflict
;;     =  - P*a_i if chose NOT-act and coordinate
(define (adj-prop-update-B a act outcome)
  (cond [act (cond [outcome (- (* P (- 1.0 a)))]
                   [else (* P (- 1.0 a))])]
        [else (cond [outcome (* P a)]
                    [else (- (* P a))])]))

;; adj-proportional-C: agent boolean boolean -> agent
;; using the (new) third proportional update rule
;; a_i =  + 0.5-(0.5-a_i)^2 if chose act and coordinate
;;     =  - 0.5-(0.5-a_i)^2 if chose act and conflict
;;     =  + 0.5-(0.5-a_i)^2 if chose NOT-act and conflict
;;     =  - 0.5-(0.5-a_i)^2 if chose NOT-act and coordinate
(define (adj-proportional-C a act outcome)
  (local ([define weight (* P 
                            (- 0.5 (sqr (- 0.5 a))) ;; this isn't what we want as when a is 0 or 1, it still gives 0.25
                            ;;(- 1 (sqr (- (* 2 a) 1))) ;; scaling the ends and middle to 0 and 1 respectively
                            )])
    (cond [act (cond [outcome (- a weight)]
                     [else (+ a weight)])]
          [else (cond [outcome (+ a weight)]
                      [else (- a weight)])])))

;; delta-i-propA2: agent boolean boolean -> agent
;; using variation of propA where adjustment is always function of gap yet-to-go with two proportions
;; delta-i = - Pn * (a_i) if choose act and conflict: decrease a_i by Pn of a_i
;;         = + Pp * (1-a_i) if choose act and coordinate: increase a_i by Pp of (1-a_i)
;;         = + Pn * (1-a_i) if choose NOT-act and conflict: decrease not-act by Pn of not-act, 
;;            or complete update for not act (1-a_i) = (1-a_i) - Pn * (1-a_i) = (1 - Pn)*(1-a_i), 
;;            and a_i =  1 - (1-Pn)*(1-a_i) = 1 - 1 + a_i + Pn - Pn*a_i = a_i + Pn*(1-a_i)
;;         = - Pp * a_i if choose NOT-act and coordinate: increas not-act by Pp of (1-not-act), 
;;            or complete update for not act = (1-a_i) + Pp * (1 - (1-a_i)) = (1-a_i) + Pp*(a_i) = 1 - a_i(1 - Pp) 
(define (delta-i-propA2 a act conflict)
  (cond [(and act conflict) (* -1 Pn a)]
        [(and act (not conflict)) (* Pp (- 1 a))]
        [(and (not act) conflict) (* Pn (- 1 a))]
        [(and (not act) (not conflict)) (* -1 Pp a)]))


;;--- Agent Interactions --------------------

;; interact : agent agent (agent boolean boolean -> agent) -> (agent agent)
;; for two agents, interact and update according to the given update function
(define (interact a1 a2 update-func)
  (local ((define a1-act (<= (random) a1))
          (define a2-act (<= (random) a2))
          (define (xor b1 b2)
            (and (or b1 b2) (not (and b1 b2))))
          (define conflict-outcome (xor a1-act a2-act))
          )
    (list
     (update-func a1 a1-act conflict-outcome)
     (update-func a2 a2-act conflict-outcome))))

;; one-iteration-f: population (number boolean boolean -> nu mber) -> population
;; randomly pair-up the given population, have each pair interact and be updated based on the given adjustf function.
(define (one-iteration-f pop adjustf)
  (foldr (lambda (p r) (append (interact (first p) (second p) adjustf) r)) empty (pair-up pop)))

;; multi-iteration-f : population N (number boolean boolean -> number) -> population
;; do n iterations of one-iteration-f
(define (multi-iteration-f pop n adjustf)
  (cond [(zero? n) pop]
        [else (multi-iteration-f (one-iteration-f pop adjustf) (sub1 n) adjustf)]))


;; build-pop-averages: pop N (number boolean boolean -> number) -> (listof number)
;; construct the sequence of average bias in the population over n iterations
(define (build-pop-averages p n f)
  (cond [(zero? n) (list (average p))]
        [else (let ((new-p (one-iteration-f p f)))
                (cons (average new-p)
                      (build-pop-averages new-p (sub1 n) f)))]))


;; sim-avg-update: population (number boolean boolean -> number) number -> number
;; for the given agent-adjustment update function, compute the avg delta for one interaction and update of the given population
;; optionally, average this over ntimes with the same population but different pairings
(define (sim-avg-update pop adjustf ntimes)
  (average (build-list ntimes
                       (lambda (_) (- (average (one-iteration-f pop adjustf))
                                      (average pop))))))


;; sim-avgstd-delta: population (number boolean boolean -> number) -> (list number number)
;; simulate a paired-up interaction and report the average update (the change in population average)
;; and also the standard deviation of the population at the start
#|
(define (sim-avgstd-delta pop adjustf)
  (list (sim-avg-delta pop adjustf 1)
        (stdv pop)))|#

;; mutli-iteration-plot : population N (number boolean boolean -> number) -> population
(define (multi-iteration-plot pop n adjustf)
  (cond [(zero? n) pop]
        [else 
         
         
         (local ([define new-pop (one-iteration-f pop adjustf)])
           (plot (function (lambda (n) (/ (foldl 0 + new-pop) (length new-pop)))))
           (multi-iteration-plot new-pop (sub1 n) adjustf))]))

;; build-filtered-list
(define (build-filtered-list n proc pred list)
  (cond
    [(zero? n) list]
    [else (local ([define next-element (proc n)])
            (cond
              [(pred next-element) (build-filtered-list (sub1 n) proc pred (cons next-element list))]
              [else (build-filtered-list n proc pred list)]))]))


;; plot-simulation : population N (number boolean boolean -> number) number -> ?
(define (plot-simulation pop n adjustf runs)
  (plot (build-list runs (lambda(_) 
                           (lines 
                            (let* ((data (build-pop-averages pop n adjustf)))
                              (build-list (length data)
                                          (lambda (m) (vector m (list-ref data m)))))
                            #:y-min 0 #:y-max 1)))))

;; plot-filtered-simulation : population N (number boolean boolean -> number) number (list -> boolean) -> ?
(define (plot-filtered-simulation pop n adjustf runs pred)
  (plot (map (lambda (l) (lines l #:y-min 0 #:y-max 1)) (build-filtered-list runs (lambda(_) 
                                                                                    (let* ((data (build-pop-averages pop n adjustf)))
                                                                                      (build-list (length data)
                                                                                                  (lambda (m) (vector m (list-ref data m))))))
                                                                             pred
                                                                             empty))))

;; average-filtered-simulation: population number (number boolean boolean -> number) number ((listof (vector number number)) -> boolean) -> (listof (vector number number))
;; create a plottable list of vector pairs averaged over 'runs' runs satisfying 'pred' starting with a given population out to n timesteps
(define (average-filtered-simulation pop n adjustf runs pred)
(local ([define trials (build-filtered-list runs (lambda(_) 
                              (let* ((data (build-pop-averages pop n adjustf)))
                                (build-list (length data)
                                            (lambda (m) (vector m (list-ref data m))))))
                       pred
                       empty)])
    (build-list n 
              (lambda (point) (vector point (/ (foldl + 0 (map (lambda (trial) (vector-ref (list-ref trial point) 1)) trials)) runs))))))

;; theoretical-model : N number/population -> (listof number)
;; called with the desired size of the list and the starting value of pbar
(define (theoretical-model n p-list)
  (cond
    [(zero? (sub1 n)) (reverse p-list)]
    [(list? p-list) (theoretical-model (sub1 n) (cons (+ (first p-list) (first (avg-pbar-delta-propB (list (first p-list))))) p-list))]
    [else (theoretical-model (sub1 n) (list (+ p-list (first (avg-pbar-delta-propB (list p-list)))) p-list))]))

;; theoretical-model-points: number number -> (listof (vector number number))
(define (theoretical-model-points n pbar0)
  (map (lambda (x y) (vector x y)) (build-list n (lambda (x) x)) (theoretical-model n pbar0)))

;; theoretical-model-with-variance:
(define (theoretical-model-with-variance n pop list)
    (cond
    [(zero? (sub1 n)) (reverse list)]
    [(empty? list) (theoretical-model-with-variance (sub1 n))]
    [else (theoretical-model-with-variance (sub1 n) (list (+ list (first (avg-pbar-delta-propB (list list)))) list))]))

;; plot-variances-over-time-with-comparison
(define (plot-variances-over-time pop n adjustf theoretical-adjustf)
  (plot (list
         (lines (map (lambda (x y) (vector x y)) (build-list n (lambda (m) m)) (variances-over-time pop n empty adjustf)) #:y-min 0 #:y-max 0.1)
         (lines (map (lambda (x y) (vector x y)) (build-list n (lambda (m) m)) (theoretical-variances-over-time pop n empty theoretical-adjustf)) #:y-min 0 #:y-max 0.1)))) 
 
;; theoretical-variances-over-time : (listof number) number (list of number) (agent number -> number) -> (list of number)
(define (theoretical-variances-over-time pop n variances adjustf)
  (cond
    [(zero? n) (reverse variances)]
    [else (theoretical-variances-over-time (map + pop (adjustf pop)) (sub1 n) (cons (variance pop) variances) adjustf)]))

;; variances-over-time : (listof number) number (list of number) (agent number -> number) -> (list of number)
(define (variances-over-time pop n variances adjustf)
  (cond [(zero? n) (reverse variances)]
        [else (variances-over-time (one-iteration-f pop adjustf) (sub1 n) (cons (variance pop) variances) adjustf)]))

;; compute-errors-propC: population number -> (pair number number)
;; For a given a population, compute (1) observed pbar from n match-ups, (2) expected pbar using theoretical delta_i, and
;; (3) expectd pbar using estimated expression based on pbar and population variance.  Return the errors between 1 and 2, and between 2 and 3.
(define (compute-errors-propC pop n)
  (local ([define psize (length pop)]
          [define pbar (/ (foldl + 0 pop) psize)]
          [define pbar1 (/ (foldl + 0 (build-list n (lambda (_) (/ (foldl + 0 (one-iteration-f pop adj-proportional-C)) psize)))) n)]
          [define pbar2 (/ (foldl + 0 (map + pop (avg-pbar-delta-propC pop))) psize)]
          [define pbar3 (+ pbar (delta-bar-propC pbar psize (variance pop)))])
            (cons (- pbar1 pbar2) (- pbar1 pbar3))))

;; plot-errors-propC: (listof population) number ((pairof number) -> number) -> image
;; plots one of the errors (determined by selector) calculated by compute-errors-propC for the given population
(define (plot-errors-propC pops n selector)
  (local ([define errors (map (lambda (pop) (compute-errors-propC pop n)) pops)])
    (plot3d (list (points3d (map (lambda (pop error) (vector (average pop) (variance pop) (selector error))) pops errors))))))

;(plot-filtered-simulation (build-list 20 (lambda (_) 0.49)) 150 adj-proportional-B 100 (lambda (l) (< (vector-ref (last l) 1) 0.49)))

#|(plot (list 
       (lines (average-filtered-simulation (build-list 20 (lambda (_) 0.49)) 200 adj-proportional-B 100 (lambda (l) (< (vector-ref (last l) 1) 0.3))) #:y-min 0 #:y-max 1)
       (lines (theoretical-model-points 200 0.49) #:y-min 0 #:y-max 1)))|#

;(plot-variances-over-time (make-rand-pop 20 0.49 2) 150 adj-proportional-C avg-pbar-delta-propC)

(define a-pop (build-list 50 (lambda (_) (make-beta-pop 20 0.5 .05))))
(plot-errors-propC a-pop 100 car)
(plot-errors-propC a-pop 100 cdr)