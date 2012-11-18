#lang racket

(require "neutils.rkt")

(define X 0.01)
(define Y 0.01)
(define CLOSE-ENOUGH 0.02)
(define CONST-INC X)  ;; used in simulation but X used in theory
(define P 0.05) ;; Proportional update factor

;; -------------------------------
;; Data Definitions

;; an agent is a number
;; a population is a (listof number)

;;--------------------------------
;; Utilities

;; cross-prod: (listof X) (listof Y) -> (listof (list X Y))
;; generate the cross product of the two lists
(define (cross-prod l1 l2)
  (cond [(empty? l2) empty]
        [else (append (map (lambda (x) (list x (car l2))) l1)
                      (cross-prod l1 (cdr l2)))]))

;; group-by: (listof (list X Y ...)) N -> (listof (listof (list X Y ...)))
;; group a list of tuples into lists of tuples which share the same value as indexed by nth
(define (group-by lot n)
  (local ((define unique-vals 
            (foldr (lambda (tuple sofar) 
                     (if (member (list-ref tuple n) sofar) sofar (cons (list-ref tuple n) sofar))) 
                   empty lot))
          )
    (map (lambda (uval) (filter (lambda (tup) (eq? uval (list-ref tup n))) lot)) unique-vals)))


;; create-pops-plus: number number (listof number) (listof number) -> (listof (list number number population))
;; create a list of augmented populations given as list of triples, where the first is the target mean,
;; the second is the spread used to generate the third is the population
(define (create-pops-plus howmany popsize means sprds)
  (foldr (lambda (mean spread sofar)
           (append (build-list howmany 
                               (lambda (_) (list mean spread (make-rand-pop popsize mean spread))))
                   sofar))
         empty
         (map first (cross-prod means sprds))
         (map second (cross-prod means sprds))))

;;--- simple utilities

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

;;--- compute average deltas based either pibar or on the estimate, pbar

;; pibar : agent number number -> number
;; given an agent, a, the average bias of a population including agent a, pbar, and the population size, n, 
;; compute the average bias not counting the given agent
(define (pibar a pbar n)
  (/ (- (* n pbar) a)
     (sub1 n)))

;; delta-i : agent number -> number
;; compute the increment/decrement to the bias of this agent based on the pibar
(define (delta-i a pibar)
  (+ (* a X a pibar)
     (* -1 (- 1 a) X (- 1 a) (- 1 pibar))
     (* (- 1 a) Y (- 1 a) pibar)
     (* -1 a Y a (- 1 pibar))))

#|
(define (delta-i-reduced a pop)
  (* X (/ 1 (sub1 (length pop)))
     (foldr (lambda (a s) (+ s (- (* 2 a) 1))) 0 (remove a pop))))

;; estimate-delta-i: (listof number) -> number
;; for a given agent a from a given full population (including a), estimate the delta-i for any given agent using pbar
(define (estimate-delta-i pop)
  (- (* X 2 (average pop)) X))
|#

;; estimate-delta-bar: (listof number) -> number
;; estimate the average delta over the given population by using only the average bias (X*(2*pbar - 1))
(define (estimate-delta-bar pop)
  (- (* X 2 (average pop)) X))

;; pbar-error: (listof number) -> number
;; given a population, determine the error between the "true" average-update based on removing each agent (i.e., pibar)
;; as compared to an estimate that ignores this consideration (i.e., pbar)
(define (pbar-error pop)
  (let ((pbar (average pop))
        (n (length pop)))
    (abs (- (average (map (lambda (a) (delta-i a (pibar a pbar n))) pop))
            (estimate-delta-bar pop)))))

;; pbar-error-wrapper: (listof number) (listof number) number number -> (listof (listof number))
;; a testing wrapper for evaluating the pbar-error:  create howmany-of-each populations of the given size with given means and spreads
;; returns listof populations where the first two numbers in the population are mean and stdv
(define (pbar-error-wrapper means spreads popsize howmany-of-each)
  (map (lambda (mu-sprd) 
         (let ((res (build-list howmany-of-each (lambda (_) (let ((nupop (make-rand-pop popsize (first mu-sprd) (second mu-sprd))))
                                                              (list (average nupop) (stdv nupop) (pbar-error nupop)))))))
           (list (first mu-sprd) (second mu-sprd)
                 (average (map first res))
                 (average (map second res))
                 (average (map third res)))))
       (cross-prod means spreads)))

;; many-pop-generate:  (listof number) (listof number) number number -> (listof (listof number))
;; this one actually generates populations that can then be used for an arbitrary purpose
(define (many-pop-generate means spreads popsize howmany-of-each)
  (foldr (lambda (mu-sprd res) 
           (let ((pops (build-list howmany-of-each (lambda (_) (let ((nupop (make-rand-pop popsize (first mu-sprd) (second mu-sprd))))
                                                                (list (first mu-sprd) (second mu-sprd) (average nupop) (stdv nupop) nupop))))))
             (append pops res)))
         empty
         (cross-prod means spreads)))

;; do-pbar-error: (listof (listof number)) string -> ...
;; run and write the pbar-error on each population with the mean and stdv found in the first two elements
(define (do-pbar-error lores fname)
  (with-output-to-file 
      (expand-user-path (string->path (string-append "~/research/NormEmerge/" fname)))
    (lambda ()
      (for-each (lambda (bundle)
                  (printf "~a~%" bundle))
                lores))))

;; do-sim-steps-to-limit: (listof number) (listof population) boolean number string -> ...
;; run each augmented-pop (from many-pops-generate) on all the limits the given number of times, write the results to file
(define (do-sim-steps-to-limit lolimits pops clip? times fname)
  (with-output-to-file 
      (expand-user-path (string->path (string-append "~/research/NormEmerge/" fname)))
    (lambda ()
      (for-each (lambda (augpop)
                  (printf "~a ~a ~a ~a ~a~%" (first augpop) (second augpop) (third augpop) (fourth augpop) (gather-sim-steps-to-limits lolimits (fifth augpop) times clip?)))
                pops))))

;;================================================================================================================
;;================================================================================================================
;; Theory

;; norm-emerge: (listof agent) -> (listof (listof number))
;; Simulate the theory.  For a given population, adapt each agent according to the theoretically computed delta
;; given the current agent bias and the population mean.  Repeat this for as many steps as needed to get close enough to either 0 or 1.
;; ****** NEED TO DOUBLE CHECK THIS ******
(define (norm-emerge pop)
  (local ((define n (length pop))

          ;; cix : agent number -> number
          ;; compute expected probability of given agent coordinating with the others
          (define (cix a pbar)
            (+ (* a (pibar a pbar)) (* (- 1 a) (- 1 (pibar a pbar)))))
          
          ;; new-pbar : (listof agent) number -> (listof agent)
          ;; compute the new biases for the given agents
          (define (new-pbar pop pbar)
            (map (lambda (a)
                   (local ((define d (constant-delta a (pibar a pbar n)))
                           ;;(define d (alt-delta a (pibar a pbar)))
                           ;;(define d (delta-i a (pibar a pbar)))
                           )
                     (if (>= d 0)
                         (min (+ a d) 1.0)
                         (max (+ a d) 0.0))))
                 pop))
          
          ;; alt-delta : agent number -> number
          ;; alternative delta computation that decellerates as converges
          (define (alt-delta a pibar)
            ;;(+ (* (- 1 a) X a pibar)
              ;; (* -1 a X (- 1 a) (- 1 pibar))
              ;; (* a Y (- 1 a) pibar)
              ;; (* -1 (- 1 a) Y a (- 1 pibar)))
            (* (+ X Y)
               (- a (sqr a))
               (- (* 2 pibar) 1))
            )
          
          ;; constant-delta : agent number -> number
          ;; compute expected delta given the
          (define (constant-delta a pibar)
            (+ (* X a pibar) (* -1 X (- 1 a) (- 1 pibar)) (* Y (- 1 a) pibar) (* Y -1 a (- 1 pibar))))

          ;; evolve : (listof agent) (listof number -> (listof number
          ;; run the math for until some criterion, returning the sequence of overall coordination
          (define (evolve pop sequence)
            (cond [(or (> (apply min pop) (- 1 CLOSE-ENOUGH))
                       (< (apply max pop) CLOSE-ENOUGH))
                   (reverse (cons (cons (average pop) pop) sequence))]
                  [else (local ((define new-pop (new-pbar pop (average pop))))
                          (evolve new-pop (cons (cons (average pop) pop) sequence)))]))
          
          )

    (evolve pop
            empty)))

;; Misc norm-emerg (simulated theory) runs -------

(define (non-biased n)
  ;;(norm-emerge (build-list n (lambda (_) (+ 0.5 (- (/ (random 1100) 10000.0) 0.05)))))
  (norm-emerge (make-rand-pop n 0.5 1))
  )

(define (one-sided n M B)
  (norm-emerge (append (build-list (- n M) (lambda (_) (+ 0.5 (- (/ (random 1000) 10000.0) 0.05))))
                       (make-list M B))))

(define (two-sided majority majbias minority minbias)
  (norm-emerge (append (make-list majority majbias)
                       (make-list minority minbias))))

;; Actual (non-simulated) Theory -----------------

;; steps-to-upper-limit: limit population -> number
;; For CONSTANT-UPDATE, compute the theoretically predicted number of steps to reach (or exceed) a limit (where limit > average(pop) > 0.5.
;; In the case of steps-to-one, use limit of 0.5.  That is, limits are in the range -.5 to .5, or 0 to .5 for the upper range.
(define (steps-to-upper-limit limit pop)
  (let* ((pbar (average pop))
         (y (+ (* 2 X) 1))
         (epsilon (- pbar 0.5)))
    (/ (- (log limit) (log epsilon))
       (log y))))

;; gather-ss-by-limit: number population number -> ???
;; calls steps-to-upper-limit, a theory based computation for CONSTANT-UPDATE analysis
(define (gather-ss-by-limit lim pop times)
  (list (+ lim 0.5)
        (* 1.0 (average (build-list times (lambda (_) (sim-steps-to-limit (+ lim 0.5) pop adjust-bias)))))
        (* 1.0 (average (build-list times (lambda (_) (sim-steps-to-limit (+ lim 0.5) pop adj-bias-no-clip)))))
        (steps-to-upper-limit lim pop)))

;; gather...: (listof number) population -> (listof (list number number))
;; get results from theory for steps to each limit for the given population
;; return list of pairs (list limit steps-to-limit).
(define (gather-theory-steps-to-limits lolimits pop)
  (map (lambda (l) (list l (new-steps-to-upper-limit l pop))) lolimits))

;; avg-delta-proportional-A: population -> number
;; using the PROPORTIONAL-A method, 1/n * (P * summate( pbar - a )), compute the average delta based on the given population values.
;; NOTE, this is not (?) making the pi_bar = pbar assumption (10/10/2011: I think it is using pbar instead of pi_bar.)
(define (avg-delta-proportional-A pop)
  (let ((pbar (average pop)))
    (* P (average (map (lambda (a) (- pbar a)) pop)))))

;; avg-delta-proportional-B: population -> number
;; compute the expected or average delta for the population under the PROPORTIANL-B-UPDATE model
;; 2*P*(2*pbar - 1)*(a - a^2)
(define (avg-delta-proportional-B pop)
  (let ((pbar (average pop))
        )
    (* 2 P (- (* 2 pbar) 1.0) (average (map (lambda (a) (- a (sqr a))) pop)))))
(define (avg-delta-proportional-B-variance pop)
  (let ((pbar (average pop))
        (sigsquared (variance pop))
        )
    (* 2 P (- (* 2 pbar) 1.0) (- pbar sigsquared (sqr pbar)))))

;; new-steps-to-lims-over-pops: (listof population) lolims -> (listof (list number number))
;; average the gather-theory-steps-to-limits over a list of populations
(define (new-steps-to-lims-over-pops manypops lims)
  (local ((define mres (map (lambda (p) (gather-theory-steps-to-limits (map (lambda (l) (+ .5 l)) lims) p)) manypops))
          (define (avg-one lopairs)
            (list (caar lopairs) (average (map second lopairs))))
          (define (avg-all lolopairs)
            (cond [(empty? (car lolopairs)) empty]
                  [else (cons (avg-one (map first lolopairs))
                              (avg-all (map cdr lolopairs)))])))
    (avg-all mres)))

;; new-steps-to-upper-limit: number population (agent boolean boolean -> agent) -> number
;; for the given update function and population, compute the number of steps to reach the given limit
(define (new-steps-to-upper-limit limit pop)
  (local ((define (help pop steps)
            (cond [(> (average pop) limit) steps]
                  [else (help (new-pop-proportional-B pop) (add1 steps))])))
    (help pop 0)))

;; new-pop-proportional-B: population -> population
;; compute the new population based on the theory
(define (new-pop-proportional-B pop)
  (let ((pbar (average pop)))
    (map (lambda (a) (+ a (* 2 P (- (* 2 pbar) 1.0) (- a (sqr a)))))
         pop)))


(define my-const-lims '(.5 .49 .48 .47 .45 .43 .4 .35 .3 .25 .2 .15 .1)) ;; '(1 .99 .98 .97 .95 .93 .9 .85 .8 .75 .7 .65 .6)
(define my-proportional-lims '(.45 .43 .4 .35 .3 .25 .2 .15 .1)) ;; took out some for proportional update

;;================================================================================================================
;;================================================================================================================
;; Empirical Simulation

;; new-gather-to-limits: (listof number) population number (agent boolean boolean -> agent) -> (listof (listof (list number number)))
;; given list of limits, a population, and a number of times to repeat the population,
;; and an adjustment, or update, function 
;; simulate the evolution of the population and as it adapts, record the number of steps to each limit as it gets there.
;; Return a list of lists, one for each time, that is a list of limit-step pairs.
;; *** initially, make this work only for the upper range ***
(define (new-gather-to-limits lolimits pop times adjustf)
  (local (;; help: (listof number) population (listof (list number number)) N -> (listof (list number number))
          (define (help lims pop limit-steps steps)
            (let ((mean (average pop)))
              (cond [(empty? lims) limit-steps]
                    [(< mean (- 1.0 (car lims))) "eject"]
                    [(>= mean (car lims)) (help (cdr lims) (one-iteration-f pop adjustf) (cons (list (car lims) steps) limit-steps) (add1 steps))]
                    [else (help lims (one-iteration-f pop adjustf) limit-steps (add1 steps))])))
          )
    (cond [(zero? times) empty]
          [else (let ((one-try (help (sort lolimits <) pop empty 0)))
                  (if (string? one-try)
                      (new-gather-to-limits lolimits pop (sub1 times) adjustf)
                      (cons one-try (new-gather-to-limits lolimits pop (sub1 times) adjustf))))])))

;; map-limit-res: (listof (listof (list number number))) -> (listof (list number number))
;; average the limit results from new-gather-to-limits
(define (map-limit-res res)
  (local ((define n (length res))
          (define (average-one lopairs)
            (list (caar lopairs) (/ (foldr + 0 (map second lopairs)) n)))
          (define (average-all lres)
            (cond [(empty? (car lres)) empty]
                  [else (cons (average-one (map car lres))
                              (average-all (map cdr lres)))])))
    (average-all res)))

;; gather-sim-steps-to-limits : (listof number) population number (agent boolean boolean -> agent) -> (listof number)
;; get results from sim for steps to each limit for given population averaged over times times.
(define (gather-sim-steps-to-limits lolimits pop times adjustf)
  (map (lambda (l) (* 1.0 (average (build-list times (lambda (_) (sim-steps-to-limit (+ l 0.5) pop adjustf))))))
       lolimits))

;; sim-steps-to-limit: number (listof agent) (agent boolean boolean -> agent) -> number
;; run a simulation and count the steps to reach a given limit
;; [updated 2/4/2011 to check the limit on both sides.  *** HOWEVER: need to consider what to do when the population crosses over
;;  the mid-point and converges to the other limit. ***]
(define (sim-steps-to-limit limit pop l-update)
  (local ((define (ssto-aux popi steps)
            (cond [(not (<= (- 1.0 limit) (average popi) limit)) steps]
                  [else
                   (when (zero? (modulo steps 100))
                     (printf "At step ~a: ~a~%" steps (average popi)))
                   (ssto-aux (mk-new-pop popi)
                                  (add1 steps))]))
          (define (mk-new-pop pop)
            (foldr (lambda (p r) (append (interact (first p) (second p) l-update) r)) empty (pair-up pop))
            ;;(one-iteration pop)
            )
          )
    (ssto-aux pop 0)))

;; population is a listof agent

;; make-population : N N number N number -> (listof agent)
;; create a population of total size tot, where f1 have bias b1 and f2 have bias b2, and the remaining (tot - (f1 + f2)) have bias 1/2
(define (make-population tot f1 b1 f2 b2)
  (append (make-list f1 b1)
          (make-list f2 b2)
          (make-list (- tot (+ f1 f2)) 1/2)))

;; Update Methods ---------------------
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
;; a_i = a_i + P*(1-a_i) if chose act and coordinate
;;     = a_i - P*(1-a_i) if chose act and conflict
;;     = a_i + P*a_i if chose NOT-act and conflict
;;     = a_i - P*a_i if chose NOT-act and coordinate
(define (adj-proportional-B a act outcome)
  (cond [act (cond [outcome (- a (* P (- 1.0 a)))]
                   [else (+ a (* P (- 1.0 a)))])]
        [else (cond [outcome (+ a (* P a))]
                    [else (- a (* P a))])]))

;; Agent Interactions --------------------

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


;; pair-up : population -> (listof (list agent agent))
;; create list of random pairs of the agents in the population
(define (pair-up pop)
  (local ((define rand-indices (shuffle-list (build-list (length pop) (lambda (x) x)) 7))
          ;; build-em : (listof N) -> (listof (pair agent agent))
          (define (build-em loi)
            (cond [(empty? loi) empty]
                  [else (cons (list (list-ref pop (first loi))
                                    (list-ref pop (second loi)))
                              (build-em (cdr (cdr loi))))]))
            )
    (build-em rand-indices)
    ))

;; sim-avg-update: population (number boolean boolean -> number) number -> number
;; for the given agent-adjustment function, compute the avg delta for one step of interaction and update over this population
;; optionally, average this over ntimes with the same population but different pairings
(define (sim-avg-update pop adjustf ntimes)
  (average (build-list ntimes
                       (lambda (_) (- (average (one-iteration-f pop adjustf))
                                      (average pop))))))

;; one-iteration-f: population (number boolean boolean -> number) -> population
;; randomly pair-up the given population, have each pair interact and be updated based on the given adjustf function.
(define (one-iteration-f pop adjustf)
  (foldr (lambda (p r) (append (interact (first p) (second p) adjustf) r)) empty (pair-up pop)))

;; one-iteration : population -> (listof agent)
;; pair up the population and perform the paired interactions, creating the new population of agents
;; preserved for backward compatibility
(define (one-iteration pop)
  (one-iteration-f pop adjust-bias))

;; multi-iteration-f : population N (number boolean boolean -> number) -> population
;; do n iterations of one-iteration-f
(define (multi-iteration-f pop n adjustf)
  (cond [(zero? n) pop]
        [else (multi-iteration-f (one-iteration-f pop adjustf) (sub1 n) adjustf)]))

;; multi-iteration : population N -> population
;; do n iterations with pop with default adjustment function adjust-bias
(define (multi-iteration pop n)
  (multi-iteration-f pop n adjust-bias))

;; bias-freq : (listof agent) -> (vectorof N)
;; Count the number of agents with each respective bias.  *** As written, this works only with the constant-update adjustment ***
(define (bias-freq pop)
  (local ((define res (make-vector (add1 (/ 1 CONST-INC)) 0))
          (define (count-em p)
            (cond [(empty? p) res]
                  [else (vector-set! res (/ (car p) CONST-INC) (add1 (vector-ref res (/ (car p) CONST-INC))))
                        (count-em (cdr p))]))
          )
    (count-em pop)))
    
;; collect-multi-iter : (listof agent) N -> (listof (vectorof number))
;; collect the frequency distributions for bias values over the number of iterations
(define (collect-multi-iter pop n)
  (local (;; run-em: N (listof agent) -> (listof (vectorof number))
          ;; run an iteration and store the results
          (define (run-em n pop)
            (cond [(zero? n) empty]
                  [else (cons (bias-freq pop)
                              (run-em (sub1 n) (one-iteration pop)))]))
          )
    (run-em n pop)))

;; merge-two-runs : (listof (vectorof number)) (listof (vectorof number)) -> (listof (vectorof number))
;; take two sequences of frequency distributions and add them
(define (merge-two-runs r1 r2)
  (map (lambda (v1 v2)
         (build-vector (add1 (/ 1 CONST-INC)) (lambda (i) (+ (vector-ref v1 i) (vector-ref v2 i)))))
       r1 r2))

;; do-multi-runs : N N (listof agent) -> (listof (vectorof number))
;; run 'rep' repititions of a simulation out to 'len' iterations with given initial population, 'pop', returning the combined frequency distribution
(define (do-multi-runs rep len pop)
  (foldr merge-two-runs
         (build-list len (lambda (_) (make-vector (add1 (/ 1 CONST-INC)) 0)))
         (build-list rep (lambda (_) (collect-multi-iter pop len)))))


;; dump-frequencies: string (listof (vectorof number)) -> ...
;; for the accumulated data as list of number-vectors, dump the numbers to a file given
(define (dump-frequencies fname list-of-vec)
  (record fname
          (map vector->list list-of-vec)))

;; record : string (listof (listof number)) -> ...
;; for data in a list of number-lists, print each value on a separate line(? -- not sure why).
(define (record fname list-of-lists)
  (with-output-to-file 
    (lambda ()
      (for-each (lambda (timeslice)
                  (for-each (lambda (agent)
                              (printf "~a~%" agent))
                            timeslice)
                  (printf "~%"))
                list-of-lists))))

;; abstract-write-to-file: string string (listof data) -> ...
;; consume the filename arg, the print constrol string, and list of data,
;; write each data item according to the control string to the file given
(define (abstract-write-to-file fname lodata)
  (with-output-to-file
      (expand-user-path (string->path (string-append "~/research/NormEmerge/" fname)))
    (lambda ()
      (for-each (lambda (data-item)
                  (for-each (lambda (datum)
                              (printf "~v~c" datum #\tab))
                            data-item)
                  (printf "~%"))
                lodata))))

;;===================================================================================================
;; Experiment Data Collection

;; create 50 populations of 100 agents for each combination of mean and spread
(define popsplus (create-pops-plus 50 100 (build-list 19 (lambda (n) (* .05 (add1 n)))) (list .25 1 4 8)))

;; getsomeres-by-mean: (listof (list number number population)) (population -> number) -> (listof (list number number))
;; for a list of poptuples consisting mean, spread and population, apply the given function which consumes a population and returns an avg-update,
;; return a list of tuples, one for each mean, consisting of mean followed by average update over all given populations having the same mean.
;; The popsplus argument is presumed to have been filtered by a particular spread.
(define (getsomeres-by-mean popsplus update-f)
  (map (lambda (grp) (list (caar grp) (average (map second grp))))
       (group-by (map (lambda (ptup) (list (first ptup) (update-f (third ptup))))
                      popsplus)
                 0)))

;; pair-up-res: (listof (list number number)) (listof (list number number)) -> (listof (list number number number number)
;; take two results-sets from getsomeres-by-mean, and zip them together as (list mean difference-between-res theory-res simulation-res)
(define (pair-up-res tres sres)
  (map (lambda (t s) (list (first t) (- (second t) (second s)) (second t) (second s))) tres sres))

#|
;; theoretical update for these populations
(define theoryres 
  (getsomeres-by-mean (filter (lambda (ptup) (= (second ptup) 8)) popsplus) avg-delta-proportional-B-variance))

;; simulated update for these same populations, 100 different pairings for each population
(define simres
  (getsomeres-by-mean (filter (lambda (ptup) (= (second ptup) 8)) popsplus) (lambda (p) (sim-avg-update p adj-proportional-B 100))))

;; paird results for .25 spread
(define prdres (pair-up-res theoryres simres))
|#

