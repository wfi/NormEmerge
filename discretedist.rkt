#lang racket

(require test-engine/racket-tests)

;; This was starter code for exploring the relationships between the number of agents
;; and the number of possible values that the population mean could take on (under a constant update regime).
;; Apparently intended to shed light on clipping issues in the constant update regime.
;; see notebook, 10/1/2011, pg 5

;; one-agent: number -> (listof number)
;; for m possible values between 0 and 1 (inclusive), build the list of possible values this agent can take on
(check-expect (one-agent 3) '(0 1/2 1))
(define (one-agent m)
  (build-list m (lambda (i) (/ i (sub1 m)))))

;; n-agents: number number -> (listof (listof number))
;; for n agents, create the list of possible values that each agent may take on if each could have one of m values
(check-expect (n-agents 2 3) '((0 1/2 1)(0 1/2 1)))
(define (n-agents n m)
  (build-list n (lambda (i) (one-agent m))))

;; mix-all: number number -> (listof (listof number))
;; for the result from n-agents, pick one value from each list in all possible combinations
(define (mix-all n m)
  (let ((base (n-agents n m))
        (iassigns (list (list))))
    (foldr (lambda (i sf) (mix-one i sf))
           iassigns
           base)))

;; mix-one: (listof number) (listof (listof number)) -> (listof (listof number))
;; for each possible value in vals, cons that value onto each collection of assignments in assigns
(define (mix-one vals assigns)
  (foldr (lambda (v sf) (append (map (lambda (p) (cons v p)) assigns)
                                sf))
         empty
         vals))

;; add-mean: (listof (listof number)) -> (listof (list number (listof number)))
;; add the mean of set of assignments
(define (add-mean loa)
  (let ((n (length (car loa))))
    (sort (map (lambda (a) (list (/ (foldr + 0 a) n) a)) loa)
          <
          #:key car)))

;; consolidate: (listof (list number (listof number))) -> (listof (cons number (listof number)))
;; consolidate agent assignments having the same mean
(define (consolidate l)
  (foldl (lambda (ds sf) 
           (if (= (car ds) (car (car sf)))
               (cons (cons (car ds) (cons (second ds) (cdr (car sf))))
                     (cdr sf))
               (cons ds sf)))
         (list (car l))
         (cdr l)))

;; count-same: (listof (list number (listof number))) -> (listof number)
;; count the number of assignments in output from add-mean, having the same mean by calling consolidate and ignoring the mean
(define (count-same l)
  (map (lambda (s) (sub1 (length s)))
       (consolidate l)))


;; prob-one: (listof (cons number (listof number))) -> (listof (listof number))
;; process output from consolidate to determine the fraction of times an agent is one for each mean
(define (prob-one l)
  (map (lambda (d) (list (car d)
                         (length (filter (lambda (x) (= 1 (car x))) (cdr d)))
                         (length (cdr d))))
       l))

;; fact: number -> number
(define (fact n)
  (cond [(zero? n) 1]
        [else (* n (fact (sub1 n)))]))

;; choose: number number -> number
;; n choose k
(define (choose n k)
  (/ (fact n) (* (fact k) (fact (- n k)))))

(build-list 10 (lambda (i) (choose 10 i)))



(test)
