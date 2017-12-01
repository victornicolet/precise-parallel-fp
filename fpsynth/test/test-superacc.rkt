#lang racket

(require "../lib.rkt"
         "./testutils.rkt")

(define (test-xsum_small_init)
  (begin
    (define the_xsum (xsum_small_init))
    (and (not (xssum-inf the_xsum))
         (not (xssum-nan the_xsum))
         (= (vector-length
             (vector-filter-not
              (lambda (x) (= x 0))
              (xssum-chunks the_xsum))) 0)
         (= (xssum-adds_untils_propagate the_xsum) XSUM_SMALL_CARRY_TERMS))))

(define (test-carry_propagate)
  (begin
    (define the_xsum (xsum_small_init))
    (xsum_carry_propagate the_xsum)
    #t))

(define (test-xsum_small_add1)
  (define the_xsum (xsum_small_init))
  (define the_vector (vector 1.0 3.0 3.24 3.45 3e-100)) ;Sum should be 7.69
  (for ([i (in-vector the_vector)])
    (set! the_xsum (xsum_small_add1 the_xsum i)))
  (xsum_small_display the_xsum))

(do_all_tests (list test-xsum_small_init
                    test-carry_propagate
                    test-xsum_small_add1))
