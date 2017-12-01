#lang racket

(require rosette "../lib.rkt" "testutils.rkt")

(define bv64 (bitvector 64))

(define (test_iashr)
  (and
   (= (iashr 1 2 bv64) 0)
   (= (iashr 2 2 bv64) 0)
   (= (iashr 5 1 bv64) 2)
   (= (iashr 16 2 bv64) 4)
   (= (iashr 1 -2 bv64) 0)))

(define (test_ilshr)
  (and
   (= (ilshr 1 2 bv64) 0)
   (= (ilshr 2 2 bv64) 0)
   (= (ilshr 4 2 bv64) 1)
   (= (ilshr 16 2 bv64) 4)
   (= (ilshr (expt 2 64) 2 bv64) 0)
   (= (ilshr (expt 2 62) 2 bv64) (expt 2 60))
   (= (ilshr 1 -2 bv64) 0)))


(define (test_ishl)
  (and
   (= (ishl (expt 2 62) 2 bv64) 0)
   (= (ishl (+ (expt 2 60) 3) 2 bv64) (* 4 (+ (expt 2 60) 3)))
   (= (ishl (expt 2 62) 2 bv64) 0)))

(define (test_i32ashr)
  (and
   (= (i32ashr 1 2) 0)
   (= (i32ashr 2 2) 0)
   (= (i32ashr 5 1) 2)
   (= (i32ashr 30 1) 15)
   (= (i32ashr 1 -2) 0)))

(define (test_i32lshr)
  (and
   (= (i32lshr 1 2) 0)
   (= (i32lshr 2 2) 0)
   (= (i32lshr 4 2) 1)
   (= (i32lshr 16 2) 4)
   (= (i32lshr (expt 2 32) 2) 0)
   (= (i32lshr (expt 2 12) 2) (expt 2 10))
   (= (i32lshr 1 -2) 0)))


(define (test_i32shl)
  (and
   (= (i32shl (expt 2 32) 2) 0)
   (= (i32shl (+ (expt 2 28) 3) 2) (* 4 (+ (expt 2 28) 3)))
   (= (i32shl (expt 2 32) 2) 0)))


(define (test_i64ashr)
  (and
   (= (i64ashr 1 2) 0)
   (= (i64ashr 2 2) 0)
   (= (i64ashr 5 1) 2)
   (= (i64ashr 16 2) 4)
   (= (i64ashr 1 -2) 0)))

(define (test_i64lshr)
  (and
   (= (i64lshr 1 2) 0)
   (= (i64lshr 2 2) 0)
   (= (i64lshr 4 2) 1)
   (= (i64lshr 16 2) 4)
   (= (i64lshr (expt 2 64) 2) 0)
   (= (i64lshr (expt 2 12) 2) (expt 2 10))
   (= (i64lshr 1 -2) 0)))


(define (test_i64shl)
  (and
   (= (i64shl (expt 2 62) 2) 0)
   (= (i64shl (+ (expt 2 60) 3) 2) (* 4 (+ (expt 2 60) 3)))
   (= (i64shl (expt 2 62) 2) 0)))

(define (test_TwoSum)
  (= 0.0 (cadr (TwoSum 1.0 1.0)))
  (not (= 0.0 (cadr (TwoSum 1.0 1e-32))))
  (= 1.0 (car (TwoSum 1.0 1e-32))))

(define (test_TwoMul)
  (= 1.0 (car (TwoMul 1.0 1.0)))
  (= 1e-32 (car (TwoMul 1e-32 1.0))))


;; Test teh extended floating point represenetations
(define (test_fpexp2_addv)
  (begin
    (define inputs (vector 1.0 1.0 1.0 -1e32 -1.0 1e32 -1.0 -1.0 1e20 1e-34))
    (define expansion (fpexp2_init 0.0))
    (fpexp2=? (fpexp2 1e20 1e-34) (fpexp2_addv expansion inputs))))

(define (test_fpexp4_addv)
  (begin
    (define inputs (vector 1.0 1.0 1.0 -1e32 -1.0 1e32 -1.0 -1.0 ))
    (define expansion (fpexp4_init 0.0))
    (fpexp4=? (fpexp4 -2.0 2.0 0.0 0.0) (fpexp4_addv expansion inputs))))

(define (test_fpexp8_addv)
  (begin
    (define inputs (vector 1.0 1.0 1.0 -1e32 -1.0 1e32 -1.0 -1.0 ))
    (define expansion (fpexp8_init 0.0))
    (fpexp8=? (fpexp8 -2.0 2.0 0.0 0.0 0.0 0.0 0.0 0.0)
       (fpexp8_addv expansion inputs))))


(define (test_IEE754_doubles)
  (and
    (= (IEE754_double_mantissa 1.0) 0)
    (= (IEE754_double_exp 1.0) 1023)
    (= (IEE754_sign 1.0) 0)
    (= (IEE754_sign -1.0) 1)
    (= (IEE754_double->real (real->IEE754_double 1.0)) 1.0)
    (= (IEE754_double->real (real->IEE754_double 3.1415)) 3.1415)
    (= (IEE754_double->real (real->IEE754_double 2e-15)) 2e-15)))



(printf "~n**** Testing bitwise_synth_ops ***~n")

(for ([test (list test_iashr
                  test_ilshr
                  test_ishl
                  test_i32ashr
                  test_i32lshr
                  test_i32shl
                  test_i64ashr
                  test_i64lshr
                  test_i64shl
                  test_TwoSum
                  test_TwoMul
                  test_fpexp2_addv
                  test_fpexp4_addv
                  test_fpexp8_addv
                  test_IEE754_doubles
                  )])
  (if (test)
      (printf "~a ok~n" test)
      (printf "Failed ~a~n" test)))
