#lang rosette

(require rosette/lib/synthax
         "bitwise_synth_ops.rkt")

(define-symbolic n integer?)

(i32ashr n 2)

;; ****************************************

(TwoSum 222e-18 1e-12)

(define a (fpexp2 0.0 0.0))
(fpexp2_add (fpexp2_add a 1e-19) -1e-19)

(define a1 (fpexp8_init 0.0))
(fpexp8_add (fpexp8_add a1 1e-19) 1e18)
