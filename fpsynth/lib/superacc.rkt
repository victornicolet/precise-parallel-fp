#lang racket

(require "bitwise_synth_ops.rkt" "utils.rkt")

(define statuses
  (list "Exact" "Inexact" "MinusInfinity" "PlusInfinity" "Overflow" "sNaN" "qNaN"))

(define debug #t)


(define XSUM_MANTISSA_BITS 52)
(define XSUM_EXP_BITS 11)
(define XSUM_MANTISSA_MASK (- (i64shl 1 XSUM_MANTISSA_BITS) 1))
(define XSUM_EXP_MASK  (- (i64shl 1 XSUM_EXP_BITS) 1))
(define XSUM_EXP_BIAS (- (i64shl 1 (- XSUM_EXP_BITS 1)) 1))
(define XSUM_SIGN_BIT (+ XSUM_MANTISSA_BITS XSUM_EXP_BITS))
(define XSUM_SIGN_MASK  (i64shl 1 XSUM_SIGN_BIT))

;; CONSTANTS DEFINING THE SMALL ACCUMULATOR FORMAT.
(define XSUM_SCHUNK_BITS 64)
(define XSUM_LOW_EXP_BITS 5)
(define XSUM_LOW_EXP_MASK (- (i64shl 1 XSUM_LOW_EXP_BITS) 1))
(define XSUM_HIGH_EXP_BITS (- XSUM_EXP_BITS  XSUM_LOW_EXP_BITS))
(define XSUM_HIGH_EXP_MASK (- (i64shl 1  XSUM_HIGH_EXP_BITS) 1))
(define XSUM_SCHUNKS (+ (i64shl 1 XSUM_HIGH_EXP_BITS) 3))
(define XSUM_LOW_MANTISSA_BITS  (i64shl 1 XSUM_LOW_EXP_BITS))
(define XSUM_HIGH_MANTISSA_BITS (- XSUM_MANTISSA_BITS XSUM_LOW_MANTISSA_BITS))
(define XSUM_LOW_MANTISSA_MASK (- (i64shl 1 XSUM_LOW_MANTISSA_BITS) 1))
(define XSUM_SMALL_CARRY_BITS (- (- XSUM_SCHUNK_BITS 1) XSUM_MANTISSA_BITS))
(define XSUM_SMALL_CARRY_TERMS  (- (i64shl 1 XSUM_SMALL_CARRY_BITS) 1))

(struct xssum (chunks inf nan adds_untils_propagate) #:transparent)

(define XSUM_LCHUNK_BITS 64)
(define XSUM_LCOUNT_BITS (- 64 XSUM_MANTISSA_BITS))
(define XSUM_LCHUNKS (i64shl 1 (+ XSUM_EXP_BITS 1)))

(struct xlsum (chunks counts chunks_used uses_used sacc))

;; Function for small accumulators
(define (xsum_small_init)
  (xssum (make-vector XSUM_SCHUNKS 0) 0 0 XSUM_SMALL_CARRY_TERMS))

(define (xsum_carry_propagate sacc)
  (begin
    (define chunks (xssum-chunks sacc))
    (define uix 0)
    (define i 0)
    ;; Set u to the uppermost non-zero chunk
    (define u (vector-last-nonzero chunks))
    (when (>= u 0)
      (begin
        (define break #f)
        (with-handlers
          ([integer? (set! break #f)])
          (do ([em 0
                   (begin
                     (set! i (vector-next-nonzero chunks i))
                     ;; If non non-zero chunk have been found, exit the loop
                     (when (= i -1) (raise 1))

                     ;; Propagate teh carry from this chunks to next non-zero
                     ;; chunk up.
                     (let* ([c (vector-ref i chunks)]
                            [chigh (i64ashr c XSUM_LOW_MANTISSA_BITS)]
                            [clow (i64and c XSUM_LOW_MANTISSA_MASK)])
                       (begin
                         (if
                          (= chigh 0)
                          ;; chigh = 0
                          (begin
                            (set! uix i)
                            (set! i (add1 i))) ;; end then
                          ;; chigh != 0
                          (begin
                            (when (= u i)
                              (if (= chigh -1)
                                  (begin (set! uix i)
                                         (raise 1)) ; Exit the loop
                                  (set! u (add1 i))))
                            (when (not (= clow 0))
                              (set! u (add1 i)))
                            ;; Change chunk[i]
                            (vector-set!
                             chunks i clow)
                            ;; Add to chunk[i+1] the high part of chunk[i]
                            (vector-set!
                             chunks (add1 i)
                             (+ (vector-ref chunks (add1 i)) chigh))
                            ;;  increment i
                            (set! i (add1 i))) ;; end else
                          ); end if chigh = 0
                         ); end begin
                       )
                     0)])
              (
               (> i u)                  ; Stop when past last non-zero chunk
               ) ; end do loop final
            ); end do loop
          ) ; end with handlers loop
        (if (< uix 0)
            (set! uix 0)
            (begin
              (do ([cchk (vector-ref chunks uix)
                         (begin
                           (vector-set! chunks (sub1 uix)
                                        (+ (vector-ref chunks (sub1 uix))
                                           (i64shl -1 XSUM_LOW_MANTISSA_BITS)))
                           (vector-set! chunks uix 0)
                           (set! uix (sub1 uix))
                           )])
                  (
                   (not (and (= cchk -1) (> uix 0))) ;; stop?
                   )
                ); end do loop
              );end begin
            ); end if uix < 0
        ); end begin
      ); end when u >= 0
    (list
     (xssum chunks
            (xssum-inf sacc)
            (xssum-nan sacc)
            (sub1 XSUM_SMALL_CARRY_TERMS))
     uix)))

(define (xsum_small_add_inf_nan sacc ivalue)
  (define mantissa (i64and ivalue XSUM_MANTISSA_MASK))
  (define new_inf (xssum-inf sacc))
  (define new_nan (xssum-nan sacc))
  (if (= mantissa 0)
      (set! new_inf
            (cond [(= (xssum-inf sacc) 0) ivalue]
                  [(not (= ivalue (xssum-inf sacc))) ivalue]
                  [else new_inf])) ; end then
      (set! new_nan
            (cond
              [(<= (i64and new_nan XSUM_MANTISSA_MASK) mantissa)
               (i64and ivalue (i64not XSUM_SIGN_MASK))]
              [else new_nan])) ; end else
      ) ;end if mantissa = 0
  (xssum
   (xssum-chunks sacc)
   new_inf
   new_nan
   (xssum-adds_untils_propagate sacc)) ; the structure returned
  )


(define (xsum_add1_no_carry sacc flt)
  (when debug (begin (printf "Adding ~a~n" flt)))
  (define chunks (xssum-chunks sacc))
  (define ivalue (double->int64 flt))
  (define flt_repr (real->IEE754_double flt))
  (define mantissa (IEE754_double-mantissa flt_repr))
  (define exponent (IEE754_double-exponent flt_repr))
  (define return #f)
  (cond
    [(and (not (= exponent 0)) (not (= exponent XSUM_EXP_MASK)))
     (begin
       (when debug (printf "Mantissa: ~a, exponent: ~a~n" mantissa exponent))
       (set! mantissa (IEE754_add_implicit_1 mantissa))
       )]
    [(= exponent 0)
     (if (= mantissa 0)
         (begin
           (when debug (printf "Mantissa and exponent are 0~n"))
           (set! return #t)
           )
         (set! exponent 1))]
    [else (begin
            (printf "Is NAN or Inf.~n")
            (set! sacc (xsum_small_add_inf_nan sacc ivalue))
            (set! return #t))])
  (when (not return)
    (let* ([low_exp (i64and exponent XSUM_LOW_EXP_MASK)]
           [high_exp (i64lshr exponent XSUM_LOW_EXP_BITS)]
           [chunk_ptr high_exp]
           [chunk0 (vector-ref chunks high_exp)]
           [chunk1 (vector-ref chunks (add1 high_exp))]
           [low_m (i64and (i64shl mantissa low_exp) XSUM_LOW_MANTISSA_MASK)]
           [high_m (i64lshr mantissa (- XSUM_LOW_MANTISSA_BITS low_exp))])
      (begin
        (when debug
          (begin(printf "[lo:~a hi:~a]~n" low_exp high_exp))
          (printf "[vlo:~a vhi:~a]~n" low_m high_m)
          (printf "[c0:~a c1:~a]~n" chunk0 chunk1)) ; end optional debug info
        (if (< ivalue 0)
            (begin
              (vector-set! chunks chunk_ptr (- chunk0 low_m))
              (vector-set! chunks (add1 chunk_ptr) (- chunk0 high_m))
              ); end then of if ivalue < 0
            (begin
              (vector-set! chunks chunk_ptr (+ chunk0 low_m))
              (vector-set! chunks (add1 chunk_ptr) (+ chunk0 high_m))
              ); end else of if ivalue < 0
            ); fi ivalue < 0
        )
      ); end let
    );end when not return
  ; return point
  (xssum chunks
         (xssum-inf sacc)
         (xssum-nan sacc)
         (xssum-adds_untils_propagate sacc))
  )

(define (xsum_small_add1 sacc flt)
  (let ([carry_propagated
         (if (= (xssum-adds_untils_propagate sacc) 0)
             (car (xsum_carry_propagate sacc))
             sacc)])
    (let ([value_added (xsum_add1_no_carry sacc flt)])
      (xssum (xssum-chunks value_added)
             (xssum-inf value_added)
             (xssum-nan value_added)
             (sub1 (xssum-adds_untils_propagate value_added))))))

(define (xsum_small_addv sacc flt len) sacc)

(define (xsum_small_round_real sacc)
  (define ivalue 0)
  (define current_acc sacc)
  (define lower_chunk 0)
  (define-values (i j e more) 0 0 0 0)
  (let ([propagated (xsum_carry_propagate sacc)])
    (set! i (cadr propagated))
    (set! current_acc (car propagated))
    (set! ivalue (vector-ref (xssum-chunk current_acc) i))
    (if (<= i 1)
        (printf ""); then i<= 1
        (printf ""); else i > 1
        ) ; end if
    ); end let after propagation
)

(define (xsum_small_round sacc)
  (when debug (printf "Rounding small accumualtor.~n"))
  (cond
    [(not (= (xssum-nan sacc) 0)) (int64->double (xssum-nan sacc))]
    [(not (= (xssum-inf sacc) 0)) (int64->double (xssum-inf sacc))]
    [else (xsum_small_round_real sacc)]))

(define (xsum_small_display sacc)
  (printf " _____________________________________________________~n")
  (printf " |   Inf:~a |  Nan:~a   |   Adds until propagate: ~a |~n\
  ----------------------------------------------------~n~a~n"
          (xssum-inf sacc)
          (xssum-nan sacc)
          (xssum-adds_untils_propagate sacc)
          (xssum-chunks sacc)))


(define (xsum_small_chunks_used sacc) sacc)

(provide xssum xssum-inf xssum-adds_untils_propagate xssum-nan xssum-chunks
         xsum_carry_propagate
         xsum_small_init xsum_small_addv xsum_small_add1 xsum_small_round
         xsum_small_display xsum_small_chunks_used
         XSUM_SMALL_CARRY_TERMS)
