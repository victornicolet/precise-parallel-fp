#lang rosette

(require "./utils.rkt"
         racket/generic)
;; Rosette doesn't give us the bitwise operations

(define bv32? (bitvector 32))
(define bv64? (bitvector 64))


(define (iashr i n t?)
  (bitvector->integer (bvashr (integer->bitvector i t?)
                              (integer->bitvector n t?))))

(define (ishl i n t?)
  (bitvector->integer (bvshl (integer->bitvector i t?)
                              (integer->bitvector n t?))))

(define (ilshr i n t?)
  (bitvector->integer (bvlshr (integer->bitvector i t?)
                              (integer->bitvector n t?))))

(define (i32ashr i n)
  (bitvector->integer (bvashr (integer->bitvector i bv32?)
                              (integer->bitvector n bv32?))))

(define (i32shl i n)
  (bitvector->integer (bvshl (integer->bitvector i bv32?)
                              (integer->bitvector n bv32?))))

(define (i32lshr i n)
  (bitvector->integer (bvlshr (integer->bitvector i bv32?)
                              (integer->bitvector n bv32?))))
(define (i32and i n)
  (bitvector->integer (bvand (integer->bitvector i bv32?)
                              (integer->bitvector n bv32?))))

(define (i32not i) (bitvector->integer (bvnot (integer->bitvector i bv32?))))

(define (i64ashr i n)
  (bitvector->integer (bvashr (integer->bitvector i bv64?)
                              (integer->bitvector n bv64?))))
(define (i64shl i n)
  (bitvector->integer (bvshl (integer->bitvector i bv64?)
                              (integer->bitvector n bv64?))))
(define (i64lshr i n)
  (bitvector->integer (bvlshr (integer->bitvector i bv64?)
                              (integer->bitvector n bv64?))))
(define (i64and i n)
  (bitvector->integer (bvand (integer->bitvector i bv64?)
                              (integer->bitvector n bv64?))))
(define (i64not i) (bitvector->integer (bvnot (integer->bitvector i bv64?))))

(define (i64sub i n)
  (bitvector->integer (bvsub (integer->bitvector i bv64?)
                              (integer->bitvector n bv64?))))

(define (xadd sor des) sor)

;; Summing two floating-point numbers introduces an error, but this error is
;; always representable as another floating point number.
(define (TwoSum a b)
  (let* ([r (+ a b)]
        [z (- r a)])
    (list r (+ (- a (- r z)) (- b z)))))


;; Mutliplying two floating point numbers is less error-prone than addition, the
;; exponents are added and the significands are mutliplied. There is no loss of
;; 'correct numbers'
(define (TwoMul a b)
  (list (* a b) 0.0))


;; Floating point expansions

;; Floating-point expansion of size 2

(struct fpexp2 (s1 s2) #:transparent)

(define (fpexp2->real expansion)
  (+ (fpexp2-s1 expansion)
     (fpexp2-s2 expansion)))

(define (fpexp2=? fp1 fp2)
  (and (= (fpexp2-s1 fp1) (fpexp2-s1 fp2))
       (= (fpexp2-s2 fp1) (fpexp2-s2 fp2))))


(define (fpexp2_add fpexp x)
  (let* ([re0 (TwoSum (fpexp2-s1 fpexp) x)]
         [re1 (TwoSum (fpexp2-s2 fpexp) (cadr re0))])
    (fpexp2 (car re0) (car re1))))


(define (fpexp2_addv fpexp v)
  (foldv fpexp2_add fpexp v))


(define (fpexp2_init x) (fpexp2 x 0.0))


(struct fpexp4 (s1 s2 s3 s4) #:transparent)

(define (fpexp4->real fp1)
  (+ (fpexp4-s1 fp1) (fpexp4-s2 fp1) (fpexp4-s3 fp1) (fpexp4-s4 fp1)))

(define (fpexp4=? fp1 fp2)
  (and (= (fpexp4-s1 fp1) (fpexp4-s1 fp2))
       (= (fpexp4-s2 fp1) (fpexp4-s2 fp2))
       (= (fpexp4-s3 fp1) (fpexp4-s3 fp2))
       (= (fpexp4-s4 fp1) (fpexp4-s4 fp2))))

(define (fpexp4_add fpexp x)
  (let* ([re0 (TwoSum (fpexp4-s1 fpexp) x)]
         [re1 (TwoSum (fpexp4-s2 fpexp) (cadr re0))]
         [re2 (TwoSum (fpexp4-s3 fpexp) (cadr re1))]
         [re3 (TwoSum (fpexp4-s4 fpexp) (cadr re2))])
    (fpexp4 (car re0) (car re1) (car re2) (car re3))))

(define (fpexp4_addv fpexp v)
  (foldv fpexp4_add fpexp v))

(define (fpexp4_init x) (fpexp4 x 0.0 0.0 0.0))


(struct fpexp8 (s1 s2 s3 s4 s5 s6 s7 s8) #:transparent)

(define (fpexp8->real fp1)
  (+ (fpexp8-s1 fp1) (fpexp8-s2 fp1) (fpexp8-s3 fp1) (fpexp8-s4 fp1)
     (fpexp8-s5 fp1) (fpexp8-s6 fp1) (fpexp8-s7 fp1) (fpexp8-s8 fp1)))


(define (fpexp8=? fp1 fp2)
  (and (= (fpexp8-s1 fp1) (fpexp8-s1 fp2))
       (= (fpexp8-s2 fp1) (fpexp8-s2 fp2))
       (= (fpexp8-s3 fp1) (fpexp8-s3 fp2))
       (= (fpexp8-s4 fp1) (fpexp8-s4 fp2))
       (= (fpexp8-s5 fp1) (fpexp8-s5 fp2))
       (= (fpexp8-s6 fp1) (fpexp8-s6 fp2))
       (= (fpexp8-s7 fp1) (fpexp8-s7 fp2))
       (= (fpexp8-s8 fp1) (fpexp8-s8 fp2))))


(define (fpexp8_add fpexp x)
  (let* ([re0 (TwoSum (fpexp8-s1 fpexp) x)]
         [re1 (TwoSum (fpexp8-s2 fpexp) (cadr re0))]
         [re2 (TwoSum (fpexp8-s3 fpexp) (cadr re1))]
         [re3 (TwoSum (fpexp8-s4 fpexp) (cadr re2))]
         [re4 (TwoSum (fpexp8-s5 fpexp) (cadr re3))]
         [re5 (TwoSum (fpexp8-s6 fpexp) (cadr re4))]
         [re6 (TwoSum (fpexp8-s7 fpexp) (cadr re5))]
         [re7 (TwoSum (fpexp8-s8 fpexp) (cadr re6))])
    (fpexp8 (car re0) (car re1) (car re2) (car re3)
            (car re4) (car re5) (car re6) (car re7))))

(define (fpexp8_addv fpexp v)
  (foldv fpexp8_add fpexp v))

(define (fpexp8_init x) (fpexp8 x 0.0 0.0 0.0
                                0.0 0.0 0.0 0.0))


;; Some additional bitvector tools
(define (bitvector-bit bvec i)
  (extract i i bvec))

(define (bitvector->string bvec size)
  (define res "")
  (for ([i (in-range 0 size)])
    (let ([bvi (bitvector->integer (concat (bv 0 1) (bitvector-bit bvec i)))])
    (set! res (string-append res (number->string bvi)))))
  res)

;; Get fp mantissas / exponent
(define X64_MANTISSA_BITS 52)
(define X64_EXP_BITS 11)
(define X64_MANTISSA_MASK (i64sub (i64shl 1 X64_MANTISSA_BITS) 1))
(define X64_EXP_MASK  (i64shl (i64sub (i64shl 1 X64_EXP_BITS) 1) X64_MANTISSA_BITS))
(define (bytes->uint64 x) (integer-bytes->integer x #f #f))
(define (uint64->bytes x) (integer->integer-bytes x 8 #f))
(define (bytes->int64 x) (integer-bytes->integer x #t #f))
(define (int64->bytes x) (integer->integer-bytes x 8 #t #f))
(define (bytes->double x) (floating-point-bytes->real x #f))
(define (double->bytes x) (real->floating-point-bytes x 8 #f))
(define (double->int64 x) (bytes->int64 (double->bytes x)))
(define (int64->double x) (bytes->double (int64->bytes x)))

(define (IEE754_double_mantissa_bv x)
  (extract X64_MANTISSA_BITS 0
           (bvand
            (bv (bytes->int64 (double->bytes x)) bv64?)
            (bv X64_MANTISSA_MASK bv64?))))

(define (IEE754_double_mantissa x)
  (bitvector->integer (IEE754_double_mantissa_bv x)))

(define (IEE754_add_implicit_1 mantissa)
  (bitvector->integer
   (bvor (bv mantissa 64) (bvshl (bv 1 64) (bv X64_MANTISSA_BITS 64)))))

(define (IEE754_double_exp_bv x)
   (extract X64_EXP_BITS 0
             (bvlshr (bvand (bv (bytes->int64 (double->bytes x)) bv64?)
                            (bv X64_EXP_MASK bv64?))
                     (bv X64_MANTISSA_BITS bv64?))))


(define (IEE754_double_exp x) (bitvector->integer (IEE754_double_exp_bv x)))

(define (IEE754_sign x)
  (bitvector->integer
   (concat (bv 0 1)
           (bitvector-bit (bv (bytes->int64 (double->bytes x)) bv64?) 63))))

(struct IEE754_double (sign mantissa exponent) #:transparent)

(define (real->IEE754_double x)
  (IEE754_double (IEE754_sign x)
                 (IEE754_double_mantissa x)
                 (IEE754_double_exp x)))

(define (IEE754_double->real x)
  (let
      ([sign (IEE754_double-sign x)]
        [mantissa (IEE754_double-mantissa x)]
        [exponent (IEE754_double-exponent x)])
    (let
        ([unsigned_real (* (+ (* mantissa (expt 2 -52)) 1.0)
                              (expt 2 (- exponent 1023)))])
      (if (= sign 0) unsigned_real (- unsigned_real)))))


(provide iashr ilshr ishl
         i32ashr i32lshr i32shl i32and i32not
         i64ashr i64lshr i64shl i64and i64not
         double->int64 int64->double
         IEE754_double_exp
         IEE754_double_exp_bv
         IEE754_double_mantissa IEE754_double_mantissa_bv IEE754_add_implicit_1
         IEE754_sign
         IEE754_double->real
         IEE754_double
         IEE754_double-sign IEE754_double-mantissa IEE754_double-exponent
         real->IEE754_double
         TwoSum TwoMul
         fpexp2 fpexp2=? fpexp2_add fpexp2_init fpexp2_addv
         fpexp4 fpexp4=? fpexp4_add fpexp4_init fpexp4_addv
         fpexp8 fpexp8=? fpexp8_add fpexp8_init fpexp8_addv)
