#lang racket


(define (foldv proc initial vect)
  (let ([x (box initial)])
    (begin
      (for ([v vect])
        (set-box! x (proc (unbox x) v)))
      (unbox x))))

(define (vector-first-rev v x)
  (begin
    (define f -1)
    (define break #f)
    (for ([i (in-range (sub1 (vector-length v)) 0 -1)])
      #:break break
      (when (x (vector-ref v i))
        (set! f i)
        (set! break #t)))
    f))

(define (vector-next v i x)
  (begin
    (define f -1)
    (define break #f)
    (for ([i (in-range i (sub1 (vector-length v)) 1)])
      #:break break
      (when (x (vector-ref v i))
        (set! f i)
        (set! break #t)))
    f))


(define (vector-last-nonzero v)
  (vector-first-rev v (lambda (x) (not (= x 0)))))

(define (vector-next-nonzero v i)
  (vector-next v i (lambda (x) (not (= x 0)))))


(provide foldv vector-first-rev vector-last-nonzero vector-next-nonzero)
