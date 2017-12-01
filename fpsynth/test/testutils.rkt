#lang racket

(define (_assert k)
  (if k (raise 0) (raise -1)))

(define (do_test e n)
  (with-handlers ([number? (when (< n 0) (printf "Test ~a failed." n))]) (e)))

(define (do_all_tests tests)
(for ([test tests])
  (if (test)
      (printf "~a ok~n" test)
      (printf "Failed ~a~n" test))))

(provide _assert do_test do_all_tests)
