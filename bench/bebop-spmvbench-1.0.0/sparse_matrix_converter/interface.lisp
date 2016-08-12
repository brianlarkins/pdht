(require 'asdf)
(pushnew "/Users/mhoemmen/pkg/cffi/" asdf:*central-registry* :test 'equal)
(asdf:oos 'asdf:load-op :cffi)

(defpackage #:sp
  (:use #:common-lisp #:cffi #:cffi-utils)
  (:export load! save! format? convert! mult destroy!))

(in-package #:sp)

;;; The Sparse Matrix Converter library
(define-foreign-library lib-sparse-matrix-converter
  (t (:default "libsparse_matrix_converter")))

;;; The BeBOP Utility Library
(define-foreign-library lib-bebop-util
  (t (:default "libbebop_util")))

(use-foreign-library lib-bebop-util)
(use-foreign-library lib-sparse-matrix-converter)


(defcfun ("sp_load" load!) :pointer
  (path :string)
  (fmt :string))

(defcfun ("sp_save" save!) :int
  (A :pointer)
  (path :string)
  (fmt :string))

(defcfun ("sp_format" format?) :void 
  (A :pointer))

(defcfun ("sp_convert" convert!) :int
  (A :pointer)
  (type :string))

(defcfun ("sp_mult" mult) :pointer
  (B :pointer)
  (A :pointer))

(defcfun ("sp_destroy" %destroy!) :void
  (A :pointer))

(defun destroy! (A)
  (progn
    (assert (and (not (null A)) (not (cffi:null-pointer-p A))))
    (%destroy! A)))
