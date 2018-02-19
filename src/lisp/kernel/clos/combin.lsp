;;;;  -*- Mode: Lisp; Syntax: Common-Lisp; Package: CLOS -*-
;;;;
;;;;  Copyright (c) 1992, Giuseppe Attardi.
;;;;  Copyright (c) 2001, Juan Jose Garcia Ripoll.
;;;;
;;;;    ECoLisp is free software; you can redistribute it and/or
;;;;    modify it under the terms of the GNU Library General Public
;;;;    License as published by the Free Software Foundation; either
;;;;    version 2 of the License, or (at your option) any later version.
;;;;
;;;;    See file '../Copyright' for full details.


(in-package "CLOS")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; COMPILING EFFECTIVE METHODS
;;;
;;; The following functions take care of transforming the forms
;;; produced by the method combinations into effective methods. In ECL
;;; effective methods are nothing but directly callable functions.
;;; Ideally, this compilation should just produce new compiled
;;; functions. However, we do not want to cons a lot of functions, and
;;; therefore we use closures.
;;;
;;; Formerly we used to keep a list of precompiled effective methods
;;; and made a structural comparison between the current method and
;;; the precompiled ones, so as to save memory. This only causes
;;; improvements in declarative combinations. For standard combinations
;;; it should be enough with a couple of different closures and hence
;;; the structural comparison is a loss of time.
;;;
;;; This is the core routine. It produces effective methods (i.e.
;;; functions) out of the forms generated by the method combinators.
;;; We consider the following cases:
;;;  1) Ordinary methods. The function of the method is extracted.
;;;  2) Functions. They map to themselves. This only happens
;;;     when these functions have been generated by previous calls
;;;     to EFFECTIVE-METHOD-FUNCTION.
;;;  3) (CALL-METHOD method rest-methods) A closure is
;;;	generated that invokes the current method while informing
;;;	it about the remaining methods.
;;;  4) (MAKE-METHOD form) A function is created that takes the
;;;	list of arguments of the generic function and evaluates
;;;	the forms in a null environment. This is the only form
;;;	that may lead to consing of new bytecodes objects. Nested
;;;	CALL-METHOD are handled via the global macro CALL-METHOD.
;;;  5) Ordinary forms are turned into lambda forms, much like
;;;	what happens with the content of MAKE-METHOD.
;;;
(defvar *avoid-compiling* nil)
(defun emf-maybe-compile (form)
  (if *avoid-compiling*
      (coerce form 'function)
      (let ((*avoid-compiling* t)
            ;; Cleavir itself uses generic functions, and we could therefore
            ;; end up here recursively, which ends quite badly.
            ;; So fall back to the bclasp compiler which does not use gfs.
            (cmp:*cleavir-compile-hook* nil))
        (compile nil form))))
;; The CLHS description of call-method and make-method is kind of hard to understand.
;; But it does say rather specifically that make-method's form is evaluated in the null
;; lexical environment plus a definition for call-method, and nothing else from CL.
;; So no call-next-method.
(defun emf-make-method (form)
  (emf-maybe-compile
   `(lambda (.method-args. .next-methods.)
      (declare (ignore .next-methods.))
      ,form)))
;; process a "method" argument to call-method
(defun emf-call-method (form)
  (cond ((method-p form) (method-function form))
        ((and (consp form) (eq (first form) 'make-method) (consp (cdr form)))
         (emf-make-method (second form)))
        (t "Invalid argument to CALL-METHOD: ~a" form)))

(defun effective-method-function (form)
  (if (consp form)
      (case (first form)
        ((make-method)
         (emf-make-method (second form)))
        ((call-method)
         (combine-method-functions
          (emf-call-method (second form))
          (mapcar #'emf-call-method (third form))))
        (otherwise (emf-make-method form)))
      (emf-make-method form)))

;;;
;;; This function is a combinator of effective methods. It creates a
;;; closure that invokes the first method while passing the information
;;; of the remaining methods. The resulting closure (or effective method)
;;; is the equivalent of (CALL-METHOD method rest-methods)
(defun combine-method-functions (method rest-methods)
  (lambda (.method-args. .next-methods.)
    (declare (ignorable .next-methods. #|no-next-methods|#)
             (core:lambda-name combine-method-functions.lambda))
    ;; TODO: Optimize this application and GF dispatch should be more efficient
    ;; .method-args. can be a valist or a regular list
    (funcall method .method-args. rest-methods)))

(defmacro call-method (method &optional rest-methods)
  `(funcall ,(emf-call-method method)
            ;; This macro is only invoked from effective-method-function, above,
            ;; in the case of having to fall back to the compiler.
            ;; Therefore .method-args. will always be bound in the environment.
            .method-args.
            ',(and rest-methods (mapcar #'emf-call-method rest-methods))))

(defun error-qualifier (m qualifier)
  (error "Standard method combination allows only one qualifier ~
          per method, either :BEFORE, :AFTER, or :AROUND; while ~
          a method with ~S was found."
	 m qualifier))

;; ----------------------------------------------------------------------
;; DEFINE-METHOD-COMBINATION
;;
;; METHOD-COMBINATION objects are instances defined in hierarchy.lsp.
;; They have slots for the name, compiler, and options. Name is obvious,
;; and the options are those provided to the thing.
;; The "compiler" is somewhat misleadingly named; it's the function that
;; outptus the effective method form.
;; These functions are stored in the global *method-combinations* hash
;; table. (the standard method on) FIND-METHOD-COMBINATION ignores the gf,
;; and makes a new METHOD-COMBINATION instance with the "compiler" looked
;; up in the hash table, and the name and options.
;; The "compiler" functions take two arguments, plus the lambda-list from
;; the define-method-combination. The first argument is the generic function
;; (used for the :generic-function option of D-M-C), the second is the sorted
;; list of applicable methods, and the rest are the method combination options.
;;

#+threads
(defparameter *method-combinations-lock* (mp:make-lock :name 'find-method-combination))
(defparameter *method-combinations* (make-hash-table :size 32 :test 'eq))


(defun search-method-combination (name)
  (mp:with-lock (*method-combinations-lock*)
    (or (gethash name *method-combinations*)
	(error "~A does not name a method combination" name))))

(defun install-method-combination (name function)
  (mp:with-lock (*method-combinations-lock*)
                (setf (gethash name *method-combinations*) function))
  name)

(defun make-method-combination (name compiler options)
  (with-early-make-instance +method-combination-slots+
    (o (find-class 'method-combination)
       :name name
       :compiler compiler
       :options options)
    o))

;; Will be upgraded into a generic function later.
(defun find-method-combination (gf method-combination-type-name method-combination-options)
  (make-method-combination method-combination-type-name
			   (search-method-combination method-combination-type-name)
			   method-combination-options
			   ))

(defun define-simple-method-combination (name &key documentation
					 identity-with-one-argument
					 (operator name))
  `(define-method-combination
     ,name (&optional (order :MOST-SPECIFIC-FIRST))
     ((around (:AROUND))
      (principal (,name) :REQUIRED t))
     ,documentation
     (let ((main-effective-method
             `(,',operator ,@(mapcar #'(lambda (x)
                                         (declare (core:lambda-name define-simple-method-combination.lambda))
                                         `(CALL-METHOD ,x NIL))
				    (if (eql order :MOST-SPECIFIC-LAST)
					(reverse principal)
					principal)))))
       (cond (around
	      `(call-method ,(first around)
		(,@(rest around) (make-method ,main-effective-method))))
	     (,(if identity-with-one-argument
		   '(rest principal)
		   t)
	      main-effective-method)
	     (t (second main-effective-method))))))

(defun define-complex-method-combination (form)
  (flet ((syntax-error ()
	   (error "~S is not a valid DEFINE-METHOD-COMBINATION form"
		  form)))
    (destructuring-bind (name lambda-list method-groups &rest body &aux
			      (group-names '())
			      (group-checks '())
			      (group-after '())
			      (generic-function '.generic-function.)
			      (method-arguments '()))
	form
      (unless (symbolp name) (syntax-error))
      (let ((x (first body)))
	(when (and (consp x) (eql (first x) :ARGUMENTS))
	  (error "Option :ARGUMENTS is not supported in DEFINE-METHOD-COMBINATION.")))
      (let ((x (first body)))
	(when (and (consp x) (eql (first x) :GENERIC-FUNCTION))
	  (setf body (rest body))
	  (unless (symbolp (setf generic-function (second x)))
	    (syntax-error))))
      (dolist (group method-groups)
	(destructuring-bind (group-name predicate &key description
					(order :most-specific-first) (required nil))
	    group
	  (if (symbolp group-name)
	      (push group-name group-names)
	      (syntax-error))
	  (let ((condition
                  (cond ((eql predicate '*) 'T)
                        ((null predicate) `(null .method-qualifiers.))
                        ((symbolp predicate)
                         `(,predicate .METHOD-QUALIFIERS.))
                        ((consp predicate)
                         (let* ((q (last predicate 0))
                                (p (copy-list (butlast predicate 0))))
                           (when (every #'symbolp p)
                             (if (eql q '*)
                                 `(every #'equal ',p .METHOD-QUALIFIERS.)
                                 `(equal ',p .METHOD-QUALIFIERS.)))))
                        (t (syntax-error)))))
	    (push `(,condition (push .METHOD. ,group-name)) group-checks))
	  (when required
	    (push `(unless ,group-name
                     ;; Effective methods can be computed in other situations than being about to call them.
                     ;; As such, compute-effective-method should not signal an error unless the computation
                     ;; is impossible. Lacking a required method is by contrast a problem that only needs to
                     ;; be signaled when the function is actually being called. So we return an error form.
                     ;; (NO-REQUIRED-METHOD is defined in fixup, but we could move it here.)
                     ;; FIXME: This form is not actually evaluated, in favor of weird bullshit.
                     ;; See compute-effective-function-maybe-optimize.
                     (return-from ,name
                       '(no-required-method ,generic-function ',group-name .arguments.)))
		  group-after))
	  (case order
	    (:most-specific-first
	     (push `(setf ,group-name (nreverse ,group-name)) group-after))
	    (:most-specific-last)
	    (otherwise
             (let ((order-var (gensym)))
               (setf group-names (append group-names (list (list order-var order)))
                     group-after (list* `(when (eq ,order-var :most-specific-first)
                                           (setf ,group-name (nreverse ,group-name)))
                                        group-after)))))))
      `(install-method-combination ',name
				   (lambda (,generic-function .methods-list. ,@lambda-list)
                                     (declare (core:lambda-name ,name))
                                     (block ,name 
                                       (let (,@group-names)
                                         (dolist (.method. .methods-list.)
                                           (let ((.method-qualifiers. (method-qualifiers .method.)))
                                             (cond ,@(nreverse group-checks)
                                                   (t (invalid-method-error .method.
                                                                            "Method qualifiers ~S are not allowed in the method~
			      combination ~S." .method-qualifiers. ',name)))))
                                         ,@group-after
                                         ,@body)))))))

(defmacro define-method-combination (name &body body)
  (if (and body (listp (first body)))
      (define-complex-method-combination (list* name body))
      (apply #'define-simple-method-combination name body)))

(defun method-combination-error (format-control &rest args)
  ;; FIXME! We should emit a more detailed error!
  (error "Method-combination error:~%~S"
	 (apply #'format nil format-control args)))

(defun invalid-method-error (method format-control &rest args)
  (error "Invalid method error for ~A~%~S"
	 method
	 (apply #'format nil format-control args)))



;;; ----------------------------------------------------------------------
;;; COMPUTE-EFFECTIVE-METHOD
;;;

(defun compute-effective-method-function (gf method-combination applicable-methods)
  (effective-method-function
   (compute-effective-method gf method-combination applicable-methods)))

;; will be upgraded into being the standard method on compute-effective-method in fixup.
(defun std-compute-effective-method (gf method-combination applicable-methods)
  (declare (type method-combination method-combination)
	   (type generic-function gf)
	   (optimize speed (safety 0)))
  ;; FIXME: early accessors here could technically be bad, if someone subclasses method-combination
  ;; On the other hand, I've never seen anyone do that. D-M-C already has arbitrary code, and
  ;; method combinations have no defined accessors - all you could do is add methods to
  ;; compute-effective-method, itself unusual because, again, arbitrary code.
  (with-early-accessors (+method-combination-slots+)
    (let* ((compiler (method-combination-compiler method-combination))
	   (options (method-combination-options method-combination)))
      (if options
	  (apply compiler gf applicable-methods options)
	  (funcall compiler gf applicable-methods)))))

(define-method-combination standard ()
    ((around (:around))
     (before (:before))
     (primary () :required t)
     (after (:after)))
  (flet ((call-methods (methods)
           (mapcar (lambda (method)
                     `(call-method ,method))
                   methods)))
    ;; We're a bit more hopeful about avoiding make-method and m-v-p1 than
    ;; the example in CLHS define-method-combination.
    ;; Performance impact is likely to be marginal at best, but why not try?
    (let* ((call-primary `(call-method ,(first primary) ,(rest primary)))
           (call-before (if before
                            `(progn ,@(call-methods before) ,call-primary)
                            call-primary))
           (call-after (if after
                           `(multiple-value-prog1 ,call-before
                              ,@(call-methods (reverse after)))
                           call-before))
           (call-around (if around
                            (if (and (null before) (null after))
                                `(call-method ,(first around)
                                              (,@(rest around)
                                               ,@primary))
                                `(call-method ,(first around)
                                              (,@(rest around)
                                               (make-method ,call-after))))
                            call-after)))
      call-around)))

(define-method-combination progn :identity-with-one-argument t)
(define-method-combination and :identity-with-one-argument t)
(define-method-combination max :identity-with-one-argument t)
(define-method-combination + :identity-with-one-argument t)
(define-method-combination nconc :identity-with-one-argument t)
(define-method-combination append :identity-with-one-argument nil)
(define-method-combination list :identity-with-one-argument nil)
(define-method-combination min :identity-with-one-argument t)
(define-method-combination or :identity-with-one-argument t)
