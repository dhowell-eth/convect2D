GFORTRAN module version '10' created from poisson_solver_module.f95
MD5:fc5a90a7d8c05be9d537982fd8aae1b6 -- If you edit this, you'll get what you deserve.

(() () () () () () () () () () () () () () () () () () () () () () () ()
() () ())

()

()

()

()

()

(2 'iteration_2dpoisson' 'poissonsolver' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION ALWAYS_EXPLICIT) (REAL 4 0 0 0
REAL ()) 3 0 (4 5 6 7) () 8 () () () 0 0)
9 'linear_interpolate' 'poissonsolver' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE IMPLICIT_PURE ALWAYS_EXPLICIT) (
UNKNOWN 0 0 0 0 UNKNOWN ()) 10 0 (11) () 0 () () () 0 0)
12 'poissonsolver' 'poissonsolver' '' 1 ((MODULE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) (UNKNOWN 0 0 0 0 UNKNOWN ()) 0 0 () ()
0 () () () 0 0)
13 'prolongate' 'poissonsolver' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
0 UNKNOWN ()) 14 0 (15 16) () 0 () () () 0 0)
17 'residue_2dpoisson' 'poissonsolver' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE IMPLICIT_PURE ALWAYS_EXPLICIT) (
UNKNOWN 0 0 0 0 UNKNOWN ()) 18 0 (19 20 21 22) () 0 () () () 0 0)
23 'restrict' 'poissonsolver' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE IMPLICIT_PURE ALWAYS_EXPLICIT) (
UNKNOWN 0 0 0 0 UNKNOWN ()) 24 0 (25 26) () 0 () () () 0 0)
27 'rms_2d' 'poissonsolver' '' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 FUNCTION IMPLICIT_PURE ALWAYS_EXPLICIT) (REAL 4 0 0 0
REAL ()) 28 0 (29) () 30 () () () 0 0)
31 'vcycle_2dpoisson' 'poissonsolver' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION RECURSIVE ALWAYS_EXPLICIT) (REAL 4
0 0 0 REAL ()) 32 0 (33 34 35) () 36 () () () 0 0)
4 'u' '' '' 3 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 4 0 0 0 REAL ()) 0 0 () (2 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0
0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
5 'f' '' '' 3 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
DUMMY) (REAL 4 0 0 0 REAL ()) 0 0 () (2 0 ASSUMED_SHAPE (CONSTANT (
INTEGER 4 0 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0 0 INTEGER
()) 0 '1') ()) 0 () () () 0 0)
6 'h' '' '' 3 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (
REAL 4 0 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
7 'alpha' '' '' 3 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(REAL 4 0 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
8 'rms_residue' '' '' 3 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 RESULT ALWAYS_EXPLICIT) (REAL 4 0 0 0 REAL ()) 0 0 () () 0 ()
() () 0 0)
11 'grid_cell' '' '' 10 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (REAL 4 0 0 0 REAL ()) 0 0 () (2 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0
0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
15 'coarse' '' '' 14 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 4 0 0 0 REAL ()) 0 0 () (2 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0
0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
16 'fine' '' '' 14 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 4 0 0 0 REAL ()) 0 0 () (2 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0
0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
19 'u' '' '' 18 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
DUMMY) (REAL 4 0 0 0 REAL ()) 0 0 () (2 0 ASSUMED_SHAPE (CONSTANT (
INTEGER 4 0 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0 0 INTEGER
()) 0 '1') ()) 0 () () () 0 0)
20 'f' '' '' 18 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
DUMMY) (REAL 4 0 0 0 REAL ()) 0 0 () (2 0 ASSUMED_SHAPE (CONSTANT (
INTEGER 4 0 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0 0 INTEGER
()) 0 '1') ()) 0 () () () 0 0)
21 'h' '' '' 18 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (
REAL 4 0 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
22 'res' '' '' 18 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 4 0 0 0 REAL ()) 0 0 () (2 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0
0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
25 'fine' '' '' 24 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 4 0 0 0 REAL ()) 0 0 () (2 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0
0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
26 'coarse' '' '' 24 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 4 0 0 0 REAL ()) 0 0 () (2 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0
0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
29 'a' '' '' 28 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
DUMMY) (REAL 4 0 0 0 REAL ()) 0 0 () (2 0 ASSUMED_SHAPE (CONSTANT (
INTEGER 4 0 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0 0 INTEGER
()) 0 '1') ()) 0 () () () 0 0)
30 'rms' '' '' 28 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 RESULT ALWAYS_EXPLICIT) (REAL 4 0 0 0 REAL ()) 0 0 () () 0 () () ()
0 0)
33 'u_f' '' '' 32 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 4 0 0 0 REAL ()) 0 0 () (2 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0
0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
34 'rhs' '' '' 32 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 4 0 0 0 REAL ()) 0 0 () (2 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0
0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
35 'h' '' '' 32 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY) (
REAL 4 0 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
36 'resv' '' '' 32 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 RESULT ALWAYS_EXPLICIT) (REAL 4 0 0 0 REAL ()) 0 0 () () 0 ()
() () 0 0)
)

('iteration_2dpoisson' 0 2 'linear_interpolate' 0 9 'poissonsolver' 0 12
'prolongate' 0 13 'residue_2dpoisson' 0 17 'restrict' 0 23 'rms_2d' 0 27
'vcycle_2dpoisson' 0 31)
