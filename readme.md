`Chop` - Round Matrix Elements to Lower Precision
==========

About
-----

`chop' is a MATLAB function for rounding the elements of a matrix to a lower
precision arithmetic with one of several forms of rounding.  Its intended
use is for simulating arithmetic of different precisions (less than double)
with various rounding modes.

The arithmetic formats supported are 
-  'b', 'bfloat16'           - bfloat16,
-  'h', 'half', 'fp16'       - IEEE half precision (the default),
-  's', 'single', 'fp32'     - IEEE single precision,
-  'd', 'double', 'fp64'     - IEEE double precision,
-  'c', 'custom',            - custom format.

Subnormal numbers can be supposrted or not.,
and in the latter case they are flushed to zero.

Several rounding modes are supported:
-     Round to nearest using round to even last bit to break ties
-     (the default).
-     Round towards plus infinity (round up).
-     Round towards minus infinity (round down).
-     Round towards zero.
-     Stochastic rounding - round to the next larger or next smaller
      floating-point number with probability proportional to
      the distance to those floating-point numbers.
-     Stochastic rounding - round to the next larger or next smaller 
      floating-point number with equal probability.

Demonstration function:
- `demo-harmonic' computes the harmonic sum in several arithmetic
   formats using all the supprte rounding modes.

Other M-file:

- `roundit' is a function for rounding a matrix to have integer entries.
  It is used by `chop'. It is not intended to be called irctly.

Test funtions:
- `test_chop' is a test function for `chop'.
- `test_roundit' is a test function for `roundit'.
Each function should print "All tests successful!".

The code was developed in MATLAB R2018b and works with versions at least
back to R2016a.

License
-------

See `license.txt` for licensing information.