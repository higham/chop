`Chop` - MATLAB code for rounding matrix elements to lower precision
==========

About
-----

`chop` is a MATLAB function for rounding the elements of a matrix to a
lower precision arithmetic with one of several forms of rounding.  Its
intended use is for simulating arithmetic of different precisions (less
than double) with various rounding modes. The input to `chop` should be
single precision or double precision and the output will have the same
type: the lower precision numbers are stored within a higher precision type.

The arithmetic formats supported are 
-  'b', 'bfloat16'           - bfloat16,
-  'h', 'half', 'fp16'       - IEEE half precision (the default),
-  's', 'single', 'fp32'     - IEEE single precision,
-  'd', 'double', 'fp64'     - IEEE double precision,
-  'c', 'custom'            - custom format.

Subnormal numbers can be supported or not,
and in the latter case they are flushed to zero.

Several rounding modes are supported:
- Round to nearest using round to even last bit to break ties
  (the default).
- Round towards plus infinity (round up).
- Round towards minus infinity (round down).
- Round towards zero.
- Stochastic rounding - round to the next larger or next smaller
  floating-point number with probability proportional to
  the distance to those floating-point numbers.
- Stochastic rounding - round to the next larger or next smaller 
  floating-point number with equal probability.

A further option causes each element of the rounded result 
to have, with a specified probability defaulting to 0.5,
a randomly chosen bit in its significand flipped. 

Demonstration function:
- `demo-harmonic` computes the harmonic series in several arithmetic
   formats using all the supported rounding modes.

Other M-file:

- `roundit` is a function for rounding a matrix to have integer entries.
  It is used by `chop` and is not intended to be called directly.

Test functions:
- `test_chop` is a test function for `chop`.
- `test_roundit` is a test function for `roundit`.

Each test function should print "All tests successful!".

The function `chop` is a successor to a function of the same name in the
[The Matrix Computation Toolbox](http://www.ma.man.ac.uk/~higham/mctoolbox/)
(also available on
[File Exchange](https://uk.mathworks.com/matlabcentral/fileexchange/2360-the-matrix-computation-toolbox)).

The test function `test_chop` needs the function 
`float_params`from the repository
[float_params](https://github.com/higham/float_params).
That function is included in this repository for convenience, but may not
be the latest version.

Requirements
---------

The code was developed in MATLAB R2018b and works with versions at least
back to R2016a.

Reference
---------

Nicholas J. Higham and Srikara Pranesh, [Simulating Low Precision
Floating-Point Arithmetic](http://eprints.maths.manchester.ac.uk/2692), MIMS Eprint 2019.4, Manchester Institute for Mathematical
Sciences, The University of Manchester, UK, March 2019.

License
-------

See `license.txt` for licensing information.