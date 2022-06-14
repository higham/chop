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
- 'q43', 'fp8-e4m3'         - NVIDIA quarter precision (4 exponent bits,
                              3 significand (mantissa) bits),
- 'q52', 'fp8-e5m2'         - NVIDIA quarter precision (5 exponent bits,
                              2 significand bits),
-  'b', 'bfloat16'          - bfloat16,
-  'h', 'half', 'fp16'      - IEEE half precision (the default),
-  's', 'single', 'fp32'    - IEEE single precision,
-  'd', 'double', 'fp64'    - IEEE double precision,
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

Optionally, each element of the rounded result has, with a specified
probability defaulting to 0.5, a randomly chosen bit in its significand
flipped.  This option is useful for simulating soft errors

A further option causes the exponent limit for the specified arithmetic to
be ignored, so overflow, underflow, or subnormal numbers will be produced
only if necessary for the data type of the input.  This option is useful
for exploring low precisions indepdent of range limitations.

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

Usage
-----

There are two main usages of `chop`.
First, one can pass `options` with every call:

```
options.format = 's'; options.round = 5; options.subnormal = 1; 
...
A(i,j) = chop(A(i,j) - chop(A(i,k) * A(k,j),options),options);
```

Here, `options.format = 's'` specifies that the precision is single,
`options.round = 5` specifies stochastic rounding, mode 1
and `options.subnormal = 1` specifies that subnormal numbers are
not flushed to zero. 
For full details of the options see the help lines in `chop.m`.

The above usage is rather tedious and produces cluttered code.
Instead we can set up the arithmetic parameters on a call of the form 
`chop([],options)` and exploit the fact that subsequent calls 
with just one input argument will reuse the previously specified `options`:

```
options.format = 's'; options.round = 5; options.subnormal = 1; 
chop([],options)
...
A(i,j) = chop(A(i,j) - chop(A(i,k)*A(k,j))); 
```

The current value of `options` is stored inside the function
(in a persistent variable, whose value is retained during the session until
the function is cleared with `clear chop` and can be obtained with 

```
[~,options] = chop
```

Requirements
---------

The code was developed in MATLAB R2018b to R2020a and works with versions
at least back to R2016a.

Reference
---------

Nicholas J. Higham and Srikara Pranesh, [Simulating Low Precision
Floating-Point Arithmetic](https://epubs.siam.org/doi/10.1137/19M1251308), 
SIAM J. Sci. Comput., 41(4):A2536-A2551, 2019.


Related 
---------

A C library [CPFloat](https://github.com/mfasi/cpfloat)
provides similar functionality to `chop`.
It includes a MEX interface, and the MEX code can be faster than `chop`.

Acknowledgements
---------

The code was written by Nick Higham and Sri Pranesh.
Max Fasi and Mantas Mikaitis have contributed improvements.
Ian McInerney contributed the quarter precision formats
and the options.randfunc code.

License
-------

See `license.txt` for licensing information.
