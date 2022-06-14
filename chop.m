function [c,options] = chop(x,options)
%CHOP    Round matrix elements to lower precision.
%   CHOP(X,options) is the matrix obtained by rounding the elements of
%   the real array X to a lower precision arithmetic with one of several
%   forms of rounding.  X should be single precision or double precision
%   and the output will have the same type.  The structure options
%   controls various aspects of the rounding.
%   1. The arithmetic format is specified by options.format, which is one of 
%       'q43', 'fp8-e4m3'       - NVIDIA quarter precision (4 exponent bits,
%                                 3 significand (mantissa) bits),
%       'q52', 'fp8-e5m2'       - NVIDIA quarter precision (5 exponent bits,
%                                 2 significand bits),
%       'b', 'bfloat16'         - bfloat16,
%       'h', 'half', 'fp16'     - IEEE half precision (the default),
%       's', 'single', 'fp32'   - IEEE single precision,
%       'd', 'double', 'fp64'   - IEEE double precision,
%       'c', 'custom'           - custom format.
%      The custom (base 2) format is defined by options.params, which is a
%      2-vector [t,emax] where t is the number of bits in the significand
%      (including the hidden bit) and emax is the maximum value of the
%      exponent.  The minimu exponent is taken to be emin = 1 - emax and
%      the IEEE floating-point number representation is assumed, so that
%      emax and the number of  bits e in the exponent are related by
%      emax = 2^(e-1) - 1.  The values of t and emax are built-in forthe 
%      other formats and will automatically be returned in options.params.
%      options.format = 'd' is intended to be used only with
%      options.subnormal = 0.
%   2. options.subnormal specifies whether subnormal numbers are supported
%      (if they are not, subnormals are flushed to zero):
%        0 = do not support subnormals (the default for bfloat16),
%        1 = support subnormals (the default for the other formats).
%   3. The form of rounding is specified by options.round:
%       1: round to nearest using round to even last bit to break ties
%          (the default);
%       2: round towards plus infinity (round up);
%       3: round towards minus infinity (round down);
%       4: round towards zero;
%       5: stochastic rounding - round to the next larger or next smaller
%          f.p. (floating-point) number with probability proportional to
%          1 minus the distance to those f.p. numbers;
%       6: stochastic rounding - round to the next larger or next smaller 
%          f.p. number with equal probability.
%      For stochastic rounding, exact f.p. numbers are not changed.
%   4. If options.flip = 1 (default 0) then each element of the rounded
%      result has, with probability options.p (default 0.5), a randomly 
%      chosen bit in its significand flipped.  This option is useful for
%      simulating soft errors.
%   5. If options.explim = 0 (default 1) then emax (the maximal
%      exponent) for the specified arithmetic is ignored, so overflow,
%      underflow, or subnormal numbers will be produced only if necessary 
%      for the data type of X.  This option is useful for exploring
%      low precisions independent of range limitations.
%   6. If options.randfunc is supplied, then in stochastic rounding (modes
%      5 and 6) the random numbers used for rounding will be generated
%      using that function. It should be a function that has a single argument
%      for the number of random numbers to generate and returns a vector of
%      the random numbers. By default, the MATLAB rand function is used.
%
%   On the first call: if options is omitted or only partially specified 
%   the defaults stated above are used.
%   On subsequent calls: if options is omitted or empty then the values
%   used in the previous call are re-used; for any missing fields the
%   default is used.
%   The chop options can also be set with CHOP([],options).
%   The options structure is stored internally in a persistent variable
%   and can be obtained with [~,options] = CHOP.

% References:
% [1] IEEE Standard for Floating-Point Arithmetic, IEEE Std 754-2019
% (revision of IEEE Std 754-2008), IEEE Computer Society, 2019.
% https://ieeexplore.ieee.org/document/8766229/
% [2] Intel Corporation, BFLOAT16---hardware numerics definition,  Nov. 2018, 
% White paper. Document number 338302-001US.
% https://software.intel.com/en-us/download/bfloat16-hardware-numerics-definition
% [3] M. Croci, M. Fasi, N. J. Higham, T. Mary, and M. Mikaitis.
% Stochastic rounding: Implementation, error analysis and applications.
% Roy. Soc. Open Sci., 9(3):1-25, 2022.

if nargin >= 1 && ~isreal(x), error('Chop requires a real input array.'), end
persistent fpopts
reset_format_settings = 0;

if isempty(fpopts) && (nargin <= 1 || (nargin == 2 && isempty(options)))
      fpopts.format = 'h'; fpopts.subnormal = 1;
      fpopts.round = 1; fpopts.flip = 0; fpopts.p = 0.5;
      fpopts.explim = 1;
      fpopts.randfunc = @(n) rand(n, 1);
      reset_format_settings = 1;
elseif nargin == 2 && ~isempty(options)
    % This is not the first call, but fpopts might have all empty fields.
    if ~isfield(options,'format') || ...
        isfield(options,'format') && isempty(options.format)
       options.format = 'h';
    end
    fpopts.format = options.format;
    if isfield(options,'subnormal') && ~isempty(options.subnormal)
         fpopts.subnormal = options.subnormal;
    else  
         if ismember(fpopts.format, {'b','bfloat16'})
             fpopts.subnormal = 0;
         else
             fpopts.subnormal = 1; 
         end
    end    
    if isfield(options,'round') && ~isempty(options.round)
       fpopts.round = options.round;
    else   
       fpopts.round = 1;
    end    
    if isfield(options,'flip') && ~isempty(options.flip)
       fpopts.flip = options.flip;
    else
       fpopts.flip = 0;
    end    
    if isfield(options,'p') && ~isempty(options.p)
       fpopts.p = options.p;
    else    
       fpopts.p = 0.5;
    end    
    if isfield(options,'explim') && ~isempty(options.explim)
       fpopts.explim = options.explim;
    else
       fpopts.explim = 1;
    end
    if isfield(options,'randfunc') && ~isempty(options.randfunc)
       fpopts.randfunc = options.randfunc;
    else
       fpopts.randfunc = @(n) rand(n, 1);
    end
    reset_format_settings = 1;
end

persistent t
persistent emax

if reset_format_settings
    if ismember(fpopts.format, {'h','half','fp16','b','bfloat16','s', ...
                                'single','fp32','d','double','fp64',...
                                'q43','fp8-e4m3','q52','fp8-e5m2'})
      if ismember(fpopts.format, {'q43','fp8-e4m3'})
           % Significand: 3 bits plus 1 hidden. Exponent: 4 bits.
           t = 4; emax = 7;
      elseif ismember(fpopts.format, {'q52','fp8-e5m2'})
           % Significand: 2 bits plus 1 hidden. Exponent: 5 bits.
           t = 3; emax = 15;
      elseif ismember(fpopts.format, {'h','half','fp16'})
           % Significand: 10 bits plus 1 hidden. Exponent: 5 bits.
           t = 11; emax = 15;
       elseif ismember(fpopts.format, {'b','bfloat16'})
           % Significand: 7 bits plus 1 hidden. Exponent: 8 bits.
           t = 8; emax = 127;  
       elseif ismember(fpopts.format, {'s','single','fp32'})
           % Significand: 23 bits plus 1 hidden. Exponent: 8 bits.
           t = 24; emax = 127;
       elseif ismember(fpopts.format, {'d','double','fp64'})
           % Significand: 52 bits plus 1 hidden. Exponent: 11 bits.
           t = 53; emax = 1023;
       end
       fpopts.params = [t emax];
    elseif ismember(fpopts.format, {'c','custom'})
       if nargin == 2 && ~isempty(options)
          if isfield(options,'params') && ~isempty(options.params)
             fpopts.params(1) = options.params(1);
             fpopts.params(2) = options.params(2);
             % Need "p_2 \ge 2p_1 + 2" to avoid double rounding problems
             % in round-to-nearest.
             if fpopts.round == 1
                maxfraction = isa(x,'single') * 11 + isa(x,'double') * 25;
             else
                maxfraction = isa(x,'single') * 23 + isa(x,'double') * 52;
             end
             if (fpopts.params(1) > maxfraction)
                error(['Precision of the custom format must be at most ' ...
                       '%d if working in %s.'], maxfraction, class(x));
             end
          end
       elseif ~isfield(fpopts,'params') || isempty(fpopts.params)
          error('Must specify options.params with options.format = ''c''.')
       end
       t = fpopts.params(1); emax = fpopts.params(2);
    else
       error('Unrecognized format.')  
    end
end

if nargout == 2, options = fpopts; end
if nargin == 0 || isempty(x), if nargout >= 1, c = []; end, return, end
 
emin = 1-emax;            % Exponent of smallest normalized number.
xmin = 2^emin;            % Smallest positive normalized number.
emins = emin + 1 - t;     % Exponent of smallest positive subnormal number.
xmins = 2^emins;          % Smallest positive subnormal number.
xmax = 2^emax * (2-2^(1-t));

% Use the representation:
% x = 2^e * d_1.d_2...d_{t-1} * s, s = 1 or -1.

c = x;
[~,e] = log2(abs(x));
e = e - 1;
ktemp = (e < emin & e >= emins);
if fpopts.explim
   k_sub = find(ktemp); k_norm = find(~ktemp);
else 
  k_sub = []; k_norm = 1:length(ktemp(:)); % Do not limit exponent.
end   

c(k_norm) = pow2(roundit(pow2(x(k_norm), t-1-e(k_norm)), fpopts), ...
                 e(k_norm)-(t-1));
if ~isempty(k_sub)
   t1 = t - max(emin-e(k_sub),0);
   c(k_sub) = pow2(roundit( pow2(x(k_sub), t1-1-e(k_sub)), fpopts ), ...
                   e(k_sub)-(t1-1));
end

if fpopts.explim
  switch(fpopts.round)
    case {1,6}
   % Any number larger than xboundary rounds to inf [1, p. 16].
   xboundary = 2^emax * (2-(1/2)*2^(1-t));
   c(find(x >= xboundary)) = inf;   % Overflow to +inf.
   c(find(x <= -xboundary)) = -inf; % Overflow to -inf.
    case 2
      c(find(x > xmax)) = inf;
      c(find(x < -xmax & x ~= -inf)) = -xmax;
    case 3
      c(find(x > xmax & x ~= inf)) = xmax;
      c(find(x < -xmax)) = -inf;
    case {4,5}
      c(find(x > xmax & x ~= inf)) = xmax;
      c(find(x < -xmax & x ~= -inf)) = -xmax;
  end
  % Round to smallest representable number or flush to zero.
  if fpopts.subnormal == 0
    min_rep = xmin;
  else
    min_rep = xmins;
  end
  k_small = abs(c) < min_rep;
  switch(fpopts.round)
    case 1
      if fpopts.subnormal == 0
        k_round = k_small & abs(c) >= min_rep/2;
      else
        k_round = k_small & abs(c) > min_rep/2;
      end
      c(k_round) = sign(c(k_round)) * min_rep;
      c(k_small & ~k_round) = 0;
    case 2
      k_round = k_small & c > 0 & c < min_rep;
      c(k_round) = min_rep;
      c(k_small & ~k_round) = 0;
    case 3
      k_round = k_small & c < 0 & c > -min_rep;
      c(k_round) = -min_rep;
      c(k_small & ~k_round) = 0;
    case {4,5,6}
      c(k_small) = 0;
  end
end
end