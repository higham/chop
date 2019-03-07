function [c,options] = chop(x,options)
%CHOP    Round matrix elements to lower precision.
%   CHOP(X,options) is the matrix obtained by rounding the elements of
%   the array X to a lower precision arithmetic with one of several
%   forms of rounding.  The arithmetic format is specified by
%   options.format, which is one of 
%     'b', 'bfloat16'           - bfloat16,
%     'h', 'half', 'fp16'       - IEEE half precision (the default),
%     's', 'single', 'fp32'     - IEEE single precision,
%     'd', 'double', 'fp64'     - IEEE double precision,
%     'c', 'custom',            - custom format.
%   In the last case the format is defined by 
%   options.params, which is a 2-vector [t, emax] where t is the
%   number of bits in the significand (including the hidden bit) and
%   emax is the maximum value of the exponent.  The values of t and emax
%   are built-in for b, h, s, d and will automatically be returned in
%   options.params. 
%   options.subnormal specifies whether subnormal numbers are supported
%   (if they are not, subnormals are flushed to zero):
%      0 = do not support subnormals (the default for bfloat16),
%      1 = support subnormals (the default for fp16, fp32 and fp64).
%   The form of rounding is specified by options.round:
%     1: round to nearest using round to even last bit to break ties
%        (the default);
%     2: round towards plus infinity (round up);
%     3: round towards minus infinity (round down);
%     4: round towards zero;
%     5: stochastic rounding - round to the next larger or next smaller
%        f.p. (floating-point) number with probability proportional to
%        the distance to those f.p. numbers;
%     6: stochastic rounding - round to the next larger or next smaller 
%        f.p. number with equal probability.
%   For stochastic rounding, exact f.p. numbers are not changed.
%   If options.flip = 1 (default 0) then each element of the rounded result 
%   has, with probability options.p (default 0.5), a randomly chosen bit
%   in its significand flipped. 
%   On the first call: if options is omitted or only partially specified 
%   the defaults stated above are used.
%   On subseqeuent calls: if options is omitted or empty then the values used 
%   in the previous call are re-used.  For any missing fields the
%   default is used.
%   The options structure is stored internally in a persistent variable
%   and can be obtained with [~,options] = CHOP.

% References:
% [1] IEEE Standard for Floating-Point Arithmetic, IEEE Std 754-2008 (revision 
% of IEEE Std 754-1985), 58, IEEE Computer Society, 2008; pages 8,
% 13. https://ieeexplore.ieee.org/document/461093
% [2] Intel Corporation, BFLOAT16---hardware numerics definition,  Nov. 2018, 
% White paper. Document number 338302-001US.
% https://software.intel.com/en-us/download/bfloat16-hardware-numerics-definition

persistent fpopts

if isempty(fpopts) && (nargin <= 1 || (nargin == 2 && isempty(options)))
      fpopts.format = 'h'; fpopts.subnormal = 1;
      fpopts.round = 1; fpopts.flip = 0; fpopts.p = 0.5;
elseif nargin == 2 && ~isempty(options)
    % This is not the first call, but fpopts might have all empty fields.
    if isfield(options,'format') && isempty(options.format)
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
end

if ismember(fpopts.format, {'h','half','fp16','b','bfloat16','s', ...
                            'single','fp32','d','double','fp64'})
   if ismember(fpopts.format, {'h','half','fp16'})
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
      end    
   elseif ~isfield(fpopts,'params') || isempty(fpopts.params)
      error('Must specify options.params with options.format = ''c''.')
   end    
   t = fpopts.params(1); emax = fpopts.params(2);
end   

if nargout == 2, options = fpopts; end
if nargin == 0 || isempty(x), if nargout >= 1, c = []; end, return, end

if fpopts.flip == 1, fpopts.t = t; end
    
emin = 1-emax;            % Exponent of smallest normalized number.
xmin = 2^emin;            % Smallest positive normalized number.
emins = emin + 1 - t;     % Exponent of smallest positive subnormal number.
xmins = 2^emins;          % Smallest positive subnormal number.
xmax = 2^emax * (2-2^(1-t));

% Use the representation:
% x = 2^e * d_1.d_2...d_{t-1} * s, s = 1 or -1.

c = x;
e = floor(log2(abs(x)));
ktemp = (e < emin & e >= emins);
k_sub = find(ktemp);
k_norm = find(~ktemp);

c(k_norm) = pow2(roundit(pow2(x(k_norm), t-1-e(k_norm)), fpopts), ...
                 e(k_norm)-(t-1));
if ~isempty(k_sub)
   t1 = t - max(emin-e(k_sub),0);
   c(k_sub) = pow2(roundit( pow2(x(k_sub), t1-1-e(k_sub)), fpopts ), ...
                   e(k_sub)-(t1-1));
   if fpopts.subnormal == 0
     c(k_sub) = 0;   % Flush subnormals to zero.
   end  
end  

% Any number large than xboundary rounds to inf [1, p. 16].
xboundary = 2^emax * (2-(1/2)*2^(1-t));
c(find(x >= xboundary)) = inf;   % Overflow to +inf.
c(find(x <= -xboundary)) = -inf; % Overflow to -inf.
c(find(abs(x) < xmins)) = 0;     % Underflow to zero.
