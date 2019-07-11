function [u,xmins,xmin,xmax,p,emins,emin,emax] = float_params(prec)
%FLOAT_PARAMS   Parameters for floating-point arithmetic.
%   [u,xmins,xmin,xmax,p,emin,emax] = FLOAT_PARAMS(prec) returns 
%     u:     the unit roundoff,
%     xmins: the smallest positive (subnormal) floating-point number,
%     xmin:  the smallest positive normalized floating-point number,
%     xmax:  the largest floating-point number,
%     p:     the number of binary digits in the significand (including the
%            implicit leading bit),
%     emins  exponent of xmins,
%     emin:  exponent of xmin,
%     emax:  exponent of xmax.
%   where prec is one of 
%    'b', 'bfloat16'           - bfloat16,
%    'h', 'half', 'fp16'       - IEEE half precision,
%    's', 'single', 'fp32'     - IEEE single precision,
%    'd', double', 'fp64'      - IEEE double precision (the default),
%    'q', 'quadruple', 'fp128' - IEEE quadruple precision.
%   For all these arithmetics the floating-point numbers have the form
%   s * 2^e * d_0.d_1d_2...d_{t-1} where s = 1 or -1, e is the exponent
%   and each d_i is 0 or 1, with d_0 = 1 for normalized numbers.
%   With no input and output arguments, FLOAT_PARAMS prints a table showing
%   all the parameters for all five precisions.
%   Note: xmax and xmin are not representable in double precision for
%   'quadruple'.

%   Author: Nicholas J. Higham.

% References:
% [1] IEEE Standard for Floating-Point Arithmetic, IEEE Std 754-2008 (revision 
%      of IEEE Std 754-1985), 58, IEEE Computer Society, 2008; pages 8,
%      13. https://ieeexplore.ieee.org/document/4610935
% [2] Intel Corporation, BFLOAT16---hardware numerics definition,  Nov. 2018, 
%      White paper. Document number 338302-001US.
%      https://software.intel.com/en-us/download/bfloat16-hardware-numerics-definition}.
% [3] https://stackoverflow.com/questions/44873802/what-is-tf-bfloat16-truncated-16-bit-floating-point
% [4] https://en.wikipedia.org/wiki/Bfloat16_floating-point_format

if nargin < 1 && nargout < 1
   precs = 'bhsdq';
   fprintf(['        u        xmins       xmin       xmax    p    emin    ' ...
            'emax\n'])
   fprintf('    -------------------------------------------------------------\n')
   for j = 1:length(precs)
      [u,xmins,xmin,xmax,p,emin,emax] = float_params(precs(j));
      fprintf('%s: %9.2e  %9.2e  %9.2e  %9.2e  %3.0f  %5.0f  %6.0f\n',...
               precs(j),u,xmins,xmin,xmax,p,emin,emax)
   end
   clear u, return
end   

if nargin < 1, prec = 'd'; end

if ismember(prec, {'h','half','fp16'})
    % Significand: 10 bits plus 1 hidden. Exponent: 5 bits.
    p = 11; emax = 15;
elseif ismember(prec, {'b','bfloat16'})
    % Significand: 7 bits plus 1 hidden. Exponent: 8 bits.
    p = 8; emax = 127;  
elseif ismember(prec, {'s','single','fp32'})
    % Significand: 23 bits plus 1 hidden. Exponent: 8 bits.
    p = 24; emax = 127;
elseif ismember(prec, {'d','double','fp64'})
    % Significand: 52 bits plus 1 hidden. Exponent: 11 bits.
    p = 53; emax = 1023;
elseif ismember(prec, {'q','quadruple','fp128'})
    % Significand: 112 bits plus 1 hidden. Exponent: 15 bits.
    p = 113; emax = 16383;
else 
    error('Unrecognized argument')
end

emin = 1-emax;            % Exponent of smallest normal number.
emins = emin + 1 - p;     % Exponent of smallest subnormal number.
xmins = 2^emins;
xmin = 2^emin;
xmax = 2^emax * (2-2^(1-p));
u = 2^(-p);
