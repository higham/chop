%DEMO_HARMONIC  
%   This function evaluates the harmonic sum \sum_{i=1}^{\infty} 1/i
%   in different arithmetics and with different rounding modes.

rng(1)
fprintf('Format  Round mode      Sum      No. terms\n')
fprintf('------------------------------------------\n')
for p = 0:2

clear options

switch p
  case 0
    prec = 'custom'; 
    % "fp8": significand 4 bits plus 1 hidden, exponent 3 bits.
    t = 5; emax = 3;
    options.params = [t emax];
  case 1, prec = 'bfloat16';
  case 2, prec = 'fp16';
end

options.format = prec;

for i = 1:6

    options.round = i;
    % Initialize: subsequent calls chop(x) reuse options.
    chop([],options)

    s = 0; n = 1;

    while true
        sold = s;
        s = chop(s + chop(1/n));
        if s == sold, break, end;
        n = n + 1;
    end 
    prec1 = [prec char(32*ones(1,7))];
    fprintf('%s     %1.0f       %9.4e     %g\n',prec1(1:8),i, s,n)

end
end