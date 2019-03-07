function test_roundit
%TEST_ROUNDIT Test the roundit function.

clear roundit options assert_eq
n = 0;
rng(1)

y = roundit(pi);     % Check no options case.
assert_eq(y,3)
y = roundit(pi,[]);  % Check empty options case.
assert_eq(y,3)

options.flip = 0;
options.round = 1;
A = [0 1.1 1.5; 1.9 2.4 0.5];
y = roundit(A,options);
assert_eq(y,[0 1 2; 2 2 0])
y = roundit(A',options);
assert_eq(y,[0 1 2; 2 2 0]')

A = [0 -1.1 -1.5; -1.9 -2.4 -0.5];
y = roundit(A,options);
assert_eq(y,[0 -1 -2; -2 -2 0])

options.round = 2;
A = [0 1.1 1.5; 1.9 2.4 0.5];
y = roundit(A,options);
assert_eq(y,[0 2 2; 2 3 1])

A = [0 -1.1 -1.5; -1.9 -2.4 -0.5];
options.round = 2;
y = roundit(A,options);
assert_eq(y,[0 -1 -1; -1 -2 0])

options.round = 3;
A = [0 1.1 1.5; 1.9 2.4 0.5];
y = roundit(A,options);
assert_eq(y,[0 1 1; 1 2 0])

options.round = 3;
A = [0 -1.1 -1.5; -1.9 -2.4 -0.5];
y = roundit(A,options);
assert_eq(y,[0 -2 -2; -2 -3 -1])

options.round = 4;
A = [0 -1.1 -1.5; -2.9 -2 -0.5; 0.5 1.5 3];
y = roundit(A,options);
assert_eq(y,[0 -1 -1; -2 -2 0; 1 2 3])

options.round = 5;
tol = 1e-6;
A = [0 1+tol 1-tol; -tol -1-tol -1+tol];
y = roundit(A,options);
assert_eq(y,[0 1 1; 0 -1 -1])  % Test succeeds with high probability.

options.round = 5;
A = [-0.5 0.5]; Y = [];
for i = 1:1e3
   Y = [Y; roundit(A,options)];
end
assert_eq(true, norm(mean(Y) - 0.5*[-1 1]) <= 5e-2)

z = 2; nsamp = 1e4;
options.round = 5;
for h = [0.1 0.2 0.3 0.4 0.6 0.7 0.8 0.9]
   x = z + h;
   a = 0; b = 0; 
   for i = 1:nsamp
       y = roundit(x,options);
       if y == z
           a = a + 1;
       elseif y == z + 1
           b = b + 1;
       else
           error('Unknown rounding.')
       end
   end
   assert_eq(true, (abs( a/nsamp-(1-h) ) <= 5e-2));
   assert_eq(true, (abs( b/nsamp - h ) <= 5e-2));
end

options.round = 6;
A = [0.2 0.5 0.9]; Y = [];
for i = 1:1e4
   Y = [Y; roundit(A,options)];
end
ratio = sum(Y(:) == 1)/ sum(Y(:) == 0);
assert_eq(true, abs(ratio-1) < 0.1) % Heuristic test.

% Test bit flip case.
options.round = 1;
m = 500;
A = 15*rand(m);
y1 = roundit(A,options); % Now y1 is an integer matrix.
options.flip = 1;
options.t = 4;  % Integers with modulus on [0,15].
y2 = roundit(y1,options); % No rounding, but bits flipped.
prop_changed = sum(sum(y1 ~= y2)) / m^2;
assert_eq(true, abs(prop_changed - 0.5) < 5e-2) % Heuristic test.
options.p = 0.1;
y2 = roundit(y1,options); % No rounding, but bits flipped.
prop_changed = sum(sum(y1 ~= y2)) / m^2;
assert_eq(true, abs(prop_changed - options.p) < 5e-2) % Heuristic test.

% Make sure zero is bit-flipped correctly.
options = rmfield(options,'p');
options.flip = 1; nsamp = 1e4;
j = 0; options.t = 3; options.p = 0.5;
for i = 1:nsamp
    y = roundit(0,options); 
    if y ~= 0, j = j +1; end
end
assert_eq(true, abs(j / nsamp - options.p) <= 1e-2)

fprintf('All tests successful!\n')

%%%%%%%%%%%%%%%%%%%%%%%
function assert_eq(a,b)
n = n+1;
if ~isequal(a,b)
   a, b
   error('Failure')
end
fprintf('Test %g succeeded.\n', n )
end

end