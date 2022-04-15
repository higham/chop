function test_chop
%TEST_CHOP Test the chop function.
%   The tests are for single precision and fp16.

clear chop fp options options2 assert_eq

rng(1)

n = 0;
uh = 2^(-11);  % Unit roundoff for fp16.
pi_h = 6432*uh; % fp16(pi)

% Check handling of defaults and persistent variable.
fp.format = 'bfloat16'; [c,options] = chop(pi,fp);
assert_eq(fp.format,options.format)
assert_eq(options.subnormal,0)

fp.format = []; [c,options] = chop(pi,fp);
assert_eq(options.format,'h')  % Check default;

fp.subnormal = 0; [c,options] = chop(pi,fp);
assert_eq(options.subnormal,0)

fp.subnormal = []; [c,options] = chop(pi,fp);
assert_eq(options.subnormal,1)  % Check default;

fp.round = []; [c,options] = chop(pi,fp);
assert_eq(options.round,1)  % Check default.

fp.flip = []; [c,options] = chop(pi,fp);
assert_eq(options.flip,0)  % Check no default.

fp.explim = []; [c,options] = chop(pi,fp);
assert_eq(options.explim,1)  % Check default.

fp.explim = 0; [c,options] = chop(pi,fp);
assert_eq(options.explim,0)  % Check no default.

clear chop fp options
fp.flip = 1; [~,options] = chop([],fp);
assert_eq(options.format,'h')
assert_eq(options.round,1)
assert_eq(options.subnormal,1)

clear chop fp options
% check all default options
fp.format = []; fp.subnormal = [];
fp.round = []; fp.flip = [];
fp.p = [];
[c,options] = chop(pi,fp);
assert_eq(options.format,'h')
assert_eq(options.subnormal,1)
assert_eq(options.round,1)
assert_eq(options.flip,0)
assert_eq(options.p,0.5)
% % Takes different path from previous test since fpopts exists.
% fp.subnormal = 0;
% fp.format = []; [c,options] = chop(pi,fp);
% assert_eq(options.format,'h')

% Check flip output.
clear chop fp
fp.flip = 1; fp.format = 'd';
c = ones(8,1);
d = chop(c,fp); assert_eq(norm(d-c,1)>0,true);
d = chop(c',fp); assert_eq(norm(d-c',1)>0,true);
fp.p = 0; % No bits flipped.
d = chop(c,fp); assert_eq(d,d);
fp.p = 1; % All bits flipped.
d = chop(c,fp); assert_eq(all(d ~= c),true);

clear chop
[~,fp] = chop;
assert_eq(fp.subnormal,1)
assert_eq(fp.format,'h')
[c,options] = chop(pi);
assert_eq(options.format,'h')
assert_eq(options.subnormal,1)
assert_eq(options.round,1)
assert_eq(options.flip,0)
assert_eq(options.p,0.5)

fp.format = 'd'; [c,options] = chop(pi,fp);
assert_eq(options.format,'d')
assert_eq(options.subnormal,1)
assert_eq(options.params, [53 1023])
[~,fp] = chop;
assert_eq(fp.format,'d')
assert_eq(fp.subnormal,1)
assert_eq(fp.params, [53 1023])

clear fp
fp.format = 'bfloat16'; [c,options] = chop(pi,fp);
assert_eq(options.format,'bfloat16')
assert_eq(options.subnormal,0)
assert_eq(options.params, [8 127])
[~,fp] = chop;
assert_eq(fp.format,'bfloat16')
assert_eq(fp.subnormal,0)
assert_eq(fp.params, [8 127])

clear chop
[~,fp] = chop;
fp.format = 'b'; [c,options] = chop(pi,fp);
assert_eq(options.subnormal,1) % No subnormals only if that field was empty.

% Check these usages do not give an error.
c = chop([]);
chop([]);
chop([],fp);
chop(1,fp);
c = chop(1,fp);

% Test matrix.
options.format = 'b';
A = magic(4);
C = chop(A,options);
assert_eq(A,C);
B = A + randn(size(A))*1e-12;
C = chop(B,options);
assert_eq(A,C);
A2 = hilb(6); C = chop(A2);

options.format = 'c';
options.params = [8 127];  % bfloat16
C1 = chop(A,options);
assert_eq(A,C1);
C2 = chop(B,options);
assert_eq(A,C2);
assert_eq(C,chop(A2));

clear options
options.format = 'c';
options.params = [11 15];  % h
options2.format = 'h';
A = hilb(6);
[X1,opt] = chop(A,options);
[X2,opt2] = chop(A,options2);
assert_eq(X1,X2)
% assert_eq(chop(A,options),chop(A,options2));

% Row vector
clear options
options.format = 'h';
A = -10:10;
C = chop(A,options);
assert_eq(A,C);
B = A + randn(size(A))*1e-12;
C = chop(B,options);
assert_eq(A,C);

% Column vector
options.format = 's';
A = (-10:10)';
C = chop(A,options);
assert_eq(A,C);
B = A + A.*rand(size(A))*1e-14;  % Keep 0 as 0.
C = chop(B,options);
assert_eq(A,C);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main loop: test single, half and quarter formats.
for i = 1:4
clear chop fp options

if i == 1
   % Single precision tests.
   [u,xmins,xmin,xmax,p,emins,emin,emax] = float_params('single');
   options.format = 's';
elseif i == 2
   % Half precision tests.
   [u,xmins,xmin,xmax,p,emins,emin,emax] = float_params('half');
   options.format = 'h';
elseif i == 3
   % Half precision tests.
   [u,xmins,xmin,xmax,p,emins,emin,emax] = float_params('fp8-e4m3');
   options.format = 'fp8-e4m3';
elseif i == 4
   % Half precision tests.
   [u,xmins,xmin,xmax,p,emins,emin,emax] = float_params('fp8-e5m2');
   options.format = 'fp8-e5m2';
end
options.subnormal = 0;

x = pi;
if i == 1
   y = double(single(x));
elseif i == 2
   y = pi_h; % double(fp16(x));
elseif i == 3
   y = 3.25;
elseif i == 4
   y = 3.0;
end
c = chop(x,options);
assert_eq(c,y);
x = -pi;
c = chop(x,options);
assert_eq(c,-y);

% Next number power of 2.
y = 2^10;
if i == 1
   dy = double(eps(single(y)));
elseif i == 2
   dy = 2*y*uh; % double(eps(fp16(y)));
elseif i == 3
   y = 2^4;
   dy = 2*y*u;
elseif i == 4
   y = 2^4;
   dy = 2*y*u;
end
x = y + dy;
c = chop(x,options);
assert_eq(c,x)

% Number just before a power of 2.
x = y - dy;
c = chop(x,options);
assert_eq(c,x)

% Next number power of 2.
y = 2^(-4);
if i == 1
   dy = double(eps(single(y)));
elseif i == 2
   dy = 2*y*uh; % double(eps(fp16(y)));
elseif i == 3
   dy = 2*y*u;
elseif i == 4
   dy = 2*y*u;
end
x = y + dy;
c = chop(x,options);
assert_eq(c,x)

% Check other rounding options
for rmode = 1:6
    options.round = rmode;
    x = y + (dy*10^(-3));
    c = chop(x,options);
    if options.round == 2
        assert_eq(c,y+dy) % Rounding up.
    elseif options.round >= 5
        % Check rounded either up or down.
        if c ~= y+dy
           assert_eq(c,y);
        end
    else
        assert_eq(c,y);
    end
end

% Overflow tests.
for j = 1:6
  options.round = j;
  x = xmax;
  c = chop(x,options);
  assert_eq(c,x)
end

% Infinities tests.
for j = 1:6
  options.round = j;
  x = inf;
  c = chop(x,options);
  assert_eq(c,x)
  c = chop(-x,options);
  assert_eq(c,-x)
end

% IEEE 754-2019, page 27: rule for rounding to infinity.
% Round to nearest
options.round = 1; % reset the rounding mode to default
x = 2^emax * (2-(1/2)*2^(1-p));  % Round to inf.
c = chop(x,options);
assert_eq(c,inf)
c = chop(-x,options);
assert_eq(c,-inf)

x = 2^emax * (2-(3/4)*2^(1-p));  % Round to realmax.
c = chop(x,options);
assert_eq(c,xmax)
c = chop(-x,options);
assert_eq(c,-xmax)

% Round towards plus infinity
options.round = 2;
x = 2^emax * (2-(1/2)*2^(1-p));
c = chop(x,options);
assert_eq(c,inf)
c = chop(-x,options);
assert_eq(c,-xmax)

% Round towards minus infinity
options.round = 3;
c = chop(x,options);
assert_eq(c,xmax)
c = chop(-x,options);
assert_eq(c,-inf)

% Round towards zero
options.round = 4;
c = chop(x,options);
assert_eq(c,xmax)
c = chop(-x,options);
assert_eq(c,-xmax)

% Round to nearest.
options.round = 1; % reset the rounding mode to default
if i == 2
   x = 1 + 2^(-11);
   c = chop(x,options);
   assert_eq(c,1)
end

% Underflow tests.
if i == 1
    delta = double(eps(single(1)));
else
    delta = 2*uh; % double(eps(fp16(1)));
end

options.subnormal = 1;
c = chop(xmin,options); assert_eq(c,xmin)
x = [xmins xmin/2 xmin 0 xmax 2*xmax 1-delta/5 1+delta/4];
c = chop(x,options);
c_expected = [x(1:5) inf 1 1];
assert_eq(c,c_expected)

options.subnormal = 0;
c = chop(xmin,options); assert_eq(c,xmin)
x = [xmins xmin/2 xmin 0 xmax 2*xmax 1-delta/5 1+delta/4];
c = chop(x,options);
c_expected = [0 xmin x(3:5) inf 1 1];
assert_eq(c,c_expected)

% Smallest normal number and spacing between the subnormal numbers.
y = xmin; delta = xmin*2^(1-p);
x = y - delta; % The largest subnormal number.
options.subnormal = 1;
c = chop(x,options);
assert_eq(c,x)
% Round up if subnormals are not supported.
options.subnormal = 0;
c = chop(x,options);
assert_eq(c,xmin)
% Flush subnormals to zero if subnormals are not supported.
options.subnormal = 0;
c = chop(xmins,options);
assert_eq(c,0)

options.subnormal = 1;
x = xmins*8;  % A subnormal number.
c = chop(x,options);
assert_eq(c,x)

% Numbers smaller than smaller representable number.
options.subnormal = 0;
x = xmin / 2;
c = chop(x,options);
assert_eq(c,xmin)
x = -xmin / 2;
c = chop(x,options);
assert_eq(c,-xmin)
x = xmin / 4;
c = chop(x,options);
assert_eq(c,0)
x = -xmin / 4;
c = chop(x,options);
assert_eq(c,0)

options.subnormal = 1;
x = xmins / 2;
c = chop(x,options);
assert_eq(c,0)
x = -xmins / 2;
c = chop(x,options);
assert_eq(c,0)
x = xmins / 4;
c = chop(x,options);
assert_eq(c,0)
x = -xmins / 4;
c = chop(x,options);
assert_eq(c,0)

% Do not limit exponent.
options.explim = 0;
x = xmin/2;  c = chop(x,options); assert_eq(c,x)
x = -xmin/2;  c = chop(x,options); assert_eq(c,x)
x = xmax*2;  c = chop(x,options); assert_eq(c,x)
x = -xmax*2;  c = chop(x,options); assert_eq(c,x)
x = xmins/2; c = chop(x,options); assert_eq(c,x)
x = -xmins/2; c = chop(x,options); assert_eq(c,x)
A = [pi -pi; pi -pi];
C = chop(A,options);
options.explim = 1;
assert_eq(C,chop(A,options));

% Round towards plus infinity
options.round = 2;
options.subnormal = 0;
x = xmin / 2;
c = chop(x,options);
assert_eq(c,xmin)
x = -xmin / 2;
c = chop(x,options);
assert_eq(c,0)

options.subnormal = 1;
x = xmins / 2;
c = chop(x,options);
assert_eq(c,xmins)
x = -xmins / 2;
c = chop(x,options);
assert_eq(c,0)
x = xmins / 4;
c = chop(x,options);
assert_eq(c,xmins)
x = -xmins / 4;
c = chop(x,options);
assert_eq(c,0)

% Round towards minus infinity
options.round = 3;
options.subnormal = 0;
x = xmin / 2;
c = chop(x,options);
assert_eq(c,0)
x = -xmin / 2;
c = chop(x,options);
assert_eq(c,-xmin)

options.subnormal = 1;
x = xmins / 2;
c = chop(x,options);
assert_eq(c,0)
x = -xmins / 2;
c = chop(x,options);
assert_eq(c,-xmins)
x = xmins / 4;
c = chop(x,options);
assert_eq(c,0)
x = -xmins / 4;
c = chop(x,options);
assert_eq(c,-xmins)

% Round towards zero.
options.round = 4;
options.subnormal = 0;
x = xmin / 2;
c = chop(x,options);
assert_eq(c,0)
x = -xmin / 2;
c = chop(x,options);
assert_eq(c,0)

options.subnormal = 1;
x = xmins / 2;
c = chop(x,options);
assert_eq(c,0)
x = -xmins / 2;
c = chop(x,options);
assert_eq(c,0)
x = xmins / 4;
c = chop(x,options);
assert_eq(c,0)
x = -xmins / 4;
c = chop(x,options);
assert_eq(c,0)

end % for i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear options

% Test rounding with CHOP versus native rounding.
options.format = 's';
m = 100; y = zeros(3,n); z = y;
for i = 1:m
x = randn;
options.round = 2; y(i,1) = chop(x,options);
options.round = 3; y(i,2) = chop(x,options);
options.round = 4; y(i,3) = chop(x,options);
% Use undocumented function to set rounding mode in MATLAB.
feature('setround',inf), z(i,1) = single(x);
feature('setround',-inf), z(i,2) = single(x);
feature('setround',0), z(i,3) = single(x);
end
assert_eq(y,z)
feature('setround',0.5) % Back to round to nearest.

% Double precision tests.
[u,xmins,xmin,xmax,p,emins,emin,emax] = float_params('d');
options.format = 'd';
x = [1e-309 1e-320 1 1e306];  % First two entries are subnormal.
c = chop(x,options);
assert_eq(c,x)
options.subnormal = 0;
c = chop(x,options);
assert_eq(c,[0 0 x(3:4)])

options.format = 'd'; options.subnormal = 0; chop([],options)
a = chop(pi); assert_eq(a,pi)
options.format = 'd'; options.subnormal = 1; chop([],options)
a = chop(pi); assert_eq(a,pi)

x = pi^2;
clear options
options.format = 'd';
y = chop(x,options);  % Should not change x.
assert_eq(x,y);
options.round = 2;
y = chop(x,options);  % Should not change x.
assert_eq(x,y);
options.round = 3;
y = chop(x,options);  % Should not change x.
assert_eq(x,y);
options.round = 4;
y = chop(x,options);  % Should not change x.
assert_eq(x,y);

% Test on single inputs.
clear options
ps = single(pi);
pd = double(ps);
options.format = 'b';
ys = chop(ps,options);
assert_eq(isa(ys,'single'),true)
yd = chop(pd);
assert_eq(double(ys),yd)

options.format = 'h'; options.round = 2;
as = single(rand(n,1)); ad = double(as);
delta = single(rand(n,1));
cd = chop(ad + 1e-5*double(delta),options);
cs = chop(as + 1e-5*delta,options);
assert_eq(cd,double(cs));

options.format = 'c';
options.params = [11 5];
temp1 = chop(single(pi),options);
options.format = 'h';
temp2 = chop(single(pi),options);
assert_eq(temp1,temp2)

% Test base 2 logarithm
options.format = 'h';
options.round = 4;
x = single(2^-3 * (sum(2.^(-[0:23]))));
assert_eq(chop(x,options), single(2^-3 * (sum(2.^(-[0:10])))))

x = 2^-3 * (sum(2.^(-[0:52])));
assert_eq(chop(x,options), 2^-3 * (sum(2.^(-[0:10]))))

options.format = 's';
x = single(2^-3 * (sum(2.^(-[0:23]))));
assert_eq(chop(x,options), x)

x = 2^-3 * (sum(2.^(-[0:52])));
assert_eq(chop(x,options), 2^-3 * (sum(2.^(-[0:23]))))

options.format = 'd';
x = 2^-3 * (sum(2.^(-[0:52])));
assert_eq(chop(x,options), x)

options.round = 1;
temp = 0;
try
    options.format = 'c';
    options.params = [12 5];
    temp = chop(single(pi),options); % Error - double rounding!
catch
end
assert_eq(temp,0)
try
    options.format = 'c';
    options.params = [26 9];
    temp = chop(pi,options); % Error - double rounding!
catch
end
assert_eq(temp,0)
try
    temp = chop(complex(1,1)); % Error - complex data!
catch
end
assert_eq(temp,0)

fprintf('All tests successful!\n')

clear chop fp options options2 assert_eq

%%%%%%%%%%%%%%%%%%%%%%%
function assert_eq(a,b)
% if isempty(n), n = 0; end  % First call.
n = n+1;
if ~isequal(a,b)
   error('Failure')
end
fprintf('Test %g succeeded.\n', n )
end

end