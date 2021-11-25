syms x
syms c
n = 5;
n1=n+1;
x0=0;
range1 = [-0.5, 0.5];
range2 = [-0.1, 0.1];
range3 = [-0.001, 0.001];

f(x) = atan(x)*sin(x);
dxn = diff(f(x), n1);
xn = @(x) x^(n1);

%ser = series(dxn, x, n1) + rn;
%ser
func = subs(dxn, x, c) / factorial(n1) * x^6;

neg_function = -(abs(dxn));
neg_lab_function = matlabFunction(neg_function);

xn_lam_neg =@(x) -(abs(xn(x)));

Ms = [];
Ns = [];

 [xmin, fmin] = fminbnd(neg_lab_function, range1(1), range2(2));
 Ms(1) = (abs(fmin));
 [xmin, fmin] = fminbnd(neg_lab_function, range2(1), range2(2));
  Ms(2) = (abs(fmin));
 [xmin, fmin] = fminbnd(neg_lab_function, range3(1), range3(2));
  Ms(3) = (abs(fmin));

 [xmin, fmin] = fminbnd(xn_lam_neg, range1(1), range1(2));
  Ns(1) = (abs(fmin));
 [xmin, fmin] = fminbnd(xn_lam_neg, range2(1), range2(2));
  Ns(2) = (abs(fmin));
 [xmin, fmin] = fminbnd(xn_lam_neg, range3(1), range3(2));
  Ns(3) = (abs(fmin));
upper_bounds = [];

for i= 1:3
  upper_bound = Ms(i) / factorial(n1) * Ns(i);
 % upper_bound
  upper_bounds(i) = (upper_bound);
end
format long G
upper_bounds


