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


func = subs(dxn, x, c) / factorial(n1) * x^6;
func
%func =
%-x^6*((atan(c)*sin(c))/720 - cos(c)/(18*(c^2 + 1)^2) - cos(c)/(5*(c^2 + 1)^3) - 
%cos(c)/(120*(c^2 + 1)) + (c*sin(c))/(24*(c^2 + 1)^2) + (c*sin(c))/(2*(c^2 + 1)^3) + 
%(c*sin(c))/(c^2 + 1)^4 + (2*c^2*cos(c))/(9*(c^2 + 1)^3) + (12*c^2*cos(c))/(5*(c^2 + 1)^4) - 
%(16*c^4*cos(c))/(5*(c^2 + 1)^5) - (c^3*sin(c))/(c^2 + 1)^4 - (16*c^3*sin(c))/(3*(c^2 + 1)^5) + 
%(16*c^5*sin(c))/(3*(c^2 + 1)^6))
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
  upper_bounds(i) = (upper_bound);
end
format long G
upper_bounds
% upper_bounds =
%       0.00411999405243267      2.63003356583299e-07      2.12372576358679e-19

