function [x,minf] = minNT(f,x0,var,eps)
format long;
if nargin == 3
    eps = 1.0e-12;
end
[row,column] = size(x0);
if row ~= 1 || column == 1
    error('x0 must be row vector')
end

tol = 1;
x0 = transpose(x0);
gradf = jacobian(f,var);
hesf = hessian(f,var);

while tol >= eps
    v  = subs(gradf,var,x0);
    pv = subs(hesf,var,x0);
    p = -inv(pv)*transpose(v);        % because v must be a column vector
    p = double(p);
    x1 = x0 + p;
    x0 = x1;
    tol = norm(v);
end
x = x1;
minf = subs(f,var,x);
format short;