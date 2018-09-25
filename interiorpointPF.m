function [x,minf] = interiorpointPF(f,x0,inequality,r0,c,variable,eps)
% inequality is column vector
if ~isvector(x0) || ~isvector(inequality) || ~isvector(variable)
    error('Input must be a vector')
end
gx0 = subs(inequality,variable,x0);
if gx0 <= 0
    disp('初始点必须满足不等式约束！');
    x = NaN;
    minf = NaN;
    return;
end

if r0 <= 0
    disp('初始障碍因子必须大于0！');
    x = NaN;
    minf = NaN;
    return;
end

if c >= 1 || c <= 0
    disp('缩小系数必须在[0,1]之间！');
    x = NaN;
    minf = NaN;
    return;
end

if nargin == 7
    eps = 1.0e-12;
end

FE = 0;
for i=1:length(inequality)
    FE = FE + 1/inequality(i);
end

[row,column] = size(x0);
if row ~= 1 || column == 1
    error('x0 must be row vector')
end
x1 = transpose(x0);

while true
    FF = r0*FE;
    SumF = f + FF ;
    x2 = minNT(SumF,transpose(x1),variable);
    if norm(x2 - x1) < eps && r0*FE < eps
        x = x2;
        break;
    else
        r0 = c*r0;
        x1 = x2;
    end
end
minf = subs(f,variable,x);