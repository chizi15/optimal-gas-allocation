function [ f ] = f7_benefit( t )
global n coe RpgL RpoL RpwL Po Pg Pw Pig;
y = 0;
for i = 1:n+1
    y = y - coe(7,i)*t^(n+1-i);            % between 'y' and 'coe' it must be '-', because ...
    % ...the f1 is to solve the minimum
end
if y >= 0
    error('y must be < 0')
end
% watch out!!! y is negetive!!!
f = -( (-y*Po*RpoL(7)) + (-y*Pg*RpgL(7)/10000) - (-y*Pw*RpwL(7)) - Pig*t );       % units of f, y, RpgL and t are...
...dollar/d, STB/d, scf/STB and mmscf/d respectively
    if f >= 0
    disp('income of 7# well is negetive')
    end
end