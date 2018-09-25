function [ y ] = f4_production( t )
global n coe;
y = 0;
for i = 1:n+1
    y = y - coe(4,i)*t^(n+1-i);       % between 'y' and 'coe' it must be '-', because ...
        % ...the f1 is to solve the minimum
end    
if y >= 0
    error('y must be < 0')
end
end