function [ y ] = f1_production( t )
global coe n;
y = 0;
for i = 1:n+1
    y = y - coe(1,i)*t^(n+1-i);            % between 'y' and 'coe' it must be '-', because ...
    % ...the f1_production is to solve the minimum
end
if y >= 0
    error('y must be < 0')
    % matlab 中创建任何函数都是在其给定的定义域范围内创建的，如果没有给定定义域，则默认为R；
    % 这里虽然没给出自变量t的范围，但在调用它的fmincon函数中定义了的，所以其自变量范围是【XL（1），XR（1）】
end
end