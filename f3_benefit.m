function [ f ] = f3_benefit( t )
global n coe RpgL RpoL RpwL Po Pg Pw Pig;
y = 0;
for i = 1:n+1
    y = y - coe(3,i)*t^(n+1-i);            % between 'y' and 'coe' it must be '-', because ...
    % ...the f1 is to solve the minimum
end
if y >= 0
    error('y must be < 0')
end
% watch out!!! y is negetive!!!
netgasproduction = -y*RpgL(3)/10000 + t - t;     % unit: 10e4 m3/d
if netgasproduction < 0
    error('RpgL(3) is so small that gas production of 3# well is smaller than gas injection of that');
else
    f = -( (-y*Po*RpoL(3)) + (netgasproduction*Pg) - (-y*Pw*RpwL(3)) - Pig*t );       % units of f, y, RpgL and t are...
...yuan/d, m3/d, m3/m3 and 10e4/d respectively
end
if f >= 0
    disp('income of 3# well is negetive')
end
end