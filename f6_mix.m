function [ f ] = f6_mix( t )
global n coe RpgL RpoL RpwL Po Pg Pw Pig w1 w2 ymin;
%% 
f1 = 0;
for i = 1:n+1
    f1 = f1 - coe(6,i)*t^(n+1-i);            % between 'y' and 'coe' it must be '-', because ...
        % ...the f1 is to solve the minimum
end
if f1 >= 0
    error('f1 must be < 0')
end

%% 
% watch out!!! f1 is negetive!!!
f2 = -( (-f1*Po*RpoL(6)) + (-f1*Pg*RpgL(6)/10000) - (-f1*Pw*RpwL(6)) - Pig*t );       % units of f, y, RpgL and t are...
...dollar/d, STB/d, scf/STB and mmscf/d respectively
if f2 >= 0
    disp('income of 6# well is negetive')
end

%% 
f = -( w1*f1/ymin(1,6) + w2*f2/ymin(2,6) );
if f >= 0
    disp('mix index of 6# well is negetive')
end
end