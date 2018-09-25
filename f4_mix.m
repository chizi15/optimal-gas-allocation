function [ f ] = f4_mix( t )
global n coe RpgL RpoL RpwL Po Pg Pw Pig w1 w2 ymin;
%% 
f1 = 0;
for i = 1:n+1
    f1 = f1 - coe(4,i)*t^(n+1-i);            % between 'y' and 'coe' it must be '-', because ...
        % ...the f1 is to solve the minimum
end
if f1 >= 0
    error('f1 must be < 0')
end

%% 
% watch out!!! f1 is negetive!!!
netgasproduction = -f1*RpgL(4)/10000 +t - t;     % unit: 10e4 m3/d
if netgasproduction < 0
    error('RpgL(4) is so small that gas production of 4# well is smaller than gas injection of that');
end
f2 = -( (-f1*Po*RpoL(4)) + (netgasproduction*Pg) - (-f1*Pw*RpwL(4)) - Pig*t );       % units of f, y, RpgL and t are...
...dollar/d, STB/d, scf/STB and mmscf/d respectively
if f2 >= 0
    disp('income of 4# well is negetive')
end

%% 
f = -( w1*f1/ymin(1,4) + w2*f2/ymin(2,4) ); 
if f >= 0
    disp('mix index of 4# well is negetive')
end
end