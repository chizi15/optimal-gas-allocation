function [ f ] = f7_mix( t )
global n coe RpgL RpoL RpwL Po Pg Pw Pig w1 w2 ymin;
%%
f1 = 0;
for i = 1:n+1
    f1 = f1 - coe(7,i)*t^(n+1-i);            % between 'y' and 'coe' it must be '-', because ...
    % ...the f1 is to solve the minimum
end
if f1 >= 0
    error('f1 must be < 0')
end

% watch out!!! f1 is negetive!!!
f2 = -( (-f1*Po*RpoL(7)) + (-f1*Pg*RpgL(7)/10000) - (-f1*Pw*RpwL(7)) - Pig*t );       % units of f, y, RpgL and t are...
...dollar/d, STB/d, scf/STB and mmscf/d respectively
    if f2 >= 0
    disp('income of 7# well is negetive')
    end
    
    %%
    f = -( w1*f1/ymin(1,7) + w2*f2/ymin(2,7) );    % f1,f2,ymin_production,ymin_benefit are negetive, so f is negetive
    if f >= 0
        disp('mix index of 7# well is negetive')
    end
end