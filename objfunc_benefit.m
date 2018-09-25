function [ f ] = objfunc_benefit(t)
% t(i) is gas injection rate of each well,10^4 m^3/d; y is total liquid
% production rate,m^3/d; f is total economic value, yuan.
global wellquantities n coe RpgL RpoL RpwL Po Pg Pw Pig;
a = zeros(wellquantities,1);
y = 0;
for i = 1:wellquantities
    for j = 1:n+1
        y = y - coe(i,j)*t(i)^(n+1-j);            % between 'f' and 'coe' it must be '-', because ...
        % ...the objfunc is to solve the minimum
        a(i) = y;         % a(i) is always negetive
    end
    if y >= 0
        error('y must be < 0')
    end
    
    RpwL(i) = 1 - RpoL(i);
    if i == 1
        netgasproduction = -a(i)*RpgL(i)/10000 + t(i) - t(i);     % unit: 10e4 m3/d
        if netgasproduction < 0
            error('net gas production of 1# well can not be negative anyway');
        else
            f = -( ( -a(i)*Po*RpoL(i) ) + ( netgasproduction*Pg ) - ( -a(i)*Pw*RpwL(i) ) - Pig*t(i) );
        end
    elseif i > 1
        netgasproduction = ( (-a(i)) - (-a(i-1)) )*RpgL(i)/10000 + t(i) - t(i);
        if netgasproduction < 0
            error('net gas production of any well can not be negative anyway. if this error occurs, initial points x0 are wrong possible.');
        else
            f = -( (-f) + ( (-a(i)) - (-a(i-1)) )*Po*RpoL(i) + netgasproduction*Pg - ( (-a(i)) - (-a(i-1)) )*Pw*RpwL(i) - Pig*t(i) );
        end
    end
    if f >= 0
    disp('total income of all wells is negetive');
    end
end
end