function [ f ] = objfunc_mix( t )
% t(i) is gas injection rate of each well
global wellquantities n coe ytolrealmax_production ytolrealmax_benefit RpgL RpoL RpwL Po Pg Pw Pig w1 w2;
%%
% f1 is the sum of negetive dimensionless liquid production of all wells
f1 = 0;
for i = 1:wellquantities
    for j = 1:n+1
        f1 = f1 - coe(i,j)*t(i)^(n+1-j);            % between 'f1' and 'coe' it must be '-', because ...
        % ...the objfunc is to solve the minimum
    end
    if f1 >= 0
        error('f1 must be < 0')
    end
end
%%
a = zeros(wellquantities,1);
y = 0; f2 = 0;
for i = 1:wellquantities
    for j = 1:n+1
        y = y - coe(i,j)*t(i)^(n+1-j);            % between 'f' and 'coe' it must be '-', because ...
        % ...the objfunc is to solve the minimum
        a(i) = y;
    end
    if y >= 0
        error('y must be < 0')
    end
    RpwL(i) = 1 - RpoL(i);
    if i ==1
        % f2 is the sum of negetive economic benefits of all wells
        f2 = -( ( -a(i)*Po*RpoL(i) ) + ( -a(i)*RpgL(i)/10000 - t(i) + t(i) )*Pg - ( -a(i)*Pw*RpwL(i) ) - Pig*t(i) );
    elseif i > 1
        f2 = -( ( -f2 ) + ( ( -a(i) ) - ( -a(i-1) ) )*Po*RpoL(i) +( ( ( -a(i) ) - ( -a(i-1) ) )*RpgL(i)/10000 - t(i) + t(i) )*Pg- ...
            ( ( -a(i) ) - ( -a(i-1) ) )*Pw*RpwL(i) - Pig*t(i) );
    end
end
if f2 >= 0
    disp('total income of all wells is negetive');
end
%%
    f = w1 * f1 / ytolrealmax_production + w2 * f2 / ytolrealmax_benefit ;       % f1,f2 are negetive; ytolrealmax_production and...
    ...ytolrealmax_benefit are positive, f is negetive.
    if f >= 0
    disp('mix index of all wells is negetive')
    end
end