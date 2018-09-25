function [ f ] = objfunc_production(t)
global wellquantities n coe ;
f = 0;
for i = 1:wellquantities
    for j = 1:n+1
        f = f - coe(i,j)*t(i)^(n+1-j);            % between 'f' and 'coe' it must be '-', because ...
        % ...the objfunc is to solve the minimum
    end
    if f >= 0
    error('f must be < 0')
    end
end 
end