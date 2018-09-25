function [c,ceq] = noncons(t)
global wellquantities n coe RpgL Dg DL;
%% 
% t(i) is injection air rate of each well(mmscfd)
prodgas = 0;
for i = 1:wellquantities
    for j = 1:n+1
        % prodgas is production gas rate, t(i) is injection gas rate
        prodgas = prodgas + coe(i,j)*t(i)^(n+1-j)*RpgL(i);   % between 't(i)' and 'coe'...
        %...it must be '+', because the Dg constraint is to solve the maximum
    end
    if prodgas <= 0
        error('prodgas must be > 0')
    end
    tolgas = prodgas/10000 + t(i);    % prodgas is the sum of production gas of all wells, and t(i) is gas injection rate of each well
    if tolgas <= 0 || tolgas > 10000
        error('tolgas must be > 0 and < 10000')
    end
end
c(1) = tolgas - Dg;     % tolgas is the sum of injection gas rate and production gas rate, c(1) must be <= 0
%% 
prodliquid = 0;
for i = 1:wellquantities
    for j = 1:n+1
        % prodliquid is production liquid rate
        prodliquid = prodliquid + coe(i,j)*t(i)^(n+1-j);      % between 'c(2)' and 'coe'...
        %...it must be '+', because the DL constraint is to solve the maximum
    end
    if prodliquid <= 0
        error('prodliquid must be > 0')
    end
end
c(2) = prodliquid - DL;
%% 
ceq = [];
end