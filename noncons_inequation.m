function [c,ceq] = noncons_inequation(t)
global wellquantities n coe RpgL Dg DL                        % injectiongastype;
%% RpgL
% t(i) is injection gas rate of each well(10^4 m^3/d)
prodgas = 0;
for i = 1:wellquantities
    for j = 1:n+1
        % prodgas is production gas rate, t(i) is injection gas rate
        prodgas = prodgas + coe(i,j)*t(i)^(n+1-j)*RpgL(i);   % between 't(i)' and 'coe'...
        %...it must be '+', because the Dg constraint is to solve the
        %maximum£» the unit of 'prodgas' is m3/d.
    end
    if prodgas <= 0
        error('prodgas must be > 0')
    end
    prodgas = prodgas + t(i)*10000;
end
if prodgas/10000 < sum(t)
    error('total gas production rate is smaller than total gas injectioin rate');
else 
    tolgas = prodgas/10000;
end
% if injectiongastype == 1
%     if prodgas/10000 > sum(t)
%         tolgas = prodgas/10000;
%     else
%         tolgas = sum(t);
%     end
% elseif  injectiongastype == 2
%     tolgas = prodgas/10000 + sum(t);
% end

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