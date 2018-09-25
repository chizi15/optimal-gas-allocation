% 用最小二乘法拟合三次多项式的，求得的三次多项式的系数a，b，c，d
syms a b c d ;
data(1,1) = 140; data(1,2) = 0.126;
data(2,1) = 400; data(2,2) = 0.182;
data(3,1) = 900; data(3,2) = 0.240;
data(4,1) = 1600; data(4,2) = 0.333;
M = ( a*data(1,1)^3 + b*data(1,1)^2 + c*data(1,1) + d - data(1,2) )^2 ...
+ ( a*data(2,1)^3 + b*data(2,1)^2 + c*data(2,1) + d - data(2,2) )^2 ...
+ ( a*data(3,1)^3 + b*data(3,1)^2 + c*data(3,1) + d - data(3,2) )^2 ...
+ ( a*data(4,1)^3 + b*data(4,1)^2 + c*data(4,1) + d - data(4,2) )^2;
eqna = diff(M,a) == 0;
eqnb = diff(M,b) == 0;
eqnc = diff(M,c) == 0;
eqnd = diff(M,d) == 0;
sol = solve([eqna,eqnb,eqnc,eqnd],[a,b,c,d]);
disp( sol.a ); disp( sol.b ); disp( sol.c ); disp( sol.d );
x(1) = 2000; x(2) = 1200; x(3) = 600; x(4) = 200;
y = zeros(4,1);
for i = 1:4
   y(i) = sol.a*x(i)^3 + sol.b*x(i)^2 + sol.c*x(i) + sol.d;
   dots = [x(i) y(i)];
   format long;
   disp( dots );
end

%% 
% 用最小二乘法拟合直线，求得直线的系数a，b
syms a b;
data(1,1) = 140; data(1,2) = 0.126;
data(2,1) = 400; data(2,2) = 0.182;
data(3,1) = 900; data(3,2) = 0.240;
data(4,1) = 1600; data(4,2) = 0.333;
M = ( a*data(1,1) + b - data(1,2) )^2 + ( a*data(2,1) + b - data(2,2) )^2 ...
    + ( a*data(3,1) + b - data(3,2) )^2 + ( a*data(4,1) + b - data(4,2) )^2;
eqna = diff(M,a) == 0;
eqnb = diff(M,b) == 0;
sol = solve([eqna,eqnb],[a,b]);
disp( sol.a ); disp( sol.b );
x(1) = 2000; x(2) = 1200; x(3) = 600; x(4) = 200;
y = zeros(4,1);
for i = 1:4
   y(i) = sol.a*x(i) + sol.b;
   dots = [x(i) y(i)];
   format long;
   disp( dots );
end

%% 

y1 = 300*( (1+0.12)^2.5 - 1 ) + 600*( (1.12)^1.5 - 1 ) + 400*( (1.12)^0.5 - 1);
%% 

f = 0.05; n = 30; Inv = 1000; % f为资金年变化率，f取除零外的一切实数；n为资金经历的整年数,非负整数，单位：年; Inv为投资额或贷款额，单位：万；
% y_season = f*( 1 + f )^n * ( 1 + f )^( 1/(4*2) ) / ( 4*( ( 1 + f )^( 2/(4*2) ) - 1 ) ) - ( 1 + f )^(n+0.5);   % 在最后或第一个非整计息年，每季度中均等注入（I/4）与年中一次性注入I的系数差值
% y_month = f*( 1 + f )^n*( 1 + f )^( 1/(12*2) ) / ( 12*( ( 1 + f )^( 2/(12*2) ) - 1 ) ) - ( 1 + f )^(n+0.5);  % 在最后或第一个非整计息年，每月中均等注入（I/12）与年中一次性注入I的系数差值
% y_day = f*( 1 + f )^n*( 1 + f )^( 1/(365*2) ) / ( 365*( ( 1 + f )^( 2/(365*2) ) - 1 ) ) - ( 1 + f )^(n+0.5);  % 在最后或第一个非整计息年，每日中均等注入（I/365）与年中一次性注入I的系数差值
% y_sm = y_season - y_month;
% y_md = y_month - y_day;

PF_moy = Inv*( (1+f)^(n+0.5) - 1 );
PF_s = Inv*( (1+f)^n*(1+f)^( 1/(4*2) )*f / ( 4*( (1+f)^( 2/(4*2) ) -1 ) ) - 1 );
PF_m = Inv*( (1+f)^n*(1+f)^( 1/(12*2) )*f / ( 12*( (1+f)^( 2/(12*2) ) -1 ) ) - 1 );
PF_d = Inv*( (1+f)^n*(1+f)^( 1/(365*2) )*f / ( 365*( (1+f)^( 2/(365*2) ) -1 ) ) - 1 );
D_dnm = PF_d - PF_m;
D_mns = PF_m - PF_s;
D_sny = PF_s - PF_moy;
D_mny = PF_m - PF_moy;
D_dny = PF_d - PF_moy;

%% 
x = (-1+13^0.5)/(-4);
y = x; z = 2*x;
lamida1 = 3*(1 - 13^0.5/13);
lamida2 = 5/2 - 29*3^0.5/26;
f1 = 2*x -2*lamida1*x + lamida2;

%% 
f = ( 2*( (-1+3^0.5)/2 )^2 + (2-3^0.5)^2 )^0.5;
f = ( 2*( (-1-3^0.5)/2 )^2 + (2-3^0.5)^2 )^0.5;
