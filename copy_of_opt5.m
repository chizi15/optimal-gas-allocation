tic
load ('kelameili_5.mat');
global wellquantities n coe RpgL RpoL RpwL Po Pg Pw Pig Dg DL gas_provide...
    ytolrealmax_production ytolrealmax_benefit w1 w2 ymin XL xmax objectmodel injectiongastype

%%
% gas injection & liquid production rate (10^4 m^3/d, m^3/d)
data{1,1} = mmscfd1; data{1,2} = STBd1;
data{2,1} = mmscfd2; data{2,2} = STBd2;
data{3,1} = mmscfd3; data{3,2} = STBd3;
data{4,1} = mmscfd4; data{4,2} = STBd4;
data{5,1} = mmscfd5; data{5,2} = STBd5;
[wellquantities,~] = size(data);

%%
RpoL = [0.85 0.12 0.25 0.35 0.20];  % RpoL is a row vector

for i = 1:wellquantities
    if RpoL(i) < 0 || RpoL(i) > 1
        error('RpoL(i) must be between [0,1]');
    end
end
% evaluate production gas&liquid rate(RpgL) for nonconstraints
% distribute memory for wcut in advance to enhance computation speed
RpgL = zeros(wellquantities,1);     % RpgL is a column vector
% gas&oil rate on surface standard condition, the unit of GOR and RpgL is: m3/m3
% "production gas&oil ratio on surface standard condition, scf/STB"   (not use "scf/STB")

GOR(1) = 1; GOR(2) = 1; GOR(3) = 1; GOR(4) = 1; GOR(5) = 1;  % m3/m3  notice!!! GOR means net gas production...
...and oil production ratio, excluding injection gas amount.

disp('production gas liquid ratio of each well are below:(m3/m3)')
for i = 1:wellquantities
    if GOR(i) <= 0
        error('GOR must be > 0')
    end
    RpgL(i) = RpoL(i) * GOR(i);  % RpoL is a row vector, RpgL is a column vector, no controvertial in MATLAB
    disp(RpgL(i));     % m3/m3
end
RpwL = 1 - RpoL;
m3toSTB = 0.15898;   % unit:m3/STB

RMBTODOLLAR = 7;  % unit:RBM/dollar
crudeoilprice = 60;   % unit:dollar/STB
Po = crudeoilprice * RMBTODOLLAR / m3toSTB;     % price of oil (yuan/m3)
Pw = 30;       % price of water dealing (yuan/m3)
Pg = 20000;    % price of gas (yuan/10E4 m3)
Pig = 3000;    % price of gas injectioin (yuan/10E4 m3)
Dg = 15;       % total outlet gas dealing capacity (10E4 m3/d)
DL = 100;      % total outlet liquid dealing capacity (m3/d)
x0 = 1*ones(1,wellquantities);    % initial points of fmincon (10E4 m3/d)
A = ones(1,wellquantities); gas_provide = 15;   % the constraint of the sum of gas injection of each well (10E4 m3/d)

if Po <= 0 || Pg <= 0 || Pw <= 0 || Pig <= 0 || Dg <= 0 || DL <= 0 || gas_provide <= 0
    error('Po,Pg,Pw,Pig,Dg,DL,b must be positive')
elseif Po <= Pw
    error('Po must be larger than Pw')
end

%%
n = 8;                % set fitting power

if n <= 0
    error('n must be a positive')
end

%%
objectmodeltype = ['p' 'b' 'm'];  % 'p' denote production, 'b' denote benefit, 'm' denote mix

objectmodel = 1;     
if objectmodel ~= 1 && objectmodel ~= 2 && objectmodel ~= 3
    error('objectmodel must be 1,2 or 3')
end
% notice!!!  when w1 is close to 0 or 1, because the divisor(denominater)
% 'ytolrealmax_benefit' and 'ytolrealmax_production' are large in some extent, in the process of calculation
% it will lead to some extent of spreading error, so the minmum of object function
% 'w1 * f1 / ytolrealmax_production' and 'w2 * f2 / ytolrealmax_benefit' can't be exact with
% the minimum of object function 'f1' and 'f2'
if objectmodel == 3
    w1 = 0.8 ; w2 = 1 - w1;
    if w1 < 0 || w1 > 1
        error('w1 must be between [0,1]');
    end
end

%%
solvertype = ['f' 'g' 'i'];  % 'f' denote fmincon, 'g' denote globalsearch, 'i' denote interiorpointPF

solver = 2;

if solver ~= 1 && solver ~= 2  && solver ~= 3
    error('solver must be 1,2,3')
end
% if solver == 3
%     x0 = [2 2 1 1.5 1.5 0.5];
% r0 = 0.1; c = 0.1;
% % for i = 1:wellquantities
% %    syms gas_injection(i);
% % end
% % syms gas_injection1 gas_injection2 gas_injection3 gas_injection4 gas_injection5 gas_injection6
% % sym gas_injection;
% % variable = [gas_injection(1) gas_injection(2) gas_injection(3) gas_injection(4) gas_injection(5) gas_injection(6)];
% end

%%
algorithmoffmincontype = ['i' 'S' 'a'];

algorithmoffmincon = 1;     %  algorithmoffmincon must be 1,2 or 3

if algorithmoffmincon ~= 1 && algorithmoffmincon ~= 2 && algorithmoffmincon ~= 3
    error('algorithmoffmincon must be 1, 2 or 3')
end
% if algorithmoffmincon == 1, algorithm is interior-point;
% if algorithmoffmincon == 2, algorithm is SQP;
% if algorithmoffmincon == 3, algorithm is active-set.

%%
% graph must be 1 or 2 to show 'plots' plus 'PlotFcns of fmincon' of total wells...
% or 'PlotFcns of fmincon' of a certain well

graph = 1;

if graph ~= 1 && graph ~= 2
    error('graph must be 1 or 2')
end

if graph == 2
    acertainwell = 1;
    if acertainwell ~= 1 && acertainwell ~= 2 && acertainwell ~= 3 && acertainwell ~= 4 && ...
            acertainwell ~= 5 && acertainwell ~= 6 && acertainwell ~= 7
        error('acertainwell ~= 6 must be 1,2,3,4,5,6,7')
    end
end

%%
% define 'plotfunctions' to choose what plot of 'fmincon' will show
plotfunctions_fmincontype = ['x' 'f' 'v' 'c' 's' 'o'];
% 'x' denote the value of amount of gas injection, 'f' denote funccount,
% 'v' denote fval, 'c' denote constrviolation,
% 's' denote stepsize, 'o' denote firstorderopt

plotfunctions_fmincon = 1;

if plotfunctions_fmincon ~= 1 && plotfunctions_fmincon ~= 2 && plotfunctions_fmincon ~= 3 &&...
        plotfunctions_fmincon ~= 4 && plotfunctions_fmincon ~= 5 && plotfunctions_fmincon ~= 6
    error('plotfunctions must be 1,2,3,4,5,6')
end

plotfunctions_globalsearchtype = ['b' 'f'];   % 'b' denote bestf, 'f' denote funccount

plotfunctions_globalsearch = 1;

if plotfunctions_globalsearch ~= 1 && plotfunctions_globalsearch ~= 2
    error('plotfunctions must be 1,2')
end

languagetype = ['E' 'C'];    % 'E' denote Eglish, 'C' denote Chinese
language = 2;
if language ~= 1 && language ~= 2
    error('language must be 1,2')
end

injectiongastypename = ['g' 'n'];    % 1 == 'g' denote natural gas, 2 == 'n' denote nitrogen gas

injectiongastype = 1;

if injectiongastype ~= 1 && injectiongastype ~= 2
    error('injectiongastype must be 1,2')
end

calculationerror = 0.15;
if calculationerror <= 0
    error('calculationerror must > 0')
end

constraintype = 'inequation';        % constraintype is 'inequation' or 'equation'

%%
% above this section it is parameter definition and data input, below this
% it is code calculation








%%
XL = zeros(wellquantities,1); XR = zeros(wellquantities,1); coe = zeros(wellquantities,n+1);
for i = 1:wellquantities
    coe(i,:) = polyfit(data{i,1},data{i,2},n);
    XL(i) = min(data{i,1}); XR(i) = max(data{i,1});
end

%%
options = cell(length(objectmodeltype),wellquantities+1);

if algorithmoffmincon == 1
    algorithmoffmincon = 'interior-point';
elseif algorithmoffmincon == 2
    algorithmoffmincon = 'sqp';
elseif algorithmoffmincon == 3
    algorithmoffmincon = 'active-set';
end

%%
if plotfunctions_fmincon == 1
    plotfunctions_fmincon = @optimplotx;
elseif plotfunctions_fmincon == 2
    plotfunctions_fmincon = @optimplotfunccount;
elseif plotfunctions_fmincon == 3
    plotfunctions_fmincon = @optimplotfval;
elseif plotfunctions_fmincon == 4
    plotfunctions_fmincon = @optimplotconstrviolation;
elseif plotfunctions_fmincon == 5
    plotfunctions_fmincon = @optimplotstepsize ;
elseif plotfunctions_fmincon == 6
    plotfunctions_fmincon = @optimplotfirstorderopt ;
end

if plotfunctions_globalsearch == 1
    plotfunctions_globalsearch = @gsplotbestf;
elseif plotfunctions_globalsearch == 2
    plotfunctions_globalsearch = @gsplotfunccount;
end

%%

if solver == 1
    for i = 1:length(objectmodeltype)
        for j = 1:wellquantities
            options{i,j} = optimoptions(@fmincon,'Algorithm',algorithmoffmincon);
        end
    end
end

if solver == 2
    for i = 1:length(objectmodeltype)
        for j = 1:wellquantities
            options{i,j} = optimoptions(@fmincon, 'Algorithm',algorithmoffmincon,...
            'MaxFunEvals',10000, 'MaxIter',4000, 'TolCon',1e-12, 'TolFun',1e-12, 'TolX',1e-20);
            % 'DerivativeCheck','on', 'display','iter-detailed', 'Diagnostics','on', 'FinDiffType','central', 'FunValCheck','on',
        end
    end
end

if graph == 2    % to show 'PlotFcns of fmincon' of a certain well
    options{objectmodel,acertainwell} = optimoptions(@fmincon, 'Algorithm',algorithmoffmincon,...
        'display','iter-detailed', 'DerivativeCheck','on', 'Diagnostics','on',...
        'FinDiffType','central', 'FunValCheck','on', 'MaxFunEvals',10000, 'MaxIter',4000,...
        'TolCon',1e-12, 'TolFun',1e-12, 'TolX',1e-20,'PlotFcns',plotfunctions_fmincon);
end

%%
% problem = cell{3,wellquantities+1};
% output = cell{3,wellquantities+1};
% allmins = cell{3,wellquantities+1};
% the input arguments and ouput arguments of built-in function can't be predefined as cell array
% xmax = zeros(3,wellquantities);
% ymax = zeros(3,wellquantities);
% if ymin are assigned to zeros (3,7), it will lead to f_mix erroring

%% 
exitflag = zeros(length(objectmodeltype),wellquantities+1);
ymin = zeros(length(objectmodeltype),wellquantities);
output = cell(length(objectmodeltype),wellquantities+1);
allmins = cell(length(objectmodeltype),wellquantities+1);
%solve extremum of each well during gas injection area
%% 

if objectmodel == 1     % maximum production object function
    problem{objectmodel,1} = createOptimProblem('fmincon', 'objective',@f1_production, 'x0',ones(1,1),...
        'lb',XL(1), 'ub',XR(1), 'options',options{objectmodel,1});
    problem{objectmodel,2} = createOptimProblem('fmincon', 'objective',@f2_production, 'x0',ones(1,1),...
        'lb',XL(2), 'ub',XR(2), 'options',options{objectmodel,2});
    problem{objectmodel,3} = createOptimProblem('fmincon', 'objective',@f3_production, 'x0',ones(1,1),...
        'lb',XL(3),  'ub',XR(3), 'options',options{objectmodel,3});
    problem{objectmodel,4} = createOptimProblem('fmincon', 'objective',@f4_production, 'x0',ones(1,1),...
        'lb',XL(4), 'ub',XR(4), 'options',options{objectmodel,4});
    problem{objectmodel,5} = createOptimProblem('fmincon', 'objective',@f5_production, 'x0',ones(1,1),...
        'lb',XL(5), 'ub',XR(5), 'options',options{objectmodel,5});
    
    if solver == 1 || 3
%         [xmax(objectmodel,1),ymin(objectmodel,1),exitflag(objectmodel,1),output{objectmodel,1}]...
%             = fmincon(problem{objectmodel,1});
%         [xmax(objectmodel,2),ymin(objectmodel,2),exitflag(objectmodel,2),output{objectmodel,2}]...
%             = fmincon(problem{objectmodel,2});
%         [xmax(objectmodel,3),ymin(objectmodel,3),exitflag(objectmodel,3),output{objectmodel,3}]...
%             = fmincon(problem{objectmodel,3});
%         [xmax(objectmodel,4),ymin(objectmodel,4),exitflag(objectmodel,4),output{objectmodel,4}]...
%             = fmincon(problem{objectmodel,4});
%         [xmax(objectmodel,5),ymin(objectmodel,5),exitflag(objectmodel,5),output{objectmodel,5}]...
%             = fmincon(problem{objectmodel,5});
        for i = 1:wellquantities
                [xmax(objectmodel,i),ymin(objectmodel,i),exitflag(objectmodel,i),output{objectmodel,i}]...
            = fmincon(problem{objectmodel,i});
        end
    end
    
    if solver == 2
        gs = GlobalSearch;
%         [xmax(objectmodel,1),ymin(objectmodel,1),exitflag(objectmodel,1),output{objectmodel,1},...
%             allmins{objectmodel,1}] = run(gs,problem{objectmodel,1});
%         [xmax(objectmodel,2),ymin(objectmodel,2),exitflag(objectmodel,2),output{objectmodel,2},...
%             allmins{objectmodel,2}] = run(gs,problem{objectmodel,2});
%         [xmax(objectmodel,3),ymin(objectmodel,3),exitflag(objectmodel,3),output{objectmodel,3},...
%             allmins{objectmodel,3}] = run(gs,problem{objectmodel,3});
%         [xmax(objectmodel,4),ymin(objectmodel,4),exitflag(objectmodel,4),output{objectmodel,4},...
%             allmins{objectmodel,4}] = run(gs,problem{objectmodel,4});
%         [xmax(objectmodel,5),ymin(objectmodel,5),exitflag(objectmodel,5),output{objectmodel,5},...
%             allmins{objectmodel,5}] = run(gs,problem{objectmodel,5});
         for i = 1:wellquantities
        [xmax(objectmodel,i),ymin(objectmodel,i),exitflag(objectmodel,i),output{objectmodel,i},...
            allmins{objectmodel,i}] = run(gs,problem{objectmodel,i});
        end       
    end
end

%%
if objectmodel == 2   % maximum economic benefits object function
    
    problem{objectmodel,1} = createOptimProblem('fmincon', 'objective',@f1_benefit, 'x0',ones(1,1),...
        'lb',XL(1), 'ub',XR(1), 'options',options{objectmodel,1});
    problem{objectmodel,2} = createOptimProblem('fmincon', 'objective',@f2_benefit, 'x0',ones(1,1),...
        'lb',XL(2), 'ub',XR(2), 'options',options{objectmodel,2});
    problem{objectmodel,3} = createOptimProblem('fmincon', 'objective',@f3_benefit, 'x0',ones(1,1),...
        'lb',XL(3), 'ub',XR(3), 'options',options{objectmodel,3});
    problem{objectmodel,4} = createOptimProblem('fmincon', 'objective',@f4_benefit, 'x0',ones(1,1),...
        'lb',XL(4), 'ub',XR(4), 'options',options{objectmodel,4});
    problem{objectmodel,5} = createOptimProblem('fmincon', 'objective',@f5_benefit, 'x0',ones(1,1),...
        'lb',XL(5), 'ub',XR(5), 'options',options{objectmodel,5});
    
    if solver == 1 || 3
%         [xmax(objectmodel,1),ymin(objectmodel,1),exitflag(objectmodel,1),output{objectmodel,1}] = fmincon(problem{objectmodel,1});
%         [xmax(objectmodel,2),ymin(objectmodel,2),exitflag(objectmodel,2),output{objectmodel,2}] = fmincon(problem{objectmodel,2});
%         [xmax(objectmodel,3),ymin(objectmodel,3),exitflag(objectmodel,3),output{objectmodel,3}] = fmincon(problem{objectmodel,3});
%         [xmax(objectmodel,4),ymin(objectmodel,4),exitflag(objectmodel,4),output{objectmodel,4}] = fmincon(problem{objectmodel,4});
%         [xmax(objectmodel,5),ymin(objectmodel,5),exitflag(objectmodel,5),output{objectmodel,5}] = fmincon(problem{objectmodel,5});
        for i = 1:wellquantities
            [xmax(objectmodel,i),ymin(objectmodel,i),exitflag(objectmodel,i),output{objectmodel,i}]...
                = fmincon(problem{objectmodel,i});
        end
    end
    
    if solver == 2
        gs = GlobalSearch;
%         [xmax(objectmodel,1),ymin(objectmodel,1),exitflag(objectmodel,1),output{objectmodel,1},allmins{objectmodel,1}] = run(gs,problem{objectmodel,1});
%         [xmax(objectmodel,2),ymin(objectmodel,2),exitflag(objectmodel,2),output{objectmodel,2},allmins{objectmodel,2}] = run(gs,problem{objectmodel,2});
%         [xmax(objectmodel,3),ymin(objectmodel,3),exitflag(objectmodel,3),output{objectmodel,3},allmins{objectmodel,3}] = run(gs,problem{objectmodel,3});
%         [xmax(objectmodel,4),ymin(objectmodel,4),exitflag(objectmodel,4),output{objectmodel,4},allmins{objectmodel,4}] = run(gs,problem{objectmodel,4});
%         [xmax(objectmodel,5),ymin(objectmodel,5),exitflag(objectmodel,5),output{objectmodel,5},allmins{objectmodel,5}] = run(gs,problem{objectmodel,5});
        for i = 1:wellquantities
            [xmax(objectmodel,i),ymin(objectmodel,i),exitflag(objectmodel,i),output{objectmodel,i},allmins{objectmodel,i}]...
                = run(gs,problem{objectmodel,i});         
        end
    end
end

%%
if objectmodel == 3    % maximum mixed object function
    
    problem{objectmodel-2,1} = createOptimProblem('fmincon','objective',@f1_production,'x0',ones(1,1),...
        'lb',XL(1),'ub',XR(1),'options',options{objectmodel-2,1});
    problem{objectmodel-2,2} = createOptimProblem('fmincon','objective',@f2_production,'x0',ones(1,1),...
        'lb',XL(2),'ub',XR(2),'options',options{objectmodel-2,2});
    problem{objectmodel-2,3} = createOptimProblem('fmincon','objective',@f3_production,'x0',ones(1,1),...
        'lb',XL(3),'ub',XR(3),'options',options{objectmodel-2,3});
    problem{objectmodel-2,4} = createOptimProblem('fmincon','objective',@f4_production,'x0',ones(1,1),...
        'lb',XL(4),'ub',XR(4),'options',options{objectmodel-2,4});
    problem{objectmodel-2,5} = createOptimProblem('fmincon','objective',@f5_production,'x0',ones(1,1),...
        'lb',XL(5),'ub',XR(5),'options',options{objectmodel-2,5});
    
    problem{objectmodel-1,1} = createOptimProblem('fmincon','objective',@f1_benefit,'x0',ones(1,1),...
        'lb',XL(1),'ub',XR(1),'options',options{objectmodel-1,1});
    problem{objectmodel-1,2} = createOptimProblem('fmincon','objective',@f2_benefit,'x0',ones(1,1),...
        'lb',XL(2),'ub',XR(2),'options',options{objectmodel-1,2});
    problem{objectmodel-1,3} = createOptimProblem('fmincon', 'objective',@f3_benefit, 'x0',ones(1,1),...
        'lb',XL(3), 'ub',XR(3), 'options',options{objectmodel-1,3});
    problem{objectmodel-1,4} = createOptimProblem('fmincon','objective',@f4_benefit,'x0',ones(1,1),...
        'lb',XL(4),'ub',XR(4),'options',options{objectmodel-1,4});
    problem{objectmodel-1,5} = createOptimProblem('fmincon','objective',@f5_benefit,'x0',ones(1,1),...
        'lb',XL(5),'ub',XR(5),'options',options{objectmodel-1,5});
    
    problem{objectmodel,1} = createOptimProblem('fmincon','objective',@f1_mix,'x0',ones(1,1),...
        'lb',XL(1),'ub',XR(1),'options',options{objectmodel,1});
    problem{objectmodel,2} = createOptimProblem('fmincon','objective',@f2_mix,'x0',ones(1,1),...
        'lb',XL(2),'ub',XR(2),'options',options{objectmodel,2});
    problem{objectmodel,3} = createOptimProblem('fmincon','objective',@f3_mix,'x0',ones(1,1),...
        'lb',XL(3),'ub',XR(3),'options',options{objectmodel,3});
    problem{objectmodel,4} = createOptimProblem('fmincon','objective',@f4_mix,'x0',ones(1,1),...
        'lb',XL(4),'ub',XR(4),'options',options{objectmodel,4});
    problem{objectmodel,5} = createOptimProblem('fmincon','objective',@f5_mix,'x0',ones(1,1),...
        'lb',XL(5),'ub',XR(5),'options',options{objectmodel,5});
    
    if solver == 1 || 3
%         [xmax(objectmodel-2,1),ymin(objectmodel-2,1),exitflag(objectmodel-2,1),output{objectmodel-2,1}] = fmincon(problem{objectmodel-2,1});
%         [xmax(objectmodel-2,2),ymin(objectmodel-2,2),exitflag(objectmodel-2,2),output{objectmodel-2,2}] = fmincon(problem{objectmodel-2,2});
%         [xmax(objectmodel-2,3),ymin(objectmodel-2,3),exitflag(objectmodel-2,3),output{objectmodel-2,3}] = fmincon(problem{objectmodel-2,3});
%         [xmax(objectmodel-2,4),ymin(objectmodel-2,4),exitflag(objectmodel-2,4),output{objectmodel-2,4}] = fmincon(problem{objectmodel-2,4});
%         [xmax(objectmodel-2,5),ymin(objectmodel-2,5),exitflag(objectmodel-2,5),output{objectmodel-2,5}] = fmincon(problem{objectmodel-2,5});
%         
%         [xmax(objectmodel-1,1),ymin(objectmodel-1,1),exitflag(objectmodel-1,1),output{objectmodel-1,1}] = fmincon(problem{objectmodel-1,1});
%         [xmax(objectmodel-1,2),ymin(objectmodel-1,2),exitflag(objectmodel-1,2),output{objectmodel-1,2}] = fmincon(problem{objectmodel-1,2});
%         [xmax(objectmodel-1,3),ymin(objectmodel-1,3),exitflag(objectmodel-1,3),output{objectmodel-1,3}] = fmincon(problem{objectmodel-1,3});
%         [xmax(objectmodel-1,4),ymin(objectmodel-1,4),exitflag(objectmodel-1,4),output{objectmodel-1,4}] = fmincon(problem{objectmodel-1,4});
%         [xmax(objectmodel-1,5),ymin(objectmodel-1,5),exitflag(objectmodel-1,5),output{objectmodel-1,5}] = fmincon(problem{objectmodel-1,5});
%         
%         [xmax(objectmodel,1),ymin(objectmodel,1),exitflag(objectmodel,1),output{objectmodel,1}] = fmincon(problem{objectmodel,1});
%         [xmax(objectmodel,2),ymin(objectmodel,2),exitflag(objectmodel,2),output{objectmodel,2}] = fmincon(problem{objectmodel,2});
%         [xmax(objectmodel,3),ymin(objectmodel,3),exitflag(objectmodel,3),output{objectmodel,3}] = fmincon(problem{objectmodel,3});
%         [xmax(objectmodel,4),ymin(objectmodel,4),exitflag(objectmodel,4),output{objectmodel,4}] = fmincon(problem{objectmodel,4});
%         [xmax(objectmodel,5),ymin(objectmodel,5),exitflag(objectmodel,5),output{objectmodel,5}] = fmincon(problem{objectmodel,5});
        for i = 1:wellquantities
            [xmax(objectmodel-2,i),ymin(objectmodel-2,i),exitflag(objectmodel-2,i),output{objectmodel-2,i}] = fmincon(problem{objectmodel-2,i});
            [xmax(objectmodel-1,i),ymin(objectmodel-1,i),exitflag(objectmodel-1,i),output{objectmodel-1,i}] = fmincon(problem{objectmodel-1,i});
            [xmax(objectmodel,i),ymin(objectmodel,i),exitflag(objectmodel,i),output{objectmodel,i}] = fmincon(problem{objectmodel,i});     
        end
    end
    
    if solver == 2
        gs = GlobalSearch;
        for i = 1:wellquantities
        [xmax(objectmodel-2,i),ymin(objectmodel-2,i),exitflag(objectmodel-2,i),output{objectmodel-2,i}] = run(gs,problem{objectmodel-2,i});
        [xmax(objectmodel-1,i),ymin(objectmodel-1,i),exitflag(objectmodel-1,i),output{objectmodel-1,i}] = run(gs,problem{objectmodel-1,i});
        [xmax(objectmodel,i),ymin(objectmodel,i),exitflag(objectmodel,i),output{objectmodel,i}] = run(gs,problem{objectmodel,i});        
        end
%         [xmax(objectmodel-2,1),ymin(objectmodel-2,1),exitflag(objectmodel-2,1),output{objectmodel-2,1}] = run(gs,problem{objectmodel-2,1});
%         [xmax(objectmodel-2,2),ymin(objectmodel-2,2),exitflag(objectmodel-2,2),output{objectmodel-2,2}] = run(gs,problem{objectmodel-2,2});
%         [xmax(objectmodel-2,3),ymin(objectmodel-2,3),exitflag(objectmodel-2,3),output{objectmodel-2,3}] = run(gs,problem{objectmodel-2,3});
%         [xmax(objectmodel-2,4),ymin(objectmodel-2,4),exitflag(objectmodel-2,4),output{objectmodel-2,4}] = run(gs,problem{objectmodel-2,4});
%         [xmax(objectmodel-2,5),ymin(objectmodel-2,5),exitflag(objectmodel-2,5),output{objectmodel-2,5}] = run(gs,problem{objectmodel-2,5});
%         
%         [xmax(objectmodel-1,1),ymin(objectmodel-1,1),exitflag(objectmodel-1,1),output{objectmodel-1,1}] = run(gs,problem{objectmodel-1,1});
%         [xmax(objectmodel-1,2),ymin(objectmodel-1,2),exitflag(objectmodel-1,2),output{objectmodel-1,2}] = run(gs,problem{objectmodel-1,2});
%         [xmax(objectmodel-1,3),ymin(objectmodel-1,3),exitflag(objectmodel-1,3),output{objectmodel-1,3}] = run(gs,problem{objectmodel-1,3});
%         [xmax(objectmodel-1,4),ymin(objectmodel-1,4),exitflag(objectmodel-1,4),output{objectmodel-1,4}] = run(gs,problem{objectmodel-1,4});
%         [xmax(objectmodel-1,5),ymin(objectmodel-1,5),exitflag(objectmodel-1,5),output{objectmodel-1,5}] = run(gs,problem{objectmodel-1,5});
%         
%         [xmax(objectmodel,1),ymin(objectmodel,1),exitflag(objectmodel,1),output{objectmodel,1}] = run(gs,problem{objectmodel,1});
%         [xmax(objectmodel,2),ymin(objectmodel,2),exitflag(objectmodel,2),output{objectmodel,2}] = run(gs,problem{objectmodel,2});
%         [xmax(objectmodel,3),ymin(objectmodel,3),exitflag(objectmodel,3),output{objectmodel,3}] = run(gs,problem{objectmodel,3});
%         [xmax(objectmodel,4),ymin(objectmodel,4),exitflag(objectmodel,4),output{objectmodel,4}] = run(gs,problem{objectmodel,4});
%         [xmax(objectmodel,5),ymin(objectmodel,5),exitflag(objectmodel,5),output{objectmodel,5}] = run(gs,problem{objectmodel,5});
    end
end

%%
ymax = zeros(length(objectmodeltype),wellquantities);
for i = 1:length(objectmodeltype)
    for j = 1:wellquantities
        ymax(i,j) = -ymin(i,j);
    end
end

%%
% grid on;   establish noly one blank graph with grid once in a file
if graph == 1  % && objectmodel == 1
    % figure;
    hold on;
    xmaxofdata = zeros(1,wellquantities); indices = zeros(1,wellquantities);
    for i = 1:wellquantities
        indices(i) = find( data{i,2} == max(data{i,2}) );
        xmaxofdata(i) = data{i,1}(indices(i));
        if abs( ymax(1,i) - max(data{i,2}) ) / max(data{i,2}) < calculationerror &&...
                abs( xmax(1,i) - xmaxofdata(i) ) / xmaxofdata(i) < calculationerror + 0.3
            % telling if calculation optimization value is trustable
            disp 'single well maximum production(m3/d) are';
            disp(ymax(1,i));
            disp 'single well minimum injection(10E4 m3 / d) and maximum injection(10E4 m3 / d) are';
            disp(XL(i)); disp(xmax(1,i));
            plot(data{i,1},data{i,2});
        end
    end
    if language == 1
        title('injection & production curve');
        xlabel('gas injection/(10^4 m^3/d)');
        ylabel('liquid production/(m^3/d)');
    elseif language == 2
        title('注气产液量实际曲线图');
        xlabel('注气量/(10^4 m^3/天)');
        ylabel('产液量/(m^3/天)');
    end
    hold off;
end

% elseif objectmodel == 2      % maximum economic benefits object function
% for i = 1:wellquantities
%         disp 'single well maximum economic benefits (dollar/d) are';
%         disp(ymax_benefit(i));
%         disp 'single well minimum injection(10E4 m3 / d) and maximum injection(10E4 m3 / d) are';
%         disp(XL(i)); disp(xmax_benefit(i));
%         plot(data{i,1},data{i,2});
% end
% title('injection & benefits curve');
% xlabel('gas injection (unit:10E4 m3 / d)'); ylabel('economic benefits (unit:dollar/d)');
%
% elseif objectmodel == 3      % maximum mixed object function
% for i = 1:wellquantities
%         disp 'single well maximum the mixed object function (dimensionless) are';
%         disp(ymax_mix(i));
%         disp 'single well minimum injection(10E4 m3 / d) and maximum injection(10E4 m3 / d) are';
%         disp(XL(i)); disp(xmax_mix(i));
%         plot(data{i,1},data{i,2});
% end
% title('injection & dimensionless curve');
% xlabel('gas injection (unit:10E4 m3 / d)'); ylabel('the mixed object function (unit:none)');
% end

%%
for i = 1:length(objectmodeltype)
    if solver == 1
        
        if graph == 1
            options{i,wellquantities+1} = optimoptions(@fmincon, 'Algorithm',algorithmoffmincon,...
                'Diagnostics','on', 'DerivativeCheck','on', 'display','iter-detailed',...
                'FinDiffType','central', 'FunValCheck','on', 'MaxFunEvals',10000, 'MaxIter',3000,...
                'TolCon',1e-12, 'TolFun',1e-12, 'TolX',1e-20, 'PlotFcns',plotfunctions_fmincon);
        end
        
        if graph == 2     % not to show 'PlotFcns' of 'fmincon'
            options{i,wellquantities+1} = optimoptions(@fmincon, 'Algorithm',algorithmoffmincon,...
                'Diagnostics','on', 'DerivativeCheck','on', 'display','iter-detailed',...
                'FinDiffType','central', 'FunValCheck','on', 'MaxFunEvals',10000, 'MaxIter',3000,...
                'TolCon',1e-12, 'TolFun',1e-12, 'TolX',1e-20);
        end
        
    end
    if solver == 2
        
        if graph == 1
            options{i,wellquantities+1} = optimoptions(@fmincon,'Algorithm',algorithmoffmincon,...
                'DerivativeCheck','on',...
                'FinDiffType','central','FunValCheck','on','MaxFunEvals',10000,'MaxIter',3000,...
                'TolCon',1e-12,'TolFun',1e-12,'TolX',1e-20,'PlotFcns',plotfunctions_fmincon);
        end
        
        if graph == 2     % not to show 'PlotFcns' of 'fmincon'
            options{i,wellquantities+1} = optimoptions(@fmincon,'Algorithm',algorithmoffmincon,...
                'DerivativeCheck','on',...
                'FinDiffType','central','FunValCheck','on','MaxFunEvals',10000,'MaxIter',3000,...
                'TolCon',1e-12,'TolFun',1e-12,'TolX',1e-20);
        end
        
    end
end

%%
% xmax_total = cell{1,3};
% ymin_total = cell{1,3};
% the input arguments and output
% arguments of built-in functions are not allowed to be predefined.
%%
gs = GlobalSearch('Display','iter', 'PlotFcns',plotfunctions_globalsearch);     % @gsplotfunccount   @gsplotbestf
%% 

if objectmodel == 1    
    LB = [XL(1) XL(2) XL(3) XL(4) XL(5)];
    UB_production = [( xmax(objectmodel,1)+XR(1) )/2 ( xmax(objectmodel,2)+XR(2) )/2 ...
        ( xmax(objectmodel,3)+XR(3) )/2 ( xmax(objectmodel,4)+XR(4) )/2 ( xmax(objectmodel,5)+XR(5) )/2];
%     UB_production = [xmax(objectmodel,1) xmax(objectmodel,2) xmax(objectmodel,3),...
%         xmax(objectmodel,4) xmax(objectmodel,5)];      % constraint of gas injection of each well
    
    if solver == 1
        switch constraintype
            case 'inequation'
                % use fmincon to solve optimizaiton distribution
                [xmax_total{objectmodel},ymin_total{objectmodel},exitflag(objectmodel,wellquantities+1),...
                    output{objectmodel,wellquantities+1},lambda,grad,hessian] = fmincon(@objfunc_production,...
                    x0, A, gas_provide,[], [], LB, UB_production,@noncons_inequation,options{objectmodel,wellquantities+1});
            case 'equation'
                [xmax_total{objectmodel},ymin_total{objectmodel},exitflag(objectmodel,wellquantities+1),...
                    output{objectmodel,wellquantities+1},lambda,grad,hessian] = fmincon(@objfunc_production,...
                    x0, [], [],A, gas_provide, LB, UB_production,@noncons_equation,options{objectmodel,wellquantities+1});
        end
    end
    
    if solver == 2
        switch constraintype
            case 'inequation'
                % use GlobalSearch to solve optimizaiton distribution
                problem{objectmodel,wellquantities+1} = createOptimProblem('fmincon','objective',...
                    @objfunc_production,'x0',x0,'lb',LB,'ub',UB_production,'Aineq',A,...
                    'bineq',gas_provide,'nonlcon',@noncons_inequation,'options',options{objectmodel,wellquantities+1});
            case 'equation'
                problem{objectmodel,wellquantities+1} = createOptimProblem('fmincon','objective',...
                    @objfunc_production,'x0',x0,'lb',LB,'ub',UB_production,'Aeq',A,...
                    'beq',gas_provide,'nonlcon',@noncons_equation,'options',options{objectmodel,wellquantities+1});
        end
        [xmax_total{objectmodel},ymin_total{objectmodel},exitflag(objectmodel,wellquantities+1),...
            output{objectmodel,wellquantities+1},allmins{objectmodel,wellquantities+1}]...
            = run(gs,problem{objectmodel,wellquantities+1});
    end
    
    %%
    %     if solver == 3
    %         % use interior point punishment function to solve optimization
    %         % distributioin
    % %         sym t(1) t(2) t(3) t(4) t(5) t(6);
    % %        syms gasprovideconstraint(t(1),t(2),t(3),t(4),t(5),t(6))
    % %       gasprovideconstraint(t(1),t(2),t(3),t(4),t(5),t(6)) = gas_provide - t(1)-t(2)-t(3)-t(4)-t(5)-t(6);
    % % syms t;
    % [gas_deal,liquid_deal,gasprovideconstraint,xlowbound,xupbound] = inequality_constraint(t);
    %          inequality = [gas_deal,liquid_deal,gasprovideconstraint, xlowbound, xupbound];     % inequality is a row vector
    % [xmax_total{objectmodel},ymin_total{objectmodel}] = interiorpointPF(objfunc_production(t),x0,inequality,r0,c,[t(1),t(2),t(3),t(4),t(5),t(6)]);
    %     end
    %     max gas injection rate at max liquid production rate and unconstraint max liquid production...
    %     rate of block wells
    %%
    if solver == 1 || solver == 2
        disp(exitflag(objectmodel,wellquantities+1)); disp(output{objectmodel,wellquantities+1});
    end
    if solver == 1
        disp(lambda); disp(grad); disp(hessian);
    end
    
    Ginjmax = 0; Lprodmax = 0;
    for i = 1:wellquantities
        Ginjmax = Ginjmax + xmax(objectmodel,i); Lprodmax = Lprodmax + ymax(objectmodel,i);
    end
    
    disp('the constraint of total gas provide ability of all wells(10E4 m3 / d)');
    disp(gas_provide);
    disp('unconstraint maxmum gas injction rate(10E4 m3 / d)');
    disp(Ginjmax);
    toptsum = sum(xmax_total{objectmodel});
    disp('the sum of actual gas injction rate of all wells with constraint(10E4 m3 / d)');
    disp(toptsum);
    disp('actual gas injction rate of each well with constraint(10E4 m3 / d)');
    disp(xmax_total{objectmodel});
    
    y(objectmodel,1) = f1_production( xmax_total{objectmodel}(1) );
    y(objectmodel,2) = f2_production( xmax_total{objectmodel}(2) );
    y(objectmodel,3) = f3_production( xmax_total{objectmodel}(3) );
    y(objectmodel,4) = f4_production( xmax_total{objectmodel}(4) );
    y(objectmodel,5) = f5_production( xmax_total{objectmodel}(5) );
    disp('the constraint of total liquid deal ability of all wells(m3/d)');
    disp(DL);
    disp('unconstraint maxmum liquid production rate(m3/d)');
    disp(Lprodmax);
    Y = [-y(objectmodel,1) -y(objectmodel,2) -y(objectmodel,3) -y(objectmodel,4) -y(objectmodel,5)];
    ysum = sum(Y);
    if abs(ysum + ymin_total{objectmodel})/ysum > 10^(-10)   % ysum is positive and ymin_total is negative so they can be added.
        error('the difference between total function and sum of all functions are comparatively large');
    end
    disp('the sum of actual maxmum liquid production rate with constraint(m3/d)');
    disp(ysum);
    disp('actual liquid production rate of each well with constraint(m3/d)');
    disp(Y);
    
    disp('the constraint of total gas deal ability of all wells(10E4 m3 / d)');
    disp(Dg);
    
    unconstraintnetproducegas = 0;
    for i = 1:wellquantities
        unconstraintnetproducegas = unconstraintnetproducegas + RpgL(i)*ymax(objectmodel,i)/10000;   % should be "+" because "ymax(objectmodel,i)" is positive
    end
    unconstraintgasneededtobedealt = Ginjmax + unconstraintnetproducegas;
    if unconstraintgasneededtobedealt < Ginjmax
        error('unconstraint total gas production rate is smaller than total gas injectioin rate');
    end
    disp('under unconstraint total gas rate needed to be dealt with(10E4 m3 / d)');
    disp(unconstraintgasneededtobedealt);
    
%     if injectiongastype == 1
%         if producegas > Ginjmax
%             actualgasneededtobedealt = producegas;
%         else
%             actualgasneededtobedealt = Ginjmax;
%         end
%     elseif  injectiongastype == 2
%         actualgasneededtobedealt = producegas + Ginjmax;
%     end
    
    producegas = 0;
    for i = 1:wellquantities
        producegas = producegas - ( RpgL(i)*y(objectmodel,i)...
            - xmax_total{objectmodel}(i)*10000 ) / 10000;   % should be "-" because "y(objectmodel,i)" is negative    
    end
    if producegas < toptsum
        error('actual total gas production rate is smaller than catual total gas injectioin rate');
    else
        actualgasneededtobedealt = producegas;
    end
    disp('with constraint total gas rate needed to be dealt with(10E4 m3 / d)');
    disp(actualgasneededtobedealt);
    
%     if injectiongastype == 1
%         if producegas > toptsum
%             actualgasneededtobedealt = producegas;
%         else
%             actualgasneededtobedealt = toptsum;
%         end
%     elseif  injectiongastype == 2
%         actualgasneededtobedealt = producegas + toptsum;
%     end
    
    if solver == 2
        disp(allmins{objectmodel,wellquantities+1});
    end
end

%%
if objectmodel == 2
    % use fmincon to solve optimizaiton distribution
    LB = [XL(1) XL(2) XL(3) XL(4) XL(5)];
    UB_benefit = [(xmax(objectmodel,1)+XR(1))/2 (xmax(objectmodel,2)+XR(2))/2 (xmax(objectmodel,3)+XR(3))/2,...
         (xmax(objectmodel,4)+XR(4))/2 (xmax(objectmodel,5)+XR(5))/2]; 
%     UB_benefit = [XR(1) XR(2) XR(3) XR(4) XR(5)];
%     UB_benefit = [xmax(objectmodel,1) xmax(objectmodel,objectmodel) xmax(objectmodel,3),...
%         xmax(objectmodel,4) xmax(objectmodel,5)];    % constraint of gas injection of each well
    
    if solver == 1
        switch constraintype
            case 'inequation'
                [xmax_total{objectmodel},ymin_total{objectmodel},exitflag(objectmodel,wellquantities+1),...
                    output{objectmodel,wellquantities+1},lambda,grad,hessian] = fmincon(@objfunc_benefit,...
                    x0, A,gas_provide,[], [], LB, UB_benefit,@noncons_inequation,options{objectmodel,wellquantities+1});
            case 'equation'
                [xmax_total{objectmodel},ymin_total{objectmodel},exitflag(objectmodel,wellquantities+1),...
                    output{objectmodel,wellquantities+1},lambda,grad,hessian] = fmincon(@objfunc_benefit,...
                    x0,[], [], A,gas_provide, LB, UB_benefit,@noncons_equation,options{objectmodel,wellquantities+1});
        end
    end
    
    if solver == 2
        switch constraintype
            case 'inequation'
                problem{objectmodel,wellquantities+1} = createOptimProblem('fmincon','objective',@objfunc_benefit,...
                    'x0',x0,'lb',LB,'ub',UB_benefit,'Aineq',A,'bineq',gas_provide,'nonlcon',@noncons_inequation,...
                    'options',options{objectmodel,wellquantities+1});
            case 'equation'
                problem{objectmodel,wellquantities+1} = createOptimProblem('fmincon','objective',@objfunc_benefit,...
                    'x0',x0,'lb',LB,'ub',UB_benefit,'Aeq',A,'beq',gas_provide,'nonlcon',@noncons_equation,...
                    'options',options{objectmodel,wellquantities+1});
        end
        [xmax_total{objectmodel},ymin_total{objectmodel},exitflag(objectmodel,wellquantities+1),...
            output{objectmodel,wellquantities+1},allmins{objectmodel,wellquantities+1}]...
            = run(gs,problem{objectmodel,wellquantities+1});
    end
    
    %%
    disp(exitflag(objectmodel,wellquantities+1)); disp(output{objectmodel,wellquantities+1});
    if solver == 1
        disp(lambda); disp(grad); disp(hessian);
    end
    
    % unconstraint max gas injection rate of all wells at max economic benefit and unconstraint max economic benefit of all wells
    Ginjmax = 0; Lprodmax = 0;
    for i = 1:wellquantities
        Ginjmax = Ginjmax + xmax(objectmodel,i); Lprodmax = Lprodmax + ymax(objectmodel,i);
    end
    
    disp('the constraint of total gas provide ability of all wells(10E4 m3 / d)');
    disp(gas_provide);
    disp('the sum of gas injection rate of block wells under unconstraint(10E4 m3 / d)');
    disp(Ginjmax);
    toptsum = sum(xmax_total{objectmodel});
    disp('the sum of actual gas injction rate of all wells with constraint(10E4 m3 / d)');
    disp(toptsum);
    disp('actual gas injction rate of each well under constraints(10E4 m3 / d)');
    disp(xmax_total{objectmodel});
    
    y_unconstraint(objectmodel,1) = f1_production( xmax(objectmodel,1) );    % y_unconstraint(objectmodel,i) are negative
    y_unconstraint(objectmodel,2) = f2_production( xmax(objectmodel,2) );
    y_unconstraint(objectmodel,3) = f3_production( xmax(objectmodel,3) );
    y_unconstraint(objectmodel,4) = f4_production( xmax(objectmodel,4) );
    y_unconstraint(objectmodel,5) = f5_production( xmax(objectmodel,5) );
    Y_unconstraint = [-y_unconstraint(objectmodel,1) -y_unconstraint(objectmodel,2) -y_unconstraint(objectmodel,3)...
         -y_unconstraint(objectmodel,4) -y_unconstraint(objectmodel,5)];
    ytotal_unconstraint = sum(Y_unconstraint); 
    
    y(objectmodel,1) = f1_production( xmax_total{objectmodel}(1) );
    y(objectmodel,2) = f2_production( xmax_total{objectmodel}(2) );
    y(objectmodel,3) = f3_production( xmax_total{objectmodel}(3) );
    y(objectmodel,4) = f4_production( xmax_total{objectmodel}(4) );
    y(objectmodel,5) = f5_production( xmax_total{objectmodel}(5) );
    Y = [-y(objectmodel,1) -y(objectmodel,2) -y(objectmodel,3) -y(objectmodel,4) -y(objectmodel,5)];
    ysum = sum(Y);
    
    disp('the constraint of total liquid deal ability of all wells(m3/d)');
    disp(DL);   
    disp('unconstraint maxmum liquid production rate(m3/d)');
    disp( ytotal_unconstraint );
    disp('the sum of actual maxmum liquid production rate with constraint(m3/d)');
    disp(ysum);
    disp('actual liquid production rate of each well with constraint(m3/d)');
    disp(Y);
    
    disp('the constraint of total gas deal ability of all wells(10E4 m3 / d)');
    disp(Dg);
    producegas = 0;
    for i = 1:wellquantities
        producegas = producegas - ( RpgL(i)*y_unconstraint(objectmodel,i)...
            - xmax(objectmodel,i)*10000 ) / 10000;    %  between 'producegas' and ...
        ...'RpgL(i)*y_unconstraint(objectmodel,1)/10000' must be '-', because 'y_unconstraint' is negative.
    end
    if producegas < Ginjmax
        error('unconstraint total gas production rate is smaller than unconstraint total gas injention rate');
    else
        actualgasneededtobedealt = producegas;
    end
%     if injectiongastype == 1
%         if producegas > Ginjmax
%             actualgasneededtobedealt = producegas;
%         else
%             actualgasneededtobedealt = Ginjmax;
%         end
%     elseif  injectiongastype == 2
%         actualgasneededtobedealt = producegas + Ginjmax;
%     end
    disp('under unconstraint total gas rate needed to be dealt with(10E4 m3 / d)');
    disp(actualgasneededtobedealt);
    
    producegas = 0;
    for i = 1:wellquantities
        producegas = producegas - ( RpgL(i)*y(objectmodel,i)...
             - xmax_total{objectmodel}(i)*10000 ) / 10000;   %  between 'producegas' and ...
        ...'RpgL(i)*y(objectmodel,i)/10000' must be '-', because 'y' is negative.
    end
    if producegas < toptsum
        error('actual total gas production rate is smaller than actual total gas injention rate');
    else
        actualgasneededtobedealt = producegas;
    end

%     if injectiongastype == 1
%         if producegas > toptsum
%             actualgasneededtobedealt = producegas;
%         else
%             actualgasneededtobedealt = toptsum;
%         end
%     elseif  injectiongastype == 2
%         actualgasneededtobedealt = producegas + toptsum;
%     end
    disp('with constraint total gas rate needed to be dealt with(10E4 m3 / d)');
    disp(actualgasneededtobedealt);
    
    disp('unconstraint maxmum economic benefits(dollar/d)');
    disp(Lprodmax);
    ytolrealmax_benefit = -ymin_total{objectmodel};
    disp('actual maxmum economic benefit with constraint(dollar/d)');
    disp(ytolrealmax_benefit);
    
    if solver == 2
        disp(allmins{objectmodel,wellquantities+1});
    end
end

%%

if objectmodel == 3
    % use fmincon to solve optimizaiton distribution
    % the constraints(lb,ub) of gas injection of each well
    LB = [XL(1) XL(2) XL(3) XL(4) XL(5)];
    UB_production = [(XR(1)+xmax(objectmodel-2,1))/2 (XR(2)+xmax(objectmodel-2,2))/2 (XR(3)+xmax(objectmodel-2,3))/2 ...
        (XR(4)+xmax(objectmodel-2,4))/2 (XR(5)+xmax(objectmodel-2,5))/2];
    UB_benefit = [(xmax(objectmodel-1,1)+xmax(objectmodel-2,1))/2 (xmax(objectmodel-1,2)+xmax(objectmodel-2,2))/2 ...
        (xmax(objectmodel-1,3)+xmax(objectmodel-2,3))/2 (xmax(objectmodel-1,4)+xmax(objectmodel-2,4))/2 (xmax(objectmodel-1,5)+xmax(objectmodel-2,5))/2];
    UB_mix = [(xmax(objectmodel,1)+xmax(objectmodel-2,1))/2 (xmax(objectmodel,2)+xmax(objectmodel-2,2))/2 ...
        (xmax(objectmodel,3)+xmax(objectmodel-2,3))/2 (xmax(objectmodel,4)+xmax(objectmodel-2,4))/2 (xmax(objectmodel,5)+xmax(objectmodel-2,5))/2];
%     UB_production = [xmax(objectmodel-2,1) xmax(objectmodel-2,2) xmax(objectmodel-2,3),...
%         xmax(objectmodel-2,4) xmax(objectmodel-2,5)];
%     UB_benefit = [xmax(objectmodel-1,1) xmax(objectmodel-1,2) xmax(objectmodel-1,3),...
%         xmax(objectmodel-1,4) xmax(objectmodel-1,5)];
%     UB_mix = [xmax(objectmodel,1) xmax(objectmodel,2) xmax(objectmodel,3) xmax(objectmodel,4),...
%         xmax(objectmodel,5)];
    
    if solver == 1
        switch constraintype
            case 'inequation'
                
                [xmax_total{objectmodel-2},ymin_total{objectmodel-2},exitflag(objectmodel-2,wellquantities+1),...
                    output{objectmodel-2,wellquantities+1}] = fmincon(@objfunc_production,x0,...
                    A,gas_provide,[], [], LB, UB_production,@noncons_inequation,options{objectmodel-2,wellquantities+1});
                ytolrealmax_production = -ymin_total{objectmodel-2};
                
                [xmax_total{objectmodel-1},ymin_total{objectmodel-1},exitflag(objectmodel-1,wellquantities+1),...
                    output{objectmodel-1,wellquantities+1}] = fmincon(@objfunc_benefit,x0,...
                    A,gas_provide,[], [], LB, UB_benefit,@noncons_inequation,options{objectmodel-1,wellquantities+1});
                ytolrealmax_benefit = -ymin_total{objectmodel-1};
                
                [xmax_total{objectmodel},ymin_total{objectmodel},exitflag(objectmodel,wellquantities+1),...
                    output{objectmodel,wellquantities+1},lambda,grad,hessian] = fmincon(@objfunc_mix,x0,...
                    A,gas_provide,[], [], LB, UB_mix,@noncons_inequation,options{objectmodel,wellquantities+1});
            case 'equation'
                [xmax_total{objectmodel-2},ymin_total{objectmodel-2},exitflag(objectmodel-2,wellquantities+1),...
                    output{objectmodel-2,wellquantities+1}] = fmincon(@objfunc_production,x0,...
                    [], [],A,gas_provide, LB, UB_production,@noncons_equation,options{objectmodel-2,wellquantities+1});
                ytolrealmax_production = -ymin_total{objectmodel-2};
                
                [xmax_total{objectmodel-1},ymin_total{objectmodel-1},exitflag(objectmodel-1,wellquantities+1),...
                    output{objectmodel-1,wellquantities+1}] = fmincon(@objfunc_benefit,x0,...
                    [], [],A,gas_provide, LB, UB_benefit,@noncons_equation,options{objectmodel-1,wellquantities+1});
                ytolrealmax_benefit = -ymin_total{objectmodel-1};
                
                [xmax_total{objectmodel},ymin_total{objectmodel},exitflag(objectmodel,wellquantities+1),...
                    output{objectmodel,wellquantities+1},lambda,grad,hessian] = fmincon(@objfunc_mix,x0,...
                    [], [],A,gas_provide, LB, UB_mix,@noncons_equation,options{objectmodel,wellquantities+1});
                
        end
    end
    
    if solver == 2
        switch constraintype
            case 'inequation'
                
                problem{objectmodel-2,wellquantities+1} = createOptimProblem('fmincon',...
                    'objective',@objfunc_production,'x0',x0,'lb',LB,'ub',UB_production,'Aineq',A,...
                    'bineq',gas_provide,'nonlcon',@noncons_inequation,'options',options{objectmodel-2,wellquantities+1});
                problem{objectmodel-1,wellquantities+1} = createOptimProblem('fmincon',...
                    'objective',@objfunc_benefit,'x0',x0,'lb',LB,'ub',UB_benefit,'Aineq',A,...
                    'bineq',gas_provide,'nonlcon',@noncons_inequation,'options',options{objectmodel-1,wellquantities+1});
                problem{objectmodel,wellquantities+1} = createOptimProblem('fmincon',...
                    'objective',@objfunc_mix,'x0',x0,'lb',LB,'ub',UB_mix,'Aineq',A,'bineq',gas_provide,...
                    'nonlcon',@noncons_inequation,'options',options{objectmodel,wellquantities+1});
            case 'equation'
                problem{objectmodel-2,wellquantities+1} = createOptimProblem('fmincon',...
                    'objective',@objfunc_production,'x0',x0,'lb',LB,'ub',UB_production,'Aeq',A,...
                    'beq',gas_provide,'nonlcon',@noncons_equation,'options',options{objectmodel-2,wellquantities+1});
                problem{objectmodel-1,wellquantities+1} = createOptimProblem('fmincon',...
                    'objective',@objfunc_benefit,'x0',x0,'lb',LB,'ub',UB_benefit,'Aeq',A,...
                    'beq',gas_provide,'nonlcon',@noncons_equation,'options',options{objectmodel-1,wellquantities+1});
                problem{objectmodel,wellquantities+1} = createOptimProblem('fmincon',...
                    'objective',@objfunc_mix,'x0',x0,'lb',LB,'ub',UB_mix,'Aeq',A,'beq',gas_provide,...
                    'nonlcon',@noncons_equation,'options',options{objectmodel,wellquantities+1});
                
        end
        [xmax_total{objectmodel-2},ymin_total{objectmodel-2},...
            exitflag(objectmodel-2,wellquantities+objectmodel-2),...
            output{objectmodel-2,wellquantities+1},allmins{objectmodel-2,wellquantities+1}]...
            = run(gs,problem{objectmodel-2,wellquantities+1});
        ytolrealmax_production = -ymin_total{objectmodel-2};
        
        [xmax_total{objectmodel-1},ymin_total{objectmodel-1},exitflag(objectmodel-1,wellquantities+1),...
            output{objectmodel-1,wellquantities+1},allmins{objectmodel-1,wellquantities+1}]...
            = run(gs,problem{objectmodel-1,wellquantities+1});
        ytolrealmax_benefit = -ymin_total{objectmodel-1};
        
        [xmax_total{objectmodel},ymin_total{objectmodel},exitflag(objectmodel,wellquantities+1),...
            output{objectmodel,wellquantities+1},allmins{objectmodel,wellquantities+1}]...
            = run(gs,problem{objectmodel,wellquantities+1});
    end
    
    %%
    disp(exitflag(objectmodel,wellquantities+1)); disp(output{objectmodel,wellquantities+1});
    if solver == 1
        disp(lambda); disp(grad); disp(hessian);
    end
    
    % max gas injection rate at max economic benefit and unconstraint max economic benefit of block wells
    Ginjmax = 0; Lprodmax = 0;
    for i = 1:wellquantities
        Ginjmax = Ginjmax + xmax(objectmodel,i); Lprodmax = Lprodmax + ymax(objectmodel,i);
    end
    
    disp('the constraint of total gas provide ability of all wells(10E4 m3 / d)');
    disp(gas_provide);
    disp('the sum of gas injection rate of block wells under unconstraint(10E4 m3 / d)');
    disp(Ginjmax);
    toptsum = sum(xmax_total{objectmodel});
    disp('the sum of actual gas injction rate of all wells with constraint(10E4 m3 / d)');
    disp(toptsum);
    disp('actual gas injction rate of each well with constraint(10E4 m3 / d)');
    disp(xmax_total{objectmodel});
    
    y_unconstraint(objectmodel,1) = f1_production( xmax(objectmodel,1) );
    y_unconstraint(objectmodel,2) = f2_production( xmax(objectmodel,2) );
    y_unconstraint(objectmodel,3) = f3_production( xmax(objectmodel,3) );
    y_unconstraint(objectmodel,4) = f4_production( xmax(objectmodel,4) );
    y_unconstraint(objectmodel,5) = f5_production( xmax(objectmodel,5) );
    Y_unconstraint = [-y_unconstraint(objectmodel,1) -y_unconstraint(objectmodel,2) -y_unconstraint(objectmodel,3)...
         -y_unconstraint(objectmodel,4) -y_unconstraint(objectmodel,5)];
    ytotal_unconstraint = sum(Y_unconstraint); 

    y(objectmodel,1) = f1_production( xmax_total{objectmodel}(1) );
    y(objectmodel,2) = f2_production( xmax_total{objectmodel}(2) );
    y(objectmodel,3) = f3_production( xmax_total{objectmodel}(3) );
    y(objectmodel,4) = f4_production( xmax_total{objectmodel}(4) );
    y(objectmodel,5) = f5_production( xmax_total{objectmodel}(5) );
    Y = [-y(objectmodel,1) -y(objectmodel,2) -y(objectmodel,3) -y(objectmodel,4) -y(objectmodel,5)];
    ysum = sum(Y);
    
    disp('the constraint of total liquid deal ability of all wells(m3/d)');
    disp(DL);
    disp('unconstraint maxmum liquid production rate(m3/d)');
    disp(ytotal_unconstraint);
    disp('the sum of actual maxmum liquid production rate with constraint(m3/d)');
    disp(ysum);
    disp('actual liquid production rate of each well with constraint(m3/d)');
    disp(Y);
    
    disp('the constraint of total gas deal ability of all wells(10E4 m3 / d)');
    disp(Dg);
    producegas = 0;
    for i = 1:wellquantities
        producegas = producegas - ( RpgL(i)*y_unconstraint(objectmodel,i)...
            - xmax(objectmodel,i)*10000 ) / 10000;   % "y_unconstraint" is negative
    end
    if producegas < Ginjmax
        error('uncontraint total gas production rate cannot be smaller than uncontraint total gas injection rate');
    else
        actualgasneededtobedealt = producegas;
    end
%     if injectiongastype == 1
%         if producegas > Ginjmax
%             actualgasneededtobedealt = producegas;
%         else
%             actualgasneededtobedealt = Ginjmax;
%         end
%     elseif  injectiongastype == 2
%         actualgasneededtobedealt = producegas + Ginjmax;
%     end
    disp('under unconstraint total gas rate needed to be dealt with(10E4 m3 / d)');
    disp(actualgasneededtobedealt);
    
    producegas = 0;
    for i = 1:wellquantities
        producegas = producegas - ( RpgL(i)*y(objectmodel,i)...
            - xmax_total{objectmodel}(i)*10000 ) / 10000;   %  "y" is negative
    end
    if producegas < toptsum
        error('actual total gas production rate cannot be smaller than actual total gas injection rate');
    else
        actualgasneededtobedealt = producegas;
    end
%     if injectiongastype == 1
%         if producegas > toptsum
%             actualgasneededtobedealt = producegas;
%         else
%             actualgasneededtobedealt = toptsum;
%         end
%     elseif  injectiongastype == 2
%         actualgasneededtobedealt = producegas + toptsum;
%     end
    disp('with constraint total gas rate needed to be dealt with(10E4 m3 / d)');
    disp(actualgasneededtobedealt);
    
    disp('the value of unconstraint maxmum mixed object function(dimensionless)');
    disp(Lprodmax);
    ytolrealmax_mix = -ymin_total{objectmodel};
    disp('with constraint the value of maxmum mixed object function with constraint(dimensionless)');
    disp(ytolrealmax_mix);
    disp('w1 and w2 are:');
    disp(w1); disp(w2);
    
    if solver == 2
        disp(allmins{objectmodel,wellquantities+1});
    end
end
%%
toc