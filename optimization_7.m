tic
load ('selfdefine_7.mat');
global wellquantities n coe RpgL RpoL RpwL Po Pg Pw Pig Dg DL gas_provide...
    ytolrealmax_production ytolrealmax_benefit w1 w2 ymin XL xmax objectmodel

%%
% gas injection & liquid production rate
data{1,1} = mmscfd1; data{1,2} = STBd1;
data{2,1} = mmscfd2; data{2,2} = STBd2;
data{3,1} = mmscfd3; data{3,2} = STBd3;
data{4,1} = mmscfd4; data{4,2} = STBd4;
data{5,1} = mmscfd5; data{5,2} = STBd5;
data{6,1} = mmscfd6; data{6,2} = STBd6;
data{7,1} = mmscfd7; data{7,2} = STBd7;
[wellquantities,~] = size(data);

%%
RpoL = [0.5 0.5 0.5 0.5 0.5 0.5 0.5];  % RpoL is a row vector
for i = 1:wellquantities
    if RpoL(i) < 0 || RpoL(i) > 1
        error('RpoL(i) must be between [0,1]');
    end
end
RpwL = 1 - RpoL;
Po = 100;     % price of oil (dollar/STB)
Pg = 2000;    % price of gas (dollar/mmscf)
Pw = 0.5;     % price of water dealing (dollar/STB)
Pig = 200;    % price of gas injectioin (dollar/mmscf)
Dg = 20;      % total outlet gas dealing capacity
DL = 3000;    % total outlet liquid dealing capacity
x0 = ones(1,wellquantities);    % initial points of fmincon
A = ones(1,wellquantities); gas_provide = 5;   % the constraint of the sum of gas injection of each well
if Po <= 0 || Pg <= 0 || Pw <= 0 || Pig <= 0 || Dg <= 0 || DL <= 0 || gas_provide <= 0
    error('Po,Pg,Pw,Pig,Dg,DL,b must be positive')
elseif Po <= Pw
    error('Po must be larger than Pw')
end

%%
% evaluate production gas&liquid rate(RpgL) for nonconstraints
% distribute memory for wcut in advance to enhance computation speed
wcut = zeros(wellquantities,1);
RpgL = zeros(wellquantities,1);
% gas&oil rate on surface standard condition, the unit of GOR and RpgL is: scf/STB
GOR(1) = 50; GOR(2) = 50; GOR(3) = 100; GOR(4) = 100; GOR(5) = 100; GOR(6) = 50; GOR(7) = 100;  % scf/STB
% WCut on surface standard condition
wcut(1) = 0.5; wcut(2) = 0.5; wcut(3) = 0.5; wcut(4) = 0.5; wcut(5) = 0.5; wcut(6) = 0.5; wcut(7) = 0.5;
% production gas&oil ratio on surface standard condition, scf/STB
disp('production gas liquid ratio of each well are below:(scf/STB)')
for i = 1:wellquantities
    if GOR(i) <= 0
        error('GOR must be > 0')
    end
    if wcut(i) < 0 || wcut(i) > 1
        error('wcut must be in [0,1]')
    end
    RpgL(i) = (1 - wcut(i)) * GOR(i);
    disp(RpgL(i));     % scf/STB
end

%%
n = 6;                % set fitting power
if n <= 0
    error('n must be a positive')
end
%%
objectmodeltype = ['production' 'benefit' 'mix'];
% watch out!!!
objectmodel = 1;     % objectmodel must be 1,2,3 orderly when you first run matlab !!!
if objectmodel ~= 1 && objectmodel ~= 2 && objectmodel ~= 3
    error('objectmodel must be 1,2 or 3')
end
% notice!!!  when w1 is close to 0 or 1, because the divisor
% 'ytolrealmax_benefit' and 'ytolrealmax_production' are large in some extent, in the process of calculation
% it will lead to some extent of spreading error, so the minmum of object function
% 'w1 * f1 / ytolrealmax_production' and 'w2 * f2 / ytolrealmax_benefit' can't be exact with
% the minimum of object function 'f1' and 'f2'
if objectmodel == 3
    w1 = 0.5 ; w2 = 1 - w1;
    if w1 < 0 || w1 > 1
        error('w1 must be between [0,1]');
    end
end

%%
solvertype = ['fmincon' 'globalsearch' 'interiorpointPF'];
solver = 2;
if solver ~= 1 && solver ~= 2  && solver ~= 3
    error('solver must be 1,2,3')
end
% if solver == 3
%     x0 = [2 2 1.5 1.5 1.5 0.5];
% r0 = 0.1; c = 0.1;
% % for i = 1:wellquantities
% %    syms gas_injection(i);
% % end
% % syms gas_injection1 gas_injection2 gas_injection3 gas_injection4 gas_injection5 gas_injection6
% % sym gas_injection;
% % variable = [gas_injection(1) gas_injection(2) gas_injection(3) gas_injection(4) gas_injection(5) gas_injection(6)];
% end

%%
algorithmoffmincontype = ['interior-point' 'SQP' 'active-set'];
algorithmoffmincon = 1;     %  alg must be 1,2 or 3
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
plotfunctions_fmincontype = ['x' 'funccount' 'fval' 'constrviolation' 'stepsize' 'firstorderopt'];
plotfunctions_fmincon = 2;
if plotfunctions_fmincon ~= 1 && plotfunctions_fmincon ~= 2 && plotfunctions_fmincon ~= 3 &&...
        plotfunctions_fmincon ~= 4 && plotfunctions_fmincon ~= 5 && plotfunctions_fmincon ~= 6
    error('plotfunctions must be 1,2,3,4,5,6')
end

plotfunctions_globalsearchtype = ['bestf' 'funccount'];
plotfunctions_globalsearch = 2;
if plotfunctions_fmincon ~= 1 && plotfunctions_fmincon ~= 2
    error('plotfunctions must be 1,2')
end

languagetype = ['English' 'Chinese'];
language = 2;
if language ~= 1 && language ~= 2
    error('language must be 1,2')
end

injectiongastypename = ['g' 'n'];    % 'g' denote natural gas, 'n' denote nitrogen gas
injectiongastype = 1;
if injectiongastype ~= 1 && injectiongastype ~= 2
    error('injectiongastype must be 1,2')
end
%%









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
            options{i,j} = optimoptions(@fmincon,'Algorithm',algorithmoffmincon,...
                'DerivativeCheck','on',...     % 'display','iter-detailed', 'Diagnostics','on',
                'FinDiffType','central','FunValCheck','on','MaxFunEvals',10000,'MaxIter',3000,...
                'TolCon',1e-12,'TolFun',1e-12,'TolX',1e-20);
        end
    end
end

if graph == 2    % to show 'PlotFcns of fmincon' of a certain well
    options{objectmodel,acertainwell} = optimoptions(@fmincon,'Algorithm',algorithmoffmincon,...
        'display','iter-detailed','DerivativeCheck','on','Diagnostics','on',...
        'FinDiffType','central','FunValCheck','on','MaxFunEvals',10000,'MaxIter',3000,...
        'TolCon',1e-12,'TolFun',1e-12,'TolX',1e-20,'PlotFcns',plotfunctions_fmincon);
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
%solve extremum of each well during gas injection area

if objectmodel == 1     % maximum production object function
    problem{objectmodel,1} = createOptimProblem('fmincon','objective',@f1_production,'x0',ones(1,1),...
        'lb',XL(1),'ub',XR(1),'options',options{objectmodel,1});
    problem{objectmodel,2} = createOptimProblem('fmincon','objective',@f2_production,'x0',ones(1,1),...
        'lb',XL(2),'ub',XR(2),'options',options{objectmodel,2});
    problem{objectmodel,3} = createOptimProblem('fmincon','objective',@f3_production,'x0',ones(1,1),...
        'lb',XL(3),'ub',XR(3),'options',options{objectmodel,3});
    problem{objectmodel,4} = createOptimProblem('fmincon','objective',@f4_production,'x0',ones(1,1),...
        'lb',XL(4),'ub',XR(4),'options',options{objectmodel,4});
    problem{objectmodel,5} = createOptimProblem('fmincon','objective',@f5_production,'x0',ones(1,1),...
        'lb',XL(5),'ub',XR(5),'options',options{objectmodel,5});
    problem{objectmodel,6} = createOptimProblem('fmincon','objective',@f6_production,'x0',ones(1,1),...
        'lb',XL(6),'ub',XR(6),'options',options{objectmodel,6});
    problem{objectmodel,7} = createOptimProblem('fmincon','objective',@f7_production,'x0',ones(1,1),...
        'lb',XL(7),'ub',XR(7),'options',options{objectmodel,7});
    
    if solver == 1 || 3
        [xmax(objectmodel,1),ymin(objectmodel,1),exitflag(objectmodel,1),output{objectmodel,1}] = fmincon(problem{objectmodel,1});
        [xmax(objectmodel,2),ymin(objectmodel,2),exitflag(objectmodel,2),output{objectmodel,2}] = fmincon(problem{objectmodel,2});
        [xmax(objectmodel,3),ymin(objectmodel,3),exitflag(objectmodel,3),output{objectmodel,3}] = fmincon(problem{objectmodel,3});
        [xmax(objectmodel,4),ymin(objectmodel,4),exitflag(objectmodel,4),output{objectmodel,4}] = fmincon(problem{objectmodel,4});
        [xmax(objectmodel,5),ymin(objectmodel,5),exitflag(objectmodel,5),output{objectmodel,5}] = fmincon(problem{objectmodel,5});
        [xmax(objectmodel,6),ymin(objectmodel,6),exitflag(objectmodel,6),output{objectmodel,6}] = fmincon(problem{objectmodel,6});
        [xmax(objectmodel,7),ymin(objectmodel,7),exitflag(objectmodel,7),output{objectmodel,7}] = fmincon(problem{objectmodel,7});
    end
    
    if solver == 2
        gs = GlobalSearch;
        [xmax(objectmodel,1),ymin(objectmodel,1),exitflag(objectmodel,1),output{objectmodel,1},allmins{objectmodel,1}] = run(gs,problem{objectmodel,1});
        [xmax(objectmodel,2),ymin(objectmodel,2),exitflag(objectmodel,2),output{objectmodel,2},allmins{objectmodel,2}] = run(gs,problem{objectmodel,2});
        [xmax(objectmodel,3),ymin(objectmodel,3),exitflag(objectmodel,3),output{objectmodel,3},allmins{objectmodel,3}] = run(gs,problem{objectmodel,3});
        [xmax(objectmodel,4),ymin(objectmodel,4),exitflag(objectmodel,4),output{objectmodel,4},allmins{objectmodel,4}] = run(gs,problem{objectmodel,4});
        [xmax(objectmodel,5),ymin(objectmodel,5),exitflag(objectmodel,5),output{objectmodel,5},allmins{objectmodel,5}] = run(gs,problem{objectmodel,5});
        [xmax(objectmodel,6),ymin(objectmodel,6),exitflag(objectmodel,6),output{objectmodel,6},allmins{objectmodel,6}] = run(gs,problem{objectmodel,6});
        [xmax(objectmodel,7),ymin(objectmodel,7),exitflag(objectmodel,7),output{objectmodel,7},allmins{objectmodel,7}] = run(gs,problem{objectmodel,7});
    end
end

%%
if objectmodel == 2   % maximum economic benefits object function
    
    problem{objectmodel,1} = createOptimProblem('fmincon','objective',@f1_benefit,'x0',ones(1,1),...
        'lb',XL(1),'ub',XR(1),'options',options{objectmodel,1});
    problem{objectmodel,2} = createOptimProblem('fmincon','objective',@f2_benefit,'x0',ones(1,1),...
        'lb',XL(2),'ub',XR(2),'options',options{objectmodel,2});
    problem{objectmodel,3} = createOptimProblem('fmincon','objective',@f3_benefit,'x0',ones(1,1),...
        'lb',XL(3),'ub',XR(3),'options',options{objectmodel,3});
    problem{objectmodel,4} = createOptimProblem('fmincon','objective',@f4_benefit,'x0',ones(1,1),...
        'lb',XL(4),'ub',XR(4),'options',options{objectmodel,4});
    problem{objectmodel,5} = createOptimProblem('fmincon','objective',@f5_benefit,'x0',ones(1,1),...
        'lb',XL(5),'ub',XR(5),'options',options{objectmodel,5});
    problem{objectmodel,6} = createOptimProblem('fmincon','objective',@f6_benefit,'x0',ones(1,1),...
        'lb',XL(6),'ub',XR(6),'options',options{objectmodel,6});
    problem{objectmodel,7} = createOptimProblem('fmincon','objective',@f7_benefit,'x0',ones(1,1),...
        'lb',XL(7),'ub',XR(7),'options',options{objectmodel,7});
    
    if solver == 1 || 3
        [xmax(objectmodel,1),ymin(objectmodel,1),exitflag(objectmodel,1),output{objectmodel,1}] = fmincon(problem{objectmodel,1});
        [xmax(objectmodel,2),ymin(objectmodel,2),exitflag(objectmodel,2),output{objectmodel,2}] = fmincon(problem{objectmodel,2});
        [xmax(objectmodel,3),ymin(objectmodel,3),exitflag(objectmodel,3),output{objectmodel,3}] = fmincon(problem{objectmodel,3});
        [xmax(objectmodel,4),ymin(objectmodel,4),exitflag(objectmodel,4),output{objectmodel,4}] = fmincon(problem{objectmodel,4});
        [xmax(objectmodel,5),ymin(objectmodel,5),exitflag(objectmodel,5),output{objectmodel,5}] = fmincon(problem{objectmodel,5});
        [xmax(objectmodel,6),ymin(objectmodel,6),exitflag(objectmodel,6),output{objectmodel,6}] = fmincon(problem{objectmodel,6});
        [xmax(objectmodel,7),ymin(objectmodel,7),exitflag(objectmodel,7),output{objectmodel,7}] = fmincon(problem{objectmodel,7});
    end
    
    if solver == 2
        gs = GlobalSearch;
        [xmax(objectmodel,1),ymin(objectmodel,1),exitflag(objectmodel,1),output{objectmodel,1},allmins{objectmodel,1}] = run(gs,problem{objectmodel,1});
        [xmax(objectmodel,2),ymin(objectmodel,2),exitflag(objectmodel,2),output{objectmodel,2},allmins{objectmodel,2}] = run(gs,problem{objectmodel,2});
        [xmax(objectmodel,3),ymin(objectmodel,3),exitflag(objectmodel,3),output{objectmodel,3},allmins{objectmodel,3}] = run(gs,problem{objectmodel,3});
        [xmax(objectmodel,4),ymin(objectmodel,4),exitflag(objectmodel,4),output{objectmodel,4},allmins{objectmodel,4}] = run(gs,problem{objectmodel,4});
        [xmax(objectmodel,5),ymin(objectmodel,5),exitflag(objectmodel,5),output{objectmodel,5},allmins{objectmodel,5}] = run(gs,problem{objectmodel,5});
        [xmax(objectmodel,6),ymin(objectmodel,6),exitflag(objectmodel,6),output{objectmodel,6},allmins{objectmodel,6}] = run(gs,problem{objectmodel,6});
        [xmax(objectmodel,7),ymin(objectmodel,7),exitflag(objectmodel,7),output{objectmodel,7},allmins{objectmodel,7}] = run(gs,problem{objectmodel,7});
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
    problem{objectmodel-2,6} = createOptimProblem('fmincon','objective',@f6_production,'x0',ones(1,1),...
        'lb',XL(6),'ub',XR(6),'options',options{objectmodel-2,6});
    problem{objectmodel-2,7} = createOptimProblem('fmincon','objective',@f7_production,'x0',ones(1,1),...
        'lb',XL(7),'ub',XR(7),'options',options{objectmodel-2,7});
    
    problem{objectmodel-1,1} = createOptimProblem('fmincon','objective',@f1_benefit,'x0',ones(1,1),...
        'lb',XL(1),'ub',XR(1),'options',options{objectmodel-1,1});
    problem{objectmodel-1,2} = createOptimProblem('fmincon','objective',@f2_benefit,'x0',ones(1,1),...
        'lb',XL(2),'ub',XR(2),'options',options{objectmodel-1,2});
    problem{objectmodel-1,3} = createOptimProblem('fmincon','objective',@f3_benefit,'x0',ones(1,1),...
        'lb',XL(3),'ub',XR(3),'options',options{objectmodel-1,3});
    problem{objectmodel-1,4} = createOptimProblem('fmincon','objective',@f4_benefit,'x0',ones(1,1),...
        'lb',XL(4),'ub',XR(4),'options',options{objectmodel-1,4});
    problem{objectmodel-1,5} = createOptimProblem('fmincon','objective',@f5_benefit,'x0',ones(1,1),...
        'lb',XL(5),'ub',XR(5),'options',options{objectmodel-1,5});
    problem{objectmodel-1,6} = createOptimProblem('fmincon','objective',@f6_benefit,'x0',ones(1,1),...
        'lb',XL(6),'ub',XR(6),'options',options{objectmodel-1,6});
    problem{objectmodel-1,7} = createOptimProblem('fmincon','objective',@f7_benefit,'x0',ones(1,1),...
        'lb',XL(7),'ub',XR(7),'options',options{objectmodel-1,7});
    
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
    problem{objectmodel,6} = createOptimProblem('fmincon','objective',@f6_mix,'x0',ones(1,1),...
        'lb',XL(6),'ub',XR(6),'options',options{objectmodel,6});
    problem{objectmodel,7} = createOptimProblem('fmincon','objective',@f7_mix,'x0',ones(1,1),...
        'lb',XL(7),'ub',XR(7),'options',options{objectmodel,7});
    
    if solver == 1 || 3
        [xmax(objectmodel-2,1),ymin(objectmodel-2,1),exitflag(objectmodel-2,1),output{objectmodel-2,1}] = fmincon(problem{objectmodel-2,1});
        [xmax(objectmodel-2,2),ymin(objectmodel-2,2),exitflag(objectmodel-2,2),output{objectmodel-2,2}] = fmincon(problem{objectmodel-2,2});
        [xmax(objectmodel-2,3),ymin(objectmodel-2,3),exitflag(objectmodel-2,3),output{objectmodel-2,3}] = fmincon(problem{objectmodel-2,3});
        [xmax(objectmodel-2,4),ymin(objectmodel-2,4),exitflag(objectmodel-2,4),output{objectmodel-2,4}] = fmincon(problem{objectmodel-2,4});
        [xmax(objectmodel-2,5),ymin(objectmodel-2,5),exitflag(objectmodel-2,5),output{objectmodel-2,5}] = fmincon(problem{objectmodel-2,5});
        [xmax(objectmodel-2,6),ymin(objectmodel-2,6),exitflag(objectmodel-2,6),output{objectmodel-2,6}] = fmincon(problem{objectmodel-2,6});
        [xmax(objectmodel-2,7),ymin(objectmodel-2,7),exitflag(objectmodel-2,7),output{objectmodel-2,7}] = fmincon(problem{objectmodel-2,7});
        
        [xmax(objectmodel-1,1),ymin(objectmodel-1,1),exitflag(objectmodel-1,1),output{objectmodel-1,1}] = fmincon(problem{objectmodel-1,1});
        [xmax(objectmodel-1,2),ymin(objectmodel-1,2),exitflag(objectmodel-1,2),output{objectmodel-1,2}] = fmincon(problem{objectmodel-1,2});
        [xmax(objectmodel-1,3),ymin(objectmodel-1,3),exitflag(objectmodel-1,3),output{objectmodel-1,3}] = fmincon(problem{objectmodel-1,3});
        [xmax(objectmodel-1,4),ymin(objectmodel-1,4),exitflag(objectmodel-1,4),output{objectmodel-1,4}] = fmincon(problem{objectmodel-1,4});
        [xmax(objectmodel-1,5),ymin(objectmodel-1,5),exitflag(objectmodel-1,5),output{objectmodel-1,5}] = fmincon(problem{objectmodel-1,5});
        [xmax(objectmodel-1,6),ymin(objectmodel-1,6),exitflag(objectmodel-1,6),output{objectmodel-1,6}] = fmincon(problem{objectmodel-1,6});
        [xmax(objectmodel-1,7),ymin(objectmodel-1,7),exitflag(objectmodel-1,7),output{objectmodel-1,7}] = fmincon(problem{objectmodel-1,7});
        
        [xmax(objectmodel,1),ymin(objectmodel,1),exitflag(objectmodel,1),output{objectmodel,1}] = fmincon(problem{objectmodel,1});
        [xmax(objectmodel,2),ymin(objectmodel,2),exitflag(objectmodel,2),output{objectmodel,2}] = fmincon(problem{objectmodel,2});
        [xmax(objectmodel,3),ymin(objectmodel,3),exitflag(objectmodel,3),output{objectmodel,3}] = fmincon(problem{objectmodel,3});
        [xmax(objectmodel,4),ymin(objectmodel,4),exitflag(objectmodel,4),output{objectmodel,4}] = fmincon(problem{objectmodel,4});
        [xmax(objectmodel,5),ymin(objectmodel,5),exitflag(objectmodel,5),output{objectmodel,5}] = fmincon(problem{objectmodel,5});
        [xmax(objectmodel,6),ymin(objectmodel,6),exitflag(objectmodel,6),output{objectmodel,6}] = fmincon(problem{objectmodel,6});
        [xmax(objectmodel,7),ymin(objectmodel,7),exitflag(objectmodel,7),output{objectmodel,7}] = fmincon(problem{objectmodel,7});
    end
    
    if solver == 2
        gs = GlobalSearch;
        
        [xmax(objectmodel-2,1),ymin(objectmodel-2,1),exitflag(objectmodel-2,1),output{objectmodel-2,1}] = run(gs,problem{objectmodel-2,1});
        [xmax(objectmodel-2,2),ymin(objectmodel-2,2),exitflag(objectmodel-2,2),output{objectmodel-2,2}] = run(gs,problem{objectmodel-2,2});
        [xmax(objectmodel-2,3),ymin(objectmodel-2,3),exitflag(objectmodel-2,3),output{objectmodel-2,3}] = run(gs,problem{objectmodel-2,3});
        [xmax(objectmodel-2,4),ymin(objectmodel-2,4),exitflag(objectmodel-2,4),output{objectmodel-2,4}] = run(gs,problem{objectmodel-2,4});
        [xmax(objectmodel-2,5),ymin(objectmodel-2,5),exitflag(objectmodel-2,5),output{objectmodel-2,5}] = run(gs,problem{objectmodel-2,5});
        [xmax(objectmodel-2,6),ymin(objectmodel-2,6),exitflag(objectmodel-2,6),output{objectmodel-2,6}] = run(gs,problem{objectmodel-2,6});
        [xmax(objectmodel-2,7),ymin(objectmodel-2,7),exitflag(objectmodel-2,7),output{objectmodel-2,7}] = run(gs,problem{objectmodel-2,7});
        
        [xmax(objectmodel-1,1),ymin(objectmodel-1,1),exitflag(objectmodel-1,1),output{objectmodel-1,1}] = run(gs,problem{objectmodel-1,1});
        [xmax(objectmodel-1,2),ymin(objectmodel-1,2),exitflag(objectmodel-1,2),output{objectmodel-1,2}] = run(gs,problem{objectmodel-1,2});
        [xmax(objectmodel-1,3),ymin(objectmodel-1,3),exitflag(objectmodel-1,3),output{objectmodel-1,3}] = run(gs,problem{objectmodel-1,3});
        [xmax(objectmodel-1,4),ymin(objectmodel-1,4),exitflag(objectmodel-1,4),output{objectmodel-1,4}] = run(gs,problem{objectmodel-1,4});
        [xmax(objectmodel-1,5),ymin(objectmodel-1,5),exitflag(objectmodel-1,5),output{objectmodel-1,5}] = run(gs,problem{objectmodel-1,5});
        [xmax(objectmodel-1,6),ymin(objectmodel-1,6),exitflag(objectmodel-1,6),output{objectmodel-1,6}] = run(gs,problem{objectmodel-1,6});
        [xmax(objectmodel-1,7),ymin(objectmodel-1,7),exitflag(objectmodel-1,7),output{objectmodel-1,7}] = run(gs,problem{objectmodel-1,7});
        
        [xmax(objectmodel,1),ymin(objectmodel,1),exitflag(objectmodel,1),output{objectmodel,1}] = run(gs,problem{objectmodel,1});
        [xmax(objectmodel,2),ymin(objectmodel,2),exitflag(objectmodel,2),output{objectmodel,2}] = run(gs,problem{objectmodel,2});
        [xmax(objectmodel,3),ymin(objectmodel,3),exitflag(objectmodel,3),output{objectmodel,3}] = run(gs,problem{objectmodel,3});
        [xmax(objectmodel,4),ymin(objectmodel,4),exitflag(objectmodel,4),output{objectmodel,4}] = run(gs,problem{objectmodel,4});
        [xmax(objectmodel,5),ymin(objectmodel,5),exitflag(objectmodel,5),output{objectmodel,5}] = run(gs,problem{objectmodel,5});
        [xmax(objectmodel,6),ymin(objectmodel,6),exitflag(objectmodel,6),output{objectmodel,6}] = run(gs,problem{objectmodel,6});
        [xmax(objectmodel,7),ymin(objectmodel,7),exitflag(objectmodel,7),output{objectmodel,7}] = run(gs,problem{objectmodel,7});
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
if graph == 1 && objectmodel == 1
    figure;
    hold on;
    for i = 1:wellquantities
        if abs( ymax(1,i) - max(data{i,2}) ) / max(data{i,2}) < 0.05     % if calculation optimization value...
            % is trustable
            disp 'single well maximum production(STBd) are';
            disp(ymax(1,i));
            disp 'single well minimum injection(mmscfd) and maximum injection(mmscfd) are';
            disp(XL(i)); disp(xmax(1,i));
            plot(data{i,1},data{i,2});
        end
    end
    if language == 1
        title('injection & production curve');
        xlabel('gas injection/(mmscf/d)');
        ylabel('liquid production/(STB/d)');
    elseif language == 2
        %         title('注气产液量实际曲线图');
        %         xlabel('注气量/(万方/天)');
        %         ylabel('产液量/(桶/天)');
    end
    hold off;
end

% elseif objectmodel == 2      % maximum economic benefits object function
% for i = 1:wellquantities
%         disp 'single well maximum economic benefits (dollar/d) are';
%         disp(ymax_benefit(i));
%         disp 'single well minimum injection(mmscfd) and maximum injection(mmscfd) are';
%         disp(XL(i)); disp(xmax_benefit(i));
%         plot(data{i,1},data{i,2});
% end
% title('injection & benefits curve');
% xlabel('gas injection (unit:mmscfd)'); ylabel('economic benefits (unit:dollar/d)');
%
% elseif objectmodel == 3      % maximum mixed object function
% for i = 1:wellquantities
%         disp 'single well maximum the mixed object function (dimensionless) are';
%         disp(ymax_mix(i));
%         disp 'single well minimum injection(mmscfd) and maximum injection(mmscfd) are';
%         disp(XL(i)); disp(xmax_mix(i));
%         plot(data{i,1},data{i,2});
% end
% title('injection & dimensionless curve');
% xlabel('gas injection (unit:mmscfd)'); ylabel('the mixed object function (unit:none)');
% end

%%
for i = 1:length(objectmodeltype)
    if solver == 1
        if graph == 1
            options{i,wellquantities+1} = optimoptions(@fmincon,'Algorithm',algorithmoffmincon,...
                'Diagnostics','on','DerivativeCheck','on','display','iter-detailed',...
                'FinDiffType','central','FunValCheck','on','MaxFunEvals',10000,'MaxIter',3000,...
                'TolCon',1e-12,'TolFun',1e-12,'TolX',1e-20,'PlotFcns',plotfunctions_fmincon);
        end
        if graph == 2     % not to show 'PlotFcns' of 'fmincon'
            options{i,wellquantities+1} = optimoptions(@fmincon,'Algorithm',algorithmoffmincon,...
                'Diagnostics','on','DerivativeCheck','on','display','iter-detailed',...
                'FinDiffType','central','FunValCheck','on','MaxFunEvals',10000,'MaxIter',3000,...
                'TolCon',1e-12,'TolFun',1e-12,'TolX',1e-20);
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
gs = GlobalSearch('Display','iter','PlotFcns',plotfunctions_globalsearch);     % @gsplotfunccount   @gsplotbestf

if objectmodel == 1
    LB = [XL(1) XL(2) XL(3) XL(4) XL(5) XL(6) XL(7)];
    UB_production = [xmax(objectmodel,1) xmax(objectmodel,2) xmax(objectmodel,3),...
        xmax(objectmodel,4) xmax(objectmodel,5) xmax(objectmodel,6) xmax(objectmodel,7)]; % constraint of gas injection of each well
    
    if solver == 1
        % use fmincon to solve optimizaiton distribution
        [xmax_total{objectmodel},ymin_total{objectmodel},exitflag(objectmodel,wellquantities+1),...
            output{objectmodel,wellquantities+1},lambda,grad,hessian] = fmincon(@objfunc_production,...
            x0, A, gas_provide,[], [], LB, UB_production,@noncons,options{objectmodel,wellquantities+1});                        % solve mathematics mode
    end
    if solver == 2
        % use GlobalSearch to solve optimizaiton distribution
        problem{objectmodel,wellquantities+1} = createOptimProblem('fmincon','objective',...
            @objfunc_production,'x0',x0,'lb',LB,'ub',UB_production,'Aineq',A,...
            'bineq',gas_provide,'nonlcon',@noncons,'options',options{objectmodel,wellquantities+1});
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
    Ginjmax = 0; Lprodmax = 0;
    for i = 1:wellquantities
        Ginjmax = Ginjmax + xmax(objectmodel,i); Lprodmax = Lprodmax + ymax(objectmodel,i);
    end
    if solver == 1 || solver == 2
        disp(exitflag(objectmodel,wellquantities+1)); disp(output{objectmodel,wellquantities+1});
    end
    if solver == 1
        disp(lambda); disp(grad); disp(hessian);
    end
    disp('the constraint of total gas provide ability of all wells(mmscfd)');
    disp(gas_provide);
    disp('unconstraint maxmum gas injction rate(mmscfd)');
    disp(Ginjmax);
    toptsum = sum(xmax_total{objectmodel});
    disp('the sum of actual gas injction rate of all wells(mmscfd)');
    disp(toptsum);
    disp('actual gas injction rate of each well(mmscfd)');
    disp(xmax_total{objectmodel});
    disp('the constraint of total gas deal ability of all wells(mmscfd)');
    disp(Dg);
    producegas = 0;
    for i = 1:wellquantities
        producegas = producegas + RpgL(i)*ymax(objectmodel,i)/10000;
    end
    if injectiongastype == 1
        if producegas > Ginjmax
            actualgasneededtobedealt = producegas;
        else
            actualgasneededtobedealt = Ginjmax;
        end
    elseif  injectiongastype == 2
        actualgasneededtobedealt = producegas + Ginjmax;
    end
    disp('under unconstraint total gas rate needed to be dealt with(mmscfd)');
    disp(actualgasneededtobedealt);
    disp('the constraint of total liquid deal ability of all wells(STBd)');
    disp(DL);
    disp('unconstraint maxmum liquid production rate(STBd)');
    disp(Lprodmax);
    ytolrealmax_production = -ymin_total{objectmodel};
    disp('the sum of actual maxmum liquid production rate(STBd)');
    disp(ytolrealmax_production);
    if solver == 2
        disp(allmins{objectmodel,wellquantities+1});
    end
end

%%
if objectmodel == 2
    % use fmincon to solve optimizaiton distribution
    LB = [XL(1) XL(2) XL(3) XL(4) XL(5) XL(6) XL(7)];
    UB_benefit = [xmax(objectmodel,1) xmax(objectmodel,objectmodel) xmax(objectmodel,3),...
        xmax(objectmodel,4) xmax(objectmodel,5) xmax(objectmodel,6) xmax(objectmodel,7)];% constraint of gas injection of each well
    
    if solver == 1
        [xmax_total{objectmodel},ymin_total{objectmodel},exitflag(objectmodel,wellquantities+1),...
            output{objectmodel,wellquantities+1},lambda,grad,hessian] = fmincon(@objfunc_benefit,...
            x0, A,gas_provide,[], [], LB, UB_benefit,@noncons,options{objectmodel,wellquantities+1});                        % solve mathematics mode
    end
    if solver == 2
        problem{objectmodel,wellquantities+1} = createOptimProblem('fmincon','objective',@objfunc_benefit,...
            'x0',x0,'lb',LB,'ub',UB_benefit,'Aineq',A,'bineq',gas_provide,'nonlcon',@noncons,...
            'options',options{objectmodel,wellquantities+1});
        [xmax_total{objectmodel},ymin_total{objectmodel},exitflag(objectmodel,wellquantities+1),...
            output{objectmodel,wellquantities+1},allmins{objectmodel,wellquantities+1}]...
            = run(gs,problem{objectmodel,wellquantities+1});
    end
    
    % max gas injection rate at max economic benefit and unconstraint max economic benefit of block wells
    Ginjmax = 0; Lprodmax = 0;
    for i = 1:wellquantities
        Ginjmax = Ginjmax + xmax(objectmodel,i); Lprodmax = Lprodmax + ymax(objectmodel,i);
    end
    disp(exitflag(objectmodel,wellquantities+1)); disp(output{objectmodel,wellquantities+1});
    if solver == 1
        disp(lambda); disp(grad); disp(hessian);
    end
    disp('the constraint of total gas provide ability of all wells(mmscfd)');
    disp(gas_provide);
    disp('the sum of gas injection rate of block wells under unconstraint max economic benefits(mmscfd)');
    disp(Ginjmax);
    toptsum = sum(xmax_total{objectmodel});
    disp('the sum of actual gas injction rate of all wells(mmscfd)');
    disp(toptsum);
    disp('actual gas injction rate of each well under constraints(mmscfd)');
    disp(xmax_total{objectmodel});
    disp('the constraint of total gas deal ability of all wells(mmscfd)');
    disp(Dg);
    %     producegas = 0;
    %     for i = 1:wellquantities
    %         producegas = producegas + RpgL(i)*ymax(objectmodel,i)/10000;
    %     end
    %     actualgasneededtobedealt = toptsum + producegas;
    %     disp('gas rate needed to be dealt actually');
    %     disp(actualgasneededtobedealt);
    disp('the constraint of total liquid deal ability of all wells(STBd)');
    disp(DL);
    disp('unconstraint maxmum economic benefits(dollar/d)');
    disp(Lprodmax);
    ytolrealmax_benefit = -ymin_total{objectmodel};
    disp('actual maxmum economic benefit(dollar/d)');
    disp(ytolrealmax_benefit);
    if solver == 2
        disp(allmins{objectmodel,wellquantities+1});
    end
end
%%

if objectmodel == 3
    % use fmincon to solve optimizaiton distribution
    % the constraints(lb,ub) of gas injection of each well
    LB = [XL(1) XL(2) XL(3) XL(4) XL(5) XL(6) XL(7)];
    UB_production = [xmax(objectmodel-2,1) xmax(objectmodel-2,2) xmax(objectmodel-2,3),...
        xmax(objectmodel-2,4) xmax(objectmodel-2,5) xmax(objectmodel-2,6) xmax(objectmodel-2,7)];
    UB_benefit = [xmax(objectmodel-1,1) xmax(objectmodel-1,2) xmax(objectmodel-1,3),...
        xmax(objectmodel-1,4) xmax(objectmodel-1,5) xmax(objectmodel-1,6) xmax(objectmodel-1,7)];
    UB_mix = [xmax(objectmodel,1) xmax(objectmodel,2) xmax(objectmodel,3) xmax(objectmodel,4),...
        xmax(objectmodel,5) xmax(objectmodel,6) xmax(objectmodel,7)];
    
    if solver == 1
        [xmax_total{objectmodel-2},ymin_total{objectmodel-2},exitflag(objectmodel-2,wellquantities+1),...
            output{objectmodel-2,wellquantities+1}] = fmincon(@objfunc_production,x0,...
            A,gas_provide,[], [], LB, UB_production,@noncons,options{objectmodel-2,wellquantities+1});
        ytolrealmax_production = -ymin_total{objectmodel-2};
        
        [xmax_total{objectmodel-1},ymin_total{objectmodel-1},exitflag(objectmodel-1,wellquantities+1),...
            output{objectmodel-1,wellquantities+1}] = fmincon(@objfunc_benefit,x0,...
            A,gas_provide,[], [], LB, UB_benefit,@noncons,options{objectmodel-1,wellquantities+1});
        ytolrealmax_benefit = -ymin_total{objectmodel-1};
        
        [xmax_total{objectmodel},ymin_total{objectmodel},exitflag(objectmodel,wellquantities+1),...
            output{objectmodel,wellquantities+1},lambda,grad,hessian] = fmincon(@objfunc_mix,x0,...
            A,gas_provide,[], [], LB, UB_mix,@noncons,options{objectmodel,wellquantities+1});
    end
    
    if solver == 2
        
        problem{objectmodel-2,wellquantities+1} = createOptimProblem('fmincon',...
            'objective',@objfunc_production,'x0',x0,'lb',LB,'ub',UB_production,'Aineq',A,...
            'bineq',gas_provide,'nonlcon',@noncons,'options',options{objectmodel-2,wellquantities+1});
        [xmax_total{objectmodel-2},ymin_total{objectmodel-2},...
            exitflag(objectmodel-2,wellquantities+objectmodel-2),...
            output{objectmodel-2,wellquantities+1},allmins{objectmodel-2,wellquantities+1}]...
            = run(gs,problem{objectmodel-2,wellquantities+1});
        ytolrealmax_production = -ymin_total{objectmodel-2};
        
        problem{objectmodel-1,wellquantities+1} = createOptimProblem('fmincon',...
            'objective',@objfunc_benefit,'x0',x0,'lb',LB,'ub',UB_benefit,'Aineq',A,...
            'bineq',gas_provide,'nonlcon',@noncons,'options',options{objectmodel-1,wellquantities+1});
        [xmax_total{objectmodel-1},ymin_total{objectmodel-1},exitflag(objectmodel-1,wellquantities+1),...
            output{objectmodel-1,wellquantities+1},allmins{objectmodel-1,wellquantities+1}]...
            = run(gs,problem{objectmodel-1,wellquantities+1});
        ytolrealmax_benefit = -ymin_total{objectmodel-1};
        
        problem{objectmodel,wellquantities+1} = createOptimProblem('fmincon',...
            'objective',@objfunc_mix,'x0',x0,'lb',LB,'ub',UB_mix,'Aineq',A,'bineq',gas_provide,...
            'nonlcon',@noncons,'options',options{objectmodel,wellquantities+1});
        [xmax_total{objectmodel},ymin_total{objectmodel},exitflag(objectmodel,wellquantities+1),...
            output{objectmodel,wellquantities+1},allmins{objectmodel,wellquantities+1}]...
            = run(gs,problem{objectmodel,wellquantities+1});
    end
    
    % max gas injection rate at max economic benefit and unconstraint max economic benefit of block wells
    Ginjmax = 0; Lprodmax = 0;
    for i = 1:wellquantities
        Ginjmax = Ginjmax + xmax(objectmodel,i); Lprodmax = Lprodmax + ymax(objectmodel,i);
    end
    disp(exitflag(objectmodel,wellquantities+1)); disp(output{objectmodel,wellquantities+1});
    if solver == 1
        disp(lambda); disp(grad); disp(hessian);
    end
    disp('the constraint of total gas provide ability of all wells(mmscfd)');
    disp(gas_provide);
    disp('the sum of gas injection rate of block wells under unconstraint maximum mixed object function(mmscfd)');
    disp(Ginjmax);
    toptsum = sum(xmax_total{objectmodel});
    disp('the sum of actual gas injction rate of all wells(mmscfd)');
    disp(toptsum);
    disp('actual gas injction rate of each well with constraint maxmum mixed object function(mmscfd)');
    disp(xmax_total{objectmodel});
    disp('the constraint of total gas deal ability of all wells(mmscfd)');
    disp(Dg);
    producegas = 0;
    for i = 1:wellquantities
        producegas = producegas + RpgL(i)*ymax(1,i)/10000;
    end
    if injectiongastype == 1
        if producegas > Ginjmax
            actualgasneededtobedealt = producegas;
        else
            actualgasneededtobedealt = Ginjmax;
        end
    elseif  injectiongastype == 2
        actualgasneededtobedealt = producegas + Ginjmax;
    end
    disp('under unconstraint total gas rate needed to be dealt with(mmscfd)');
    disp(actualgasneededtobedealt);
    disp('the constraint of total liquid deal ability of all wells(STBd)');
    disp(DL);
    disp('the value of unconstraint maxmum mixed object function(dimensionless)');
    disp(Lprodmax);
    ytolrealmax_mix = -ymin_total{objectmodel};
    disp('the actual value of maxmum mixed object function(dimensionless)');
    disp(ytolrealmax_mix);
    disp('w1 and w2 are:');
    disp(w1); disp(w2);
    if solver == 2
        disp(allmins{objectmodel,wellquantities+1});
    end
end
%%
toc