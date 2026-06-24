%% - ******************** -Koopman approximation- ********************** %%
%% - ********************** -Multitank System- ************************* %%
%% - ******************** -Juan Francisco Duran S- ********************* %%
clear; close all; clc
%% - ********************** -Define Variables- ************************* %%
set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultTextFontName', 'Times New Roman');
set(groot, 'defaultLegendFontName', 'Times New Roman');

set(groot, 'defaultAxesFontSize', 18);  
set(groot, 'defaultTextFontSize', 18);
set(groot, 'defaultLegendFontSize', 18);

set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
%% - ************************* -Collect Data- ************************** %%
% load(fullfile('..','datos','Datta.mat'));
load(fullfile('..','datos','Datta.mat'));

deltaT = Datta.Time(2) - Datta.Time(1); 
Data = ((Datta{:,4:4:end}))*100;

Ubig1 = normalize((Datta{:,2:4:end}),1,'center');
Ubig2 = []; %zscore(reshape(Datta{:,3:4:end},1,[]));

numSamples = length(Data)-1;
train = floor(numSamples);

Xnext = Data(:,1)';
Ubig = Ubig1(:,1)';

Xtest = Data(:,2)';
Utest = [Ubig1(:,2)]';

k = 0:1:length(Xnext)-1;
k2 = 0:1:length(Xtest)-1;

figure(1)
subplot(211)
plot(k,Xnext',k,Ubig,'LineWidth',2)
title('Training Set','FontName', 'Times', 'FontSize',15,'FontWeight','normal')
xlabel('Time [s]')
ylabel('states')
legend(["Frequency Deviation" ])
ax = gca;
set(ax,'FontName', 'Times', 'FontSize',15,'FontWeight','normal')
grid on

subplot(212)
plot(k2,Xtest',k2,Utest,'LineWidth',2)
title('Test Set','FontName', 'Times', 'FontSize',15,'FontWeight','normal')
xlabel('Time [s]')
ylabel('states')
legend(["Frequency Deviation"])
ax = gca;
set(ax,'FontName', 'Times', 'FontSize',15,'FontWeight','normal')
grid on

%% - ************************ -BASIS FUNCTION- ************************* %%
X = Xnext(:,1:end-1);
Y = Xnext(:,2:end);
U = Ubig(:,1:end-1);
%% - ******************** -SINDY OPTIM PARAMETERS- ********************* %%
nx = optimizableVariable('nx', [0, 1], 'Transform', 'none','Type','integer');
nu = optimizableVariable('nu', [0, 6], 'Transform', 'none','Type','integer');
wx = optimizableVariable('wx', [0, 17], 'Transform', 'none','Type','integer');
wu = optimizableVariable('wu', [0, 20], 'Transform', 'none','Type','integer');
lambda = optimizableVariable('lambda', [1e-9, 1e-8], 'Transform', 'log');

cost_fn = @(parameters) sindy_opt([parameters.nx,parameters.nu,parameters.wx,parameters.wu,parameters.lambda],X,U,Xtest,Utest,Y);

results = bayesopt(cost_fn, [nx nu wx wu lambda], ...
    'Verbose', 1, ...
    'MaxObjectiveEvaluations', 20, ...
    'IsObjectiveDeterministic', true, ...
    'AcquisitionFunctionName', 'expected-improvement-plus','UseParallel',true);
%% - *********************** -SINDY OPTIMIZED- ************************* %%
parameters = table2array(results.XAtMinObjective);
n = size(Xtest,1);
Theta = build_theta(X,U,parameters(1),parameters(2),parameters(3),parameters(4))';
Y = Y';
lambda = parameters(5);
Xi = Theta\Y;

for k=1:10   
    smallinds = (abs(Xi)<lambda); 
    Xi(smallinds) = 0; 
    for ind = 1:n 
        biginds = ~smallinds(:,ind);
        Xi(biginds,ind) = Theta(:,biginds)\Y(:,ind);
    end
end
%%
clear x
Nsim = 1e3;
Xprobe = Xtest(1:Nsim);
x(1) = Xprobe(1);
xlift = build_theta(x(1),Utest(:,1),parameters(1),parameters(2),parameters(3),parameters(4))';
for i = 2: Nsim
    x(i) = xlift*Xi;
    xlift = build_theta(x(i),Utest(:,i-1),parameters(1),parameters(2),parameters(3),parameters(4))';
    % if mod(i,10) == 0
    %     xlift= build_theta(Xprobe(i),Utest(:,i),parameters(1),parameters(2),parameters(3),parameters(4))';
    % end
end
J = norm(Xprobe - x,'fro') / norm(Xprobe,'fro')


figure(2)

plot([0:Nsim-1]*deltaT,Xprobe,'-b','LineWidth',2); hold on
plot([0:Nsim-1]*deltaT,x,'--r','LineWidth',2); hold on
legend(["True" "Forecast"])
ax = gca;
set(ax,'FontName', 'Times', 'FontSize',15,'FontWeight','normal')
grid on


function J = sindy_opt(parameters,X,U,Xtest,Utest,Y)
n = size(Xtest,1);
Theta = build_theta(X,U,parameters(1),parameters(2),parameters(3),parameters(4))';
Y = Y';
lambda = parameters(5);
Xi = Theta\Y;

for k=1:10    
    smallinds = (abs(Xi)<lambda); 
    Xi(smallinds) = 0; 
    for ind = 1:n 
        biginds = ~smallinds(:,ind);
        Xi(biginds,ind) = Theta(:,biginds)\Y(:,ind);
    end
end

%% - ************************* -Comparison- **************************** %%
clear x
Nsim = 1e3;
Xprobe = Xtest(1:Nsim);
x(1) = Xprobe(1);
xlift = build_theta(x(1),Utest(:,1),parameters(1),parameters(2),parameters(3),parameters(4))';
for i = 2: Nsim
    x(i) = xlift*Xi;
    xlift = build_theta(x(i),Utest(:,i-1),parameters(1),parameters(2),parameters(3),parameters(4))';
    % if mod(i,100) == 0 && i > 1
    %     xlift= build_theta(Xprobe(i),Utest(:,i),parameters(1),parameters(2),parameters(3),parameters(4))';
    % end
end
J = norm(Xprobe - x,'fro') / norm(Xprobe,'fro');
end
%% - ********************* -SINDY algorithm FIT- *********************** %%


function theta = build_theta(X,U,nx,nu,wx,wu)

X =  X';
U = U';
k = 1;

for i = 0: nx
    theta(:,k) = X(:,1).^(i);
    k = k + 1;
end

for m = 0 :nu
    theta(:,k) = U(:,1).^(m);
    k = k + 1;
end
M = theta(:,1:(nx+1));
P = theta(:,(nx+1)+1:end);

result = zeros(size(M,1), min(size(M, 2), size(P, 2)));

for i = 1:size(result, 2)
    result(:, i) = M(:, i) .* P(:, i);
end

for k = 1:wx
    theta = [theta, sin(k*X), cos(k*X)];
end

for k = 1:wu
    theta = [theta, sin(k*U), cos(k*U)];
end

theta = [(theta) (result)];
%
% theta = [theta U]';
theta = theta';

end

