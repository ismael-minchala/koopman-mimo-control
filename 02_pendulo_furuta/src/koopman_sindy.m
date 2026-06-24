%% - ******************** -Koopman approximation- ********************** %%
%% - ********************** -Pendulum System- ************************** %%
%% - ******************** -Juan Francisco Duran S- ********************* %%
clear; close all; clc
set(groot,...
'defaulttextinterpreter','latex',...
'defaultAxesTickLabelInterpreter','latex',...
'defaultLegendInterpreter','latex')

%% - ********************** -Define Variables- ************************* %%
n = 2;                              %Number of states
m = 1;                              %Number of control inputs
deltaT = 0.001;

data = load(fullfile('..','datos','pendulum_time_series.mat'));
test = load(fullfile('..','datos','pendulum_time_series_test.mat'));

%% - ************************* -Collect Data- ************************** %%
Npoints = 10e3 + 1;
disp('Starting data collection \n')
t = data.ans(1,1:Npoints);
alpha = (data.ans(2,1:Npoints))*(pi/180);
theta = (data.ans(3,1:Npoints))*(pi/180);
theta_dot = (data.ans(4,1:Npoints));
alpha_dot = (data.ans(5,1:Npoints));
voltage = (data.ans(6,1:Npoints));

t_test = test.ans(1,:);
alpha_test = test.ans(2,:)*(pi/180);
theta_test = test.ans(3,:)*(pi/180);
theta_dot_test = test.ans(4,:);
alpha_dot_test = (test.ans(5,:));
voltage_test = test.ans(6,:);

figure(1)
ejes1(1) = subplot(211);
plot(t,alpha,t,theta,'LineWidth',2)
title('Training Set','FontName', 'Times', 'FontSize',15,'FontWeight','normal')
xlabel('Time [s]')
ylabel('Angular position [rad]')
legend(["$Pendulum$"  "$Arm$"])
ax = gca;
set(ax,'FontName', 'Times', 'FontSize',15,'FontWeight','normal')
grid on
set(groot,...
    'defaulttextinterpreter','latex',...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultLegendInterpreter','latex')

grid on
ejes1(2) = subplot(212);
plot(t,voltage,'LineWidth',2)
title('Training Set','FontName', 'Times', 'FontSize',15,'FontWeight','normal')
xlabel('Time [s]')
ylabel('Applied Voltage [V]')
ax = gca;
set(ax,'FontName', 'Times', 'FontSize',15,'FontWeight','normal')
grid on
set(groot,...
    'defaulttextinterpreter','latex',...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultLegendInterpreter','latex')

linkaxes(ejes1,'x')

H = [alpha;theta;alpha_dot;theta_dot];
H_test = [alpha_test;theta_test;alpha_dot_test;theta_dot_test]';

X = H(1:2,1:end-1)';
U = voltage(:,1:end-1)';
X = [X U];
dX = H(1:2,2:end)';
dX(:,n+1) = 0*dX(:,n);
%% - ************************ -SINDY ALGORITHM- ************************ %%
fprintf('Starting SINDy \n')
npol = 10;
usesine = 1;
Theta = build_theta(X,n+1,npol,usesine);
Nlambdas = 100;
lambda = logspace(-6,-2,Nlambdas);
error = zeros(1,Nlambdas);
Ntest = 5e3;
X_pred0 = zeros(1,2);
tic
parfor j = 1 : Nlambdas
    Xi = LSTR(lambda(j),Theta,dX,n);
    X_pred = sim_dyn(Ntest,Xi,X_pred0,npol,usesine,n,voltage_test);
    error(j) = norm(H_test(1:Ntest,1:2)- X_pred, 'fro') / norm(H_test(1:Ntest,1:2), 'fro');
end
toc
[~,index] = min(error);
Xi = LSTR(lambda(index),Theta,dX,n);
X_pred = sim_dyn(Ntest,Xi,X_pred0,npol,usesine,n,voltage_test);
%% results
figure(2)
subplot(211)
plot(t_test(1:Ntest),H_test(1:Ntest,1),t_test(1:Ntest),X_pred(:,1))
subplot(212)
plot(t_test(1:Ntest),H_test(1:Ntest,2),t_test(1:Ntest),X_pred(:,2))

function X_pred = sim_dyn(Ntest,Xi,X_pred0,npol,usesine,n,u)
X_pred(1,:) = X_pred0;
xlift = build_theta([X_pred(1,:) 0],length(X_pred(1,:))+1,npol,usesine);
for i = 1:Ntest-1
    X_pred(i+1,:) = xlift*Xi(:,1:n);
    xlift = build_theta([X_pred(i,:) u(i)],length(X_pred(i,:))+1,npol,usesine);
end
end
function Xi = LSTR(lambda,Theta,dX,n)
Xi = Theta\dX;  % initial guess: Least-squares
for k=1:100
    smallinds = (abs(Xi)<lambda);   % find small coefficients
    Xi(smallinds)=0;                % and threshold
    for ind = 1:n                   % n is state dimension
        biginds = ~smallinds(:,ind);
        % Regress dynamics onto remaining terms to find sparse Xi
        Xi(biginds,ind) = Theta(:,biginds)\dX(:,ind);
    end
end
end

function yout = build_theta(yin,nVars,polyorder,usesine)

n = size(yin,1);
% yout = zeros(n,1+nVars+(nVars*(nVars+1)/2)+(nVars*(nVars+1)*(nVars+2)/(2*3))+11);

ind = 1;
% poly order 0
yout(:,ind) = ones(n,1);
ind = ind+1;

% poly order 1
for i=1:nVars
    yout(:,ind) = yin(:,i);
    ind = ind+1;
end

if(polyorder>=2)
    % poly order 2
    for i=1:nVars
        for j=i:nVars
            yout(:,ind) = yin(:,i).*yin(:,j);
            ind = ind+1;
        end
    end
end

if(polyorder>=3)
    % poly order 3
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k);
                ind = ind+1;
            end
        end
    end
end

if(polyorder>=4)
    % poly order 4
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                for l=k:nVars
                    yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k).*yin(:,l);
                    ind = ind+1;
                end
            end
        end
    end
end

if(polyorder>=5)
    % poly order 5
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                for l=k:nVars
                    for m=l:nVars
                        yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k).*yin(:,l).*yin(:,m);
                        ind = ind+1;
                    end
                end
            end
        end
    end
end

if(usesine)
    for k=1:50
        yout = [yout sin(2*pi*k*yin) cos(2*pi*k*yin)];
    end
end
end