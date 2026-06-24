%% - ******************** -Koopman approximation- ********************** %%
%% - ********************** -Pendulum System- ************************** %%
%% - ******************** -Juan Francisco Duran S- ********************* %%
clear; close all; clc
set(groot,...
'defaulttextinterpreter','latex',...
'defaultAxesTickLabelInterpreter','latex',...
'defaultLegendInterpreter','latex')

%% - ********************** -Define Variables- ************************* %%
n = 4;                              %Number of states
m = 1;                              %Number of control inputs
deltaT = 0.001;
data = load(fullfile('..','datos','pendulum_time_series.mat'));
test = load(fullfile('..','datos','pendulum_time_series_test.mat'));
% data = load('C:\Users\FranciscoDuran\Desktop\SINDy estimation\pendulum_time_series.mat');
% test = load('C:\Users\FranciscoDuran\Desktop\SINDy estimation\pendulum_time_series_test.mat');
% data = load('AJUSTAR_RUTA'  % Ruta original de Windows — ajustar al directorio de datos);
% test = load('AJUSTAR_RUTA'  % Ruta original de Windows — ajustar al directorio de datos);

%% - ************************* -Collect Data- ************************** %%
Npoints = 10e3 + 1;
disp('Starting data collection \n')
t = data.ans(1,1:Npoints);
alpha = (data.ans(2,1:Npoints))*(pi/180);
theta = (data.ans(3,1:Npoints))*(pi/180);
theta_dot = (data.ans(4,1:Npoints));
alpha_dot = (data.ans(5,1:Npoints));
voltage = (data.ans(6,1:Npoints));

N_test = 10e3;
t_test = test.ans(1,1:N_test);
alpha_test = test.ans(2,1:N_test)*(pi/180);
theta_test = test.ans(3,1:N_test)*(pi/180);
theta_dot_test = test.ans(4,1:N_test);
alpha_dot_test = (test.ans(5,1:N_test));
voltage_test = test.ans(6,1:N_test);

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

X = H(:,1:end-1);
Y = H(:,2:end);
U = voltage(:,1:end-1);
%% - ******************** -LIFTING FUNCTIONS -alpha- ******************* %%
x_alpha = sin(alpha);
H1alpha = [hermiteH(0,x_alpha);hermiteH(1,x_alpha)];
x_theta = sin(theta);
H1theta = [hermiteH(0,x_theta);hermiteH(1,x_theta)];

%% - ************************ -SINDY ALGORITHM- ************************ %%
fprintf('Starting Disperse EDMD \n')
npol_states = 1; usesine_states = 1; num_sines_states = 1;
npol_inputs = 1; usesine_inputs = 1; num_sines_inputs = 1;


Xlift =  build_theta(X',n,npol_states,usesine_states,num_sines_states)';
Ylift = build_theta(Y',n,npol_states,usesine_states,num_sines_states)';
Ulift = build_theta(U',m,npol_inputs,usesine_inputs,num_sines_inputs)';
Nlift = size(Xlift, 1);

W = Ylift;
V = [Xlift;Ulift];

VVt = V * V';
WVt = W * V';

ABC = WVt * pinv(VVt);
Alift = ABC(1:Nlift,1:Nlift);
Blift = ABC(1:Nlift,Nlift+1:end);
Clift = X*pinv(Xlift);
if isempty(Blift)
    Blift = zeros(Nlift,1);
end
%%
alpha = optimizableVariable('alpha',[1e-5,1]);
alpha2 = optimizableVariable('alpha2',[1e-5,1]);
% lambda = optimizableVariable('lambda',[0,10]);
% lambda2 = optimizableVariable('lambda2',[0,10]);

fun = @(alpha)FO(alpha,Nlift,V,W,Xlift,t_test,voltage_test,X,n,npol_states,usesine_states,npol_inputs,usesine_inputs,num_sines_inputs,num_sines_states,H_test);

results = bayesopt(fun,[alpha alpha2],'IsObjectiveDeterministic',true,'Verbose',0,'UseParallel',true);

[Alift2,Blift2,Clift2] = lasso_regression(Nlift,results.XAtMinObjective,V,W,Xlift,X);
sys = ss(Alift,Blift,Clift,0,deltaT);
sys2 = ss(Alift2,Blift2,Clift2,0,deltaT);
%%
Ypred = sim_dymamics(t_test,Alift,Blift,Clift,Nlift,npol_states,usesine_states,npol_inputs,usesine_inputs,num_sines_inputs,num_sines_states,n,voltage_test);
Ypred2 = sim_dymamics(t_test,Alift2,Blift2,Clift2,Nlift,npol_states,usesine_states,npol_inputs,usesine_inputs,num_sines_inputs,num_sines_states,n,voltage_test);

figure(4)
subplot(211)
plot(t_test,H_test(:,1),t_test,Ypred(:,1),t_test,Ypred2(:,1))
legend(["Real" "normal EDMD"])
subplot(212)
plot(t_test,H_test(:,2),t_test,Ypred(:,2),t_test,Ypred2(:,2))
legend(["Real" "normal EDMD"])
figure(5)
subplot(211)
plot(t_test,H_test(:,1),t_test,Ypred2(:,1))
legend(["Real" "Sparse EDMD"])

subplot(212)
plot(t_test,H_test(:,2),t_test,Ypred2(:,2))
legend(["Real"  "Sparse EDMD"])

error(1) = norm(H_test - Ypred, 'fro') / norm(H_test, 'fro');
error(2) = norm(H_test - Ypred2, 'fro') / norm(H_test, 'fro');

fprintf('Normal EDMD forecast error: %f\n', error(1));
fprintf('Sparse EDMD forecast error: %f\n', error(2));
save("Optimization Result.mat")

function J = FO(alpha,Nlift,V,W,Xlift,t_test,voltage_test,X,n,npol_states,usesine_states,npol_inputs,usesine_inputs,num_sines_inputs,num_sines_states,H_test)
[Alift,Blift,Clift] = lasso_regression(Nlift,alpha,V,W,Xlift,X);
Ypred = sim_dymamics(t_test,Alift,Blift,Clift,Nlift,npol_states,usesine_states,npol_inputs,usesine_inputs,num_sines_inputs,num_sines_states,n,voltage_test);
J = norm(H_test - Ypred, 'fro') / norm(H_test, 'fro');
end

function [Alift,Blift,Clift] = lasso_regression(Nlift,alpha,V,W,Xlift,X)
for i = 1 : Nlift
    [B,FitInfo] = lasso(V',W(i,:),'Alpha',alpha.alpha);
    [~,index] = min(FitInfo.MSE);
    Alift(i,:) = B(1:Nlift,index)';
    Blift(i,:) = B(Nlift+1:end,index);
end
Q = Xlift;
Nc = size(X,1);
for i = 1 : Nc
    [B,FitInfo] = lasso(Q',X(i,:),'Alpha',alpha.alpha2);
    [~,index] = min(FitInfo.MSE);
    Clift(i,:) = B(:,index)';
end

if isempty(Blift)
    Blift = zeros(Nlift,1);
end

end

function Ypred = sim_dymamics(t_test,Alift,Blift,Clift,Nlift,npol_states,usesine_states,npol_inputs,usesine_inputs,num_sines_inputs,num_sines_states,n,voltage_test)

Nsim = length(t_test);
xlift = zeros(Nlift,1);
for i = 1: Nsim-1
    ulift = build_theta(voltage_test(i),1,npol_inputs,usesine_inputs,num_sines_inputs);
    xlift(:,i+1) = Alift*xlift(:,i) + Blift*ulift';
    yest2 = Clift*xlift(:,i);
    xlift(:,i) = build_theta(yest2',n,npol_states,usesine_states,num_sines_states)';
end

Ypred = (Clift*xlift)';
end

% function yout = build_theta(yin,nVars,polyorder,usesine,num_sines)
% 
% n = size(yin,1);
% yout = zeros(n,1+nVars+(nVars*(nVars+1)/2)+(nVars*(nVars+1)*(nVars+2)/(2*3))+11);
% 
% ind = 1;
% % poly order 0
% yout(:,ind) = ones(n,1);
% ind = ind+1;
% 
% % poly order 1
% for i=1:nVars
%     yout(:,ind) = yin(:,i);
%     ind = ind+1;
% end
% 
% if(polyorder>=2)
%     % poly order 2
%     for i=1:nVars
%         for j=i:nVars
%             yout(:,ind) = yin(:,i).*yin(:,j);
%             ind = ind+1;
%         end
%     end
% end
% 
% if(polyorder>=3)
%     % poly order 3
%     for i=1:nVars
%         for j=i:nVars
%             for k=j:nVars
%                 yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k);
%                 ind = ind+1;
%             end
%         end
%     end
% end
% 
% if(polyorder>=4)
%     % poly order 4
%     for i=1:nVars
%         for j=i:nVars
%             for k=j:nVars
%                 for l=k:nVars
%                     yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k).*yin(:,l);
%                     ind = ind+1;
%                 end
%             end
%         end
%     end
% end
% 
% if(polyorder>=5)
%     % poly order 5
%     for i=1:nVars
%         for j=i:nVars
%             for k=j:nVars
%                 for l=k:nVars
%                     for m=l:nVars
%                         yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k).*yin(:,l).*yin(:,m);
%                         ind = ind+1;
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% if(usesine)
%     for k=1:num_sines
%         yout = [yout sin(k*yin) cos(k*yin)];
%     end
% end
% end

function Phi = build_theta(X,U)
theta1 = X(2,:)';
theta2 = X(1,:)';
dtheta1 = X(4,:)';
dtheta2 = X(3,:)';
V = U';

% theta1: Posición angular del brazo (rad)
% theta2: Posición angular del péndulo (rad)
% dtheta1: Velocidad angular del brazo (rad/s)
% dtheta2: Velocidad angular del péndulo (rad/s)
% V: Voltaje aplicado al motor (V)

% Inicializar el diccionario Phi como una matriz vacía
Phi = [];

% 1. Términos lineales
Phi = [Phi, theta1, theta2, dtheta1, dtheta2, V];

% 2. Términos polinómicos (hasta segundo orden)
Phi = [Phi, theta1.^2, theta2.^2, dtheta1.^2, dtheta2.^2, V.^2];

% 3. Términos de interacción entre variables
Phi = [Phi, theta1 .* theta2, theta1 .* dtheta1, theta1 .* dtheta2, ...
    theta2 .* dtheta1, theta2 .* dtheta2, dtheta1 .* dtheta2, ...
    theta1 .* V, theta2 .* V, dtheta1 .* V, dtheta2 .* V];

% 4. Términos trigonométricos
Phi = [Phi, sin(theta1), cos(theta1), sin(theta2), cos(theta2), ...
    sin(theta1 + theta2), cos(theta1 + theta2)];

% 5. Términos de interacción más complejos (opcional)
Phi = [Phi, theta1 .* sin(theta2), theta2 .* cos(theta1), ...
    dtheta1 .* sin(theta2), dtheta2 .* cos(theta1)];

% 6. Términos cúbicos (opcional, si se necesita mayor no linealidad)
Phi = [Phi, theta1.^3, theta2.^3, dtheta1.^3, dtheta2.^3, V.^3];

% 7. Términos mixtos de mayor orden (opcional)
Phi = [Phi, theta1.^2 .* theta2, theta1 .* theta2.^2, ...
    theta1 .* dtheta1.^2, theta2 .* dtheta2.^2];
Phi = [Phi, exp(theta1), exp(-theta1), exp(theta2), exp(-theta2), ...
    exp(dtheta1), exp(-dtheta1), exp(dtheta2), exp(-dtheta2), ...
    exp(V), exp(-V)];

% 6. Términos exponenciales combinados
Phi = [Phi, exp(theta1 .* theta2), exp(-theta1 .* theta2), ...
    exp(theta1 .* dtheta1), exp(-theta1 .* dtheta1), ...
    exp(theta2 .* dtheta2), exp(-theta2 .* dtheta2), ...
    exp(V .* theta1), exp(-V .* theta1), ...
    exp(V .* theta2), exp(-V .* theta2)];

% 7. Términos mixtos con exponenciales y polinomios (opcional)
Phi = [Phi, theta1 .* exp(theta2), theta2 .* exp(theta1), ...
    dtheta1 .* exp(dtheta2), dtheta2 .* exp(dtheta1), ...
    V .* exp(theta1), V .* exp(theta2)];
% Asegurarse de que Phi es una matriz donde cada columna es un término del diccionario
% y cada fila corresponde a una observación en el tiempo.
Phi = Phi';
end


