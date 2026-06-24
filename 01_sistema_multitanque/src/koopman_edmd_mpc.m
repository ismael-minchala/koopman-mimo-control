%% - ******************** -Koopman approximation- ********************** %%
%% - ********************** -Multitank System- ************************* %%
%% - ******************** -Juan Francisco Duran S- ********************* %%
clear; close all; clc
%% - ********************** -Define Variables- ************************* %%
rng(115123)
n = 2;                              %Number of states
m = 3;                              %Number of control inputs
deltaT = 1;
a2 = 4.17e-5;  A1 = 0.04; a3 = 1.8082e-5; alpha1 = 0; alpha2 = 0;
A2 = 0.04;      q0 = 3.37102e-5;  g = 9.81;

he1 = 0.0879; he2 = 0.1153;
ue1 = 0.2686; ue2 = 0.1947; ue3 = 0.32;
alphae1 = 0; alphae2 = 0;


Acont = [-(sqrt(2)*g*(a2*ue2+alphae1))/(2*sqrt(he1*g)*A1) 0;(sqrt(2)*g*(a2*ue2+alphae1))/(2*sqrt(he1*g)*A1) -(sqrt(2)*g*(a3*ue3+alphae2))/(2*sqrt(he2*g)*A2)];
Bcont = [q0/A1 -(sqrt(2*g*he1)*a2)/A1 0;0 (sqrt(2*g*he1)*a2)/A2 -(sqrt(2*g*he2)*a2)/A2];
Ccont = [1 0;0 1];
Dcont = zeros(2,3);

sys_linc = ss(Acont,Bcont,Ccont,Dcont);
sys_lind = c2d(sys_linc,deltaT);
%% - ************************* -Collect Data- ************************** %%
disp('Starting data collection \n')
Ntraj = 1; Nsim = 300000;
Cy = [1 0;0 1];                     % Output matrix                                 % Number of delays1
ny = size(Cy,1);                    % Number of outputs

load(fullfile('..','datos','Dataset_13.mat'))

%load('AJUSTAR_RUTA'  % Ruta original de Windows — ajustar al directorio de datos)

U = [u(:,1) u(:,2) u(:,3)]';
t = 0:deltaT:length(U)-1;

fprintf('Data Collection Done \n')
figure(1)
ejes1(1) = subplot(211);
plot(t,h,'LineWidth',2)
title('Training Set','FontName', 'Times', 'FontSize',15,'FontWeight','normal')
xlabel('Time [s]')
ylabel('Heights [m]')
legend(["$h_1$"  "$h_2$"])
ax = gca;
set(ax,'FontName', 'Times', 'FontSize',15,'FontWeight','normal')
grid on
ejes1(2) = subplot(212);
plot(t,U','LineWidth',2)
title('Training Set','FontName', 'Times', 'FontSize',15,'FontWeight','normal')
xlabel('Time [s]')
ylabel('Control Inputs')
legend(["$\mu_1$"  "$\mu_2$" "$\mu_3$"])
ax = gca;
set(ax,'FontName', 'Times', 'FontSize',15,'FontWeight','normal')
grid on
set(groot,...
'defaulttextinterpreter','latex',...
'defaultAxesTickLabelInterpreter','latex',...
'defaultLegendInterpreter','latex')

linkaxes(ejes1,'x')

Ntest = 10e3;
ttest = 0:deltaT:(Ntest-1)*deltaT;
h1_ref = 0.2*input_generator(Ntest,2500);
h2_ref = 0.2*input_generator(Ntest,3000);
sigma  = 0;
%% - ************************ -BASIS FUNCTION- ************************* %%
X = h(1:end-1,:)' + randn(2,Nsim)*sigma;
Y = h(2:end,:)' + randn(2,Nsim)*sigma;
U = U(:,1:end-1);
npol = 3;

%% - ************************ -EDMD ALGORITHM- ************************* %%
fprintf('Starting Lifting \n')

Xlift = build_theta(X,npol);
Ylift = build_theta(Y,npol);
XU_lift = build_theta_cross(X,U);

Nlift = size(Xlift,1);
Crosslift = size(XU_lift,1);

W = [Ylift;XU_lift;zeros(3,300000)];
V = [Xlift;XU_lift;U];
VVt = V*V';
WVt = W*V';
ABC = WVt * pinv(VVt);
Alift = ABC(1:Nlift+Crosslift,1:Nlift+Crosslift);
Blift = ABC(1:Nlift+Crosslift,Nlift+Crosslift+1:end);
Clift = X*pinv([Xlift;XU_lift]);
fprintf( 'DMD Regression residual %f \n', norm([Ylift;XU_lift] - Alift*[Xlift;XU_lift] - Blift*U,'fro') / norm([Ylift;XU_lift],'fro') );
%% - ************************ -EDMD ALGORITHM 2- ************************* %%
fprintf('Starting Lifting \n')

Xlift = build_theta(X,npol);
Ylift = build_theta(Y,npol);

Nlift = size(Xlift,1);
Crosslift = size(XU_lift,1);

W = [Ylift];
V = [Xlift;U];
VVt = V*V';
WVt = W*V';
ABC = WVt * pinv(VVt);
Alift2 = ABC(1:Nlift,1:Nlift);
Blift2 = ABC(1:Nlift,Nlift+1:end);
Clift2 = X*pinv([Xlift]);
sys = ss(Alift,Blift,Clift,0,deltaT);
%% - ************************* -Comparison- **************************** %%
[V1,D1] = eig(Alift);
[V2,D2] = eig(Alift2);
D1 = diag(D1);
D2 = diag(D2);

zoomd1 = D1(imag(D1) == 0);
zoomd2 = D2(imag(D2) == 0);

tmax = 5000;
Nsim = tmax/deltaT;

uprbs1 = [input_generator(Nsim,5) input_generator(Nsim,1) input_generator(Nsim,10)];
%uprbs1 = [0.27*ones(tmax,1) 0.15*ones(tmax,1) 0.27*ones(tmax,1)];
Xnext3 = zeros(2,1);
xstart = zeros(2,1);
%
xlift = zeros(Nlift+Crosslift,1);
xlift2 = zeros(Nlift,1);
ylift = Clift*xlift;
t = 1:deltaT:Nsim;
for i = 1: Nsim-1
    [Xnext3(1,i+1),Xnext3(2,i+1)] = twoTank_v2(uprbs1(i,1),uprbs1(i,2),uprbs1(i,3),t(i),t(i+1),Xnext3(1,i),Xnext3(2,i),Ntraj);
    yest = Clift*xlift(:,i)+sigma*randn(2,1);
    xlift(:,i) = [build_theta(yest,npol);build_theta_cross(yest,uprbs1(i,:)')];
    xlift(:,i+1) = Alift*xlift(:,i) + Blift*uprbs1(i,:)';
    xlift2(:,i+1) = Alift2*xlift2(:,i) + Blift2*uprbs1(i,:)';
end
D3 = complex(eig(sys_lind.A));
%%
fprintf( 'Prediction residual %f \n', norm(Xnext3 - Clift*xlift,'fro') / norm(Xnext3,'fro') );

figure(2)
yest = Clift*xlift;
yest2 = Clift2*xlift2;
yreal = Cy*Xnext3;
subplot(211)
plot([0:Nsim-1]*deltaT,yreal(1,:),'-r','LineWidth',2); hold on
plot([0:Nsim-1]*deltaT,yest(1,:),':b','LineWidth',2); hold on
plot([0:Nsim-1]*deltaT,yest2(1,:),'-g','LineWidth',2); hold on
xlabel('Time [s]')
ylabel('Height in tank 1 [m]')
legend(["True"  "Koopman Estimation $\Psi(\mathbf{x},\;\mathbf{u})$" "Koopman Estimation $\Psi(\mathbf{x})$"])
ax = gca;
set(ax,'FontName', 'Times', 'FontSize',15,'FontWeight','normal')
grid on

subplot(212)
plot([0:Nsim-1]*deltaT,yreal(2,:),'-r','LineWidth',2); hold on
plot([0:Nsim-1]*deltaT,yest(2,:),':b','LineWidth',2); hold on
plot([0:Nsim-1]*deltaT,yest2(2,:),'-g','LineWidth',2); hold on
xlabel('Time [s]')
ylabel('Height in tank 2 [m]')
% legend(["True"  "Koopman Estimation $\Psi(\mathbf{x},\;\mathbf{u})$" "Koopman Estimation $\Psi(\mathbf{x})$"])
ax = gca;
set(ax,'FontName', 'Times', 'FontSize',15,'FontWeight','normal')
grid on
sgtitle('Comparison between actual and forecast dynamics','FontName', 'Times', 'FontSize',20,'FontWeight','normal')


figure(3)
plot([0:Nsim-1]*deltaT,uprbs1,'LineWidth',2); hold on
xlabel("Time [s]")
title("Control signals of the validation process")
legend(["$\mu_1$" "$\mu_2$" "$\mu_3$"])
ax = gca;
set(ax,'FontName', 'Times', 'FontSize',15,'FontWeight','normal')
grid on
set(groot,...
'defaulttextinterpreter','latex',...
'defaultAxesTickLabelInterpreter','latex',...
'defaultLegendInterpreter','latex')

figure(4);

plot(D1,'*','MarkerSize',20,'Color','b'); % Azul para D1
hold on;
plot(D2,'*','MarkerSize',20,'Color','r'); % Rojo para D2
plot(D3,'*','MarkerSize',20,'Color','black')
legend(["with $\Psi(\mathbf{x},\;\mathbf{u})$", "with $\Psi(\mathbf{x})$" "Real poles"])

xlabel('Real Axis', 'Interpreter', 'latex');
ylabel('Imaginary Axis', 'Interpreter', 'latex');
title('Estimated Koopman Eigenvalues', 'Interpreter', 'latex');

% Configuración del eje
ax = gca;
set(ax, 'FontName', 'Times', 'FontSize', 15, 'FontWeight', 'normal');

% Configuración global para interpretación LaTeX
set(groot, ...
    'defaultTextInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex');

grid on;


%% - ******************* -MPC Koopman controller- ********************** %%
linear_control = load(fullfile('..','datos','linear_control_Results.mat'));
%linear_control = load('AJUSTAR_RUTA'  % Ruta original de Windows — ajustar al directorio de datos);


Tmax = 5000;
Nsim = Tmax/deltaT;

ymin = [0;0];
ymax = [0.2;0.2];
x0 = [0;0];
%yrr = [0.1*ones(1,Nsim);0.1*ones(1,Nsim)];
%yrr = [[0.15*ones(1,floor(Nsim/2)) 0.05*ones(1,floor(Nsim/2))];[0.2*ones(1,floor(Nsim/2)) 0.1*ones(1,floor(Nsim/2))]];
yrr = [0.1*sin(2*pi*[1:Nsim] / Nsim);0.1*sin(2*pi*[1:Nsim] / Nsim)]+0.1;
yrr_bound_up = [(0.1+0.002)*sin(2*pi*[1:Nsim] / Nsim); (0.1 + 0.002)*sin(2*pi*[1:Nsim] / Nsim)]+0.1;
yrr_bound_down = [(0.1-0.002)*sin(2*pi*[1:Nsim] / Nsim); (0.1-0.002)*sin(2*pi*[1:Nsim] / Nsim)]+0.1;

Q = 10*eye(2);
R = [1e-5;1e-5;1e-5].*eye(3);

Tpred = 5;
Np = round(Tpred/deltaT);
% xlift_min = nan(Nlift+Crosslift,1);
% xlift_max = nan(Nlift+Crosslift,1);
% xlift_min([2 npol+2]) = 0;
% xlift_max([2 npol+2]) = 0.2;
xlift_min = [];
xlift_max = [];

koopmanMPC = getMPC(Alift,Blift,Clift,0,Q,R,Q,Np,[0;0;0],[1;1;1],xlift_min,xlift_max,'qpoases');
linMPC = getMPC(sys_lind.A,sys_lind.B,eye(2),0,Q,R,Q,Np,[0;0;0],[1;1;1],[0;0],[0.2;0.2],'qpoases');

u_koop = [0;0;0];
x_koop = x0;
x_lin = x0;

x_est = [build_theta(x_koop,npol);build_theta_cross(x_koop,u_koop)];
x_est2 = x0;
xlift = x_est;
XX_koop = x0; UU_koop = [];
XX_lin = x0; UU_lin = [];
z = 1;
t = 0:deltaT:Nsim-1;
for i = 1:Nsim-1
    yr = yrr(:,i+1);
    u_koop = koopmanMPC(xlift,yr);
    u_lin = linMPC(x_lin,yr);
    xlift = [build_theta(x_koop,npol);build_theta_cross(x_koop,u_koop)];
    x_est(:,i) = [build_theta(Clift*x_est(:,i),npol);build_theta_cross(Clift*x_est(:,i),u_koop)];
    x_est(:,i+1) = Alift*x_est(:,i) + Blift*u_koop;
    x_est2(:,i+1) = sys_lind.A*x_est2(:,i) + sys_lind.B*u_koop;
    [x_koop(1),x_koop(2)] = twoTank_v2(u_koop(1), u_koop(2), u_koop(3), t(i), t(i+1), x_koop(1), x_koop(2),1);
    [x_lin(1),x_lin(2)] = twoTank_v2(u_lin(1), u_lin(2), u_lin(3), t(i), t(i+1), x_lin(1), x_lin(2),1);
    XX_koop = [XX_koop x_koop];
    UU_koop = [UU_koop u_koop];
    XX_lin = [XX_lin x_lin];
    UU_lin = [UU_lin u_lin];
    e(:,i) = (x_koop - yr).^2;
    e2(:,i) = (x_lin - yr).^2;
    if i>1
        ISE(:,i) = ISE(:,i-1) + (deltaT/2)*(e(i)+e(i-1));
        ISE2(:,i) = ISE2(:,i-1) + (deltaT/2)*(e2(i)+e2(i-1));

    else
        ISE(:,i) = zeros(2,1);
        ISE2(:,i) = zeros(2,1);
    end

end

%%
figure(5)
ejes(1) = subplot(211)
plot([0:Nsim-1]*deltaT,XX_koop(1,:),'b-',[0:Nsim]*deltaT,linear_control.h1,'r-',[0:Nsim-1]*deltaT,yrr(1,:),'g--','LineWidth',2); hold on
legend(["$H_1$ Koopman" "$H_1$ Linear" "Reference"])
xlabel('Time [s]')
ylabel('Height in tank 1 [m]')
ax = gca;
set(ax,'FontName', 'Times', 'FontSize',15,'FontWeight','normal')
grid on
ejes(2) = subplot(212)
plot([0:Nsim-1]*deltaT,XX_koop(2,:),'b-',[0:Nsim]*deltaT,linear_control.h2,'r-',[0:Nsim-1]*deltaT,yrr(2,:),'g--','LineWidth',2); hold on
legend(["$H_2$ Koopman" "$H_2$ Linear" "Reference"])
xlabel('Time [s]')
ylabel('Height in tank 2 [m]')
ax = gca;
set(ax,'FontName', 'Times', 'FontSize',15,'FontWeight','normal')
sgtitle("Comparison between Linear MPC and Koopman MPC",'FontName', 'Times', 'FontSize',15,'FontWeight','normal')
grid on

set(groot, ...
    'defaultTextInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex');

linkaxes(ejes,'x')

y_est = Clift*x_est;
figure(6)
subplot(211)
plot([0:Nsim-1]*deltaT,XX_koop(1,:), '--','LineWidth',2);hold on
plot([0:Nsim-1]*deltaT,y_est(1,:),'LineWidth',2); 
plot([0:Nsim-1]*deltaT,x_est2(1,:),'LineWidth',2);
ylabel('Height in tank 1')
legend(["Closed Loop Plant Dynamic" "Koopman Prediction" "Linear Prediction"]);
ax = gca;
set(ax,'FontName', 'Times', 'FontSize',15,'FontWeight','normal')
grid on

subplot(212)
plot([0:Nsim-1]*deltaT,XX_koop(2,:), '--','LineWidth',2);hold on
plot([0:Nsim-1]*deltaT,y_est(2,:),'LineWidth',2); 
plot([0:Nsim-1]*deltaT,x_est2(2,:),'LineWidth',2);
ylabel('Height in tank 2')
ax = gca;
set(ax,'FontName', 'Times', 'FontSize',15,'FontWeight','normal')
sgtitle("Model Prediction in Closed Loop Operation",'FontName', 'Times', 'FontSize',15,'FontWeight','normal')
grid on

figure(7)

plot([0:Nsim-1]*deltaT,XX_koop(1,1:end) - y_est(1,:), [0:Nsim-1]*deltaT,XX_koop(2,1:end) - y_est(2,:),'LineWidth',2); hold on
title("Generated Residual from Actual and Estimated Outputs",'FontName', 'Times', 'FontSize',15,'FontWeight','normal')
xlabel('Time [s]')
ylabel('Height Error [m]')
legend(["Tank 1", "Tank 2"])
ax = gca;
set(ax,'FontName', 'Times', 'FontSize',15,'FontWeight','normal')
grid on

figure(8)
model_size = [0,1,2,3,4,5,10,20,50,100];
model_error = [1, 0.722, 0.0703, 0.064167, 0.049, 0.087, 0.093, 0.0609, 0.0609, 0.0609];
plot(model_size,model_error,'*','MarkerSize',8); hold on
title("Error evolution with respect to the dimension of the A matrix of the model",'FontName', 'Times', 'FontSize',15,'FontWeight','normal')
xlabel('Polynomial order [n]')
ylabel('Frobenius norm of Height Error')
ax = gca;
set(ax,'FontName', 'Times', 'FontSize',15,'FontWeight','normal')
grid on

set(groot, ...
    'defaultTextInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex');

save("Koopman_model_multitank.mat",sys);


function [h1, h2] = twoTank_v2(u1, u2, u3, t1, t2, h11, h21,Ntraj)
[~, Hsim] = ode45(@(t,h) tanksModel(t, h, u1, u2, u3), [t1 t2], [h11; h21]);
h1 = Hsim(end,1:Ntraj); h2 = Hsim(end,Ntraj+1:end);
end

function dh = tanksModel(t, h, u1, u2, u3)
a2 = 4.17e-5;  A1 = 0.04; a3 = 1.8082e-5; alpha1 = 0; alpha2 = 0;
A2 = 0.04;      q0 = 3.37102e-5;  g = 9.81;


% if t>400 && t<=1200
%     alpha1 = 0e-5;
%     alpha2 = 0e-5;
% elseif t>1200 && t<=2000
%     alpha1 = 0e-5;
%     alpha2 = 0e-5;
% elseif t>2000 && t<2500
%     alpha1 = 0e-5;
%     alpha2 = 0e-5;
% elseif t > 2500 && t<4000
%     alpha1 = 0e-5;
%     alpha2 = 0.2e-5;
% else
% if t >= 3000
%     alpha1 = 1e-5;
%     alpha2 = 0e-5;
% end


h1 = h(1);                                  % Tank 1 level
h2 = h(2);                                  % Tank 2 level
dh(1,:) = real(1/A1 * (u1*q0 - u2*a2*sqrt(2*g*h1)) - alpha1*sqrt(2*g*h1)/A1);% SE h1
dh(2,:) = real(1/A2 * (u2*a2*sqrt(2*g*h1) - u3*a3*sqrt(2*g*h2) - alpha2*sqrt(2*g*h2))); % SE for h2
dh = reshape(dh,[],1);

end

function u = input_generator(N,sampling_factor)

u = idinput([N,1,1],'rgs',[0 1e-2],[0,1]);
sampled_indices = 1:sampling_factor:length(u);
sampled_values = u(sampled_indices);
u = repelem(sampled_values,sampling_factor);
u = u(1:N);
u = min(max(u,0),1);

end

function theta = build_theta(X,n)

X =  X';
k = 1;
for i = 0:n
    for j = 0: n
        theta(:,k) = real(X(:,1).^(i).*X(:,2).^(j));
        k = k + 1;
    end
end
theta = theta';

end



function theta = build_theta_cross(X,U)

X =  X';
U = U';
k = 1;

inputs = size(U,2);
states = size(X,2);

for i = 1:inputs
    for j = 1:states
    result(:, k) = (X(:, j)) .* U(:, i);
    k = k + 1;
    end
end


theta = (result');

end

