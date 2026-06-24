%% - ***************** SIMULACION MPC TANKS MODEL *********************- %%
%% - ********* AND FAULT ESTIMATION BASED ON KALMAN FILTER ************- %%
%% - **************** JUAN FRANCISCO DURAN SIGUENZA *******************- %%
clear; close all; clc
%% - ********************** SYSTEM PARAMETERS *************************- %%
rng(115123)
a2 = 4.17e-5;  A1 = 0.04; a3 = 1.8082e-5; alpha1 = 0; alpha2 = 0;
A2 = 0.04;      q0 = 3.37102e-5;  g = 9.81;

linear_control = load('mpc_linear.mat');

%% - ******************** CONTROLLER PARAMETERS ***********************- %%
Delta_t = 1;                       %Tiempo de muestreo
Np = 5;                           %Horizonte de predicción
Nc = 5;                            %Horizonte de control
lambda1 = 0.01;                       %Factor de penalización DeltaU1
lambda2 = 0.01;                       %Factor de penalización DeltaU2
lambda3 = 0.01;                       %Factor de penalización DeltaU2

sp1 = 0.01*ones(1,Np);
sp2 = 0.01*ones(1,Np);
N_sim = 9000;                        %Tiempo de simulación
%% - ************************ RESTRICCIONES ***************************- %%
deltaumin1 = -0.01;                %DeltaU1 mínima
deltaumax1 = 0.01;                 %DeltaU1 máxima
deltaumin2 = -0.01;                %DeltaU2 mínima
deltaumax2 = 0.01;                 %DeltaU2 máxima
deltaumin3 = -0.01;                %DeltaU3 mínima
deltaumax3 = 0.01;                 %DeltaU3 máxima

umin1 = 0;                      %U1 mínimo
umax1 = 1;                      %U1 máximo
umin2 = 0;                       %U2 mínimo
umax2 = 1;                        %U2 mínimo
umin3 = 0;                       %U3 mínimo
umax3 = 1;                        %U3 mínimo

ymin1 = 0;                      %Y1 mínimo
ymax1 = 0.2;                       %Y1 máximo
ymin2 = 0;                      %Y2 mínimo
ymax2 = 0.2;                       %Y2 máximo

puntos = Nc;
%% %Modelo de predicción.
load('AJUSTAR_RUTA'  % Ruta original de Windows — ajustar al directorio de datos);
npol = 3;
Ad_cont = Alift;
Bd_cont = Blift;
Cd_cont = Clift(1:2,:);
Dd_cont = zeros(size(Clift,1),size(Blift,2));
%% Dimensiones del proceso
n1 = size(Ad_cont,1);           %número de estados
q = size(Cd_cont,1);            %número de salidas
e = size(Bd_cont,2);            %número de entradas
%% %Modelo aumentado (integrador embebido)
[Am,Bm,Cm,Dm] = mod_aumentado_mimo(Ad_cont,Bd_cont,Cd_cont,Dd_cont,n1,q);
%% %Cálculo de las matrices F y Phi
[F,Phi] = mat_f_phi(Am,Bm,Cm,Np,Nc,q,e);
%% %Adecuaciones de R y Puntos de ajuste para el caso multivariable
BarRs = [sp1;sp2];
BarRs = reshape(BarRs,q*Np,1);            %Vector de consignas
Rbase = eye(e*Nc,e*Nc);
lambda_aux=[lambda1;lambda2;lambda3];
lambda=repmat(lambda_aux,Nc,1);
R=lambda.*Rbase;
%% % Preparación de las matrices de restricciones
res_deltau_max=[deltaumax1;deltaumax2;deltaumax3];
res_deltau_min=[-deltaumin1;-deltaumin2;-deltaumin3];
resu=[umax1;umax2;umax3;-umin1;-umin2;-umin3];
resy=[ymax1;ymax2;-ymin1;-ymin2];

[dim_resu,~]=size(resu);
[dim_y,~]=size(resy);

M1_max=eye(Nc*e);               %Matriz superior I para restriccion DeltaU
M1_min=M1_max*(-1);             %Matriz inferior I para restricción DeltaU
M1=[M1_max;M1_min];

ro2=eye(e);
M2_aux=zeros(puntos*e,Nc*e);
rex=tril(ones(puntos*e,Nc*e),0);
M2_aux=repmat(ro2,puntos,Nc);
M2_max=M2_aux.*rex;
M2_min=-M2_max;
M2=[M2_max;M2_min];

M3=[Phi;-Phi];

M=[M1;M2;M3];
[d1,~]=size(M);

for k=1:puntos
    gamma1_max(k*e-1*e+1:1*k*e,1) = res_deltau_max(:,1);
end
for k=1:puntos
    gamma1_min(k*e-1*e+1:1*k*e,1) = res_deltau_min(:,1);
end
gamma1=[gamma1_max;gamma1_min];

for j=1:Np
    maty(j*q-q+1:j*q,1)=resy(1:q);
end
for j=1:Np
    maty(j*q-q+Np*q+1:j*q+Np*q,1)=resy(q+1:dim_y);
end

%% Cálculo de la matriz H
H = (Phi'*Phi+R);

%% Condiciones iniciales
xm = zeros(q,1);                     %Planta a controlar
Xf = [build_theta(zeros(2,1),zeros(3,1),npol);[0;0]];
zeta_current = Xf;
zeta_prev = zeta_current;

u=zeros(e,1);                       %u(k-1)=0
y=zeros(q,1);
ures=zeros(e,1);

%ALGORITMO ITERATIVO
%*******************************************************************
t = 0:Delta_t:(N_sim*Delta_t);
h1(1:2) = [0;0]; h2(1:2) = [0;0];
u = zeros(3,1);
N = length(t);
uplot(:,1:2) = zeros(3,2);
for kkk=2:N-1
    tic
    if kkk>3000 && kkk<4000
        sp1 = 0.2*ones(1,Np);
        sp2 = 0.2*ones(1,Np);
        BarRs = [sp1;sp2];
        BarRs=reshape(BarRs,q*Np,1);            %Vector de consignas
    elseif kkk>=4000 && kkk<7000
        sp1 = 0.1*ones(1,Np);
        sp2 = 0.1*ones(1,Np);
        BarRs = [sp1;sp2];
        BarRs=reshape(BarRs,q*Np,1);            %Vector de consignas
    end
    for i=1:puntos
        gamma2(i*e-e+1:i*e,:)=-u+resu(1:e);
        gamma2(i*e-e+e*puntos+1:i*e+e*puntos,:)=u+resu(e+1:2*e);
    end
    compy=[-F*(Xf);F*(Xf)];
    gamma3=maty+compy;
    %FORMULACION DE GAMMA TOTAL
    %********************************************************************
    gamma=[gamma1;gamma2;gamma3];
    %********************************************************************
    u = calc_accion_control(e,Nc,H,Phi,BarRs,F,Xf,u,M,gamma);
        
    uplot(:,kkk) = u;
    %******************************************************************
    [h1(kkk+1), h2(kkk+1)] = twoTank_v2(u(1), u(2), u(3), t(kkk), t(kkk+1), h1(kkk), h2(kkk));
    Xf = [build_theta([h1(kkk+1);h2(kkk+1)],u,npol) - build_theta([h1(kkk);h2(kkk)],uplot(:,kkk-1),npol);[h1(kkk+1);h2(kkk+1)]];
    aux(kkk) = toc;
end

%% Gráficas de respuesta
k=0:(N_sim-1);
figure
ejes(1) = subplot(231)
plot(k,h1(1:end-1),k,linear_control.h1(1:end-1),'LineWidth',2);grid on
xlabel('Sampling Instant')
ylabel('Height [m]')
legend(["MPC Koopman" "MPC Kalman"])
title("Tank 1")
ax = gca;
ax.FontSize = 20;
set(ax, 'FontName', 'Times', 'FontSize',20,'FontWeight','normal','defaultAxesTickLabelInterpreter','latex')

ejes(2) = subplot(232)
plot(k,h2(1:end-1),k,linear_control.h2(1:end-1),'LineWidth',2);grid on
xlabel('Sampling Instant')
ylabel('Height [m]')
legend(["MPC Koopman" "MPC Kalman"])
title("Tank 2")
ax = gca;
ax.FontSize = 20;
set(ax, 'FontName', 'Times', 'FontSize',20,'FontWeight','normal','defaultAxesTickLabelInterpreter','latex')

ejes(3) = subplot(233)
plot(k,uplot(1,:),k,linear_control.uplot(1,:),'LineWidth',2);grid on
xlabel('Sampling Instant')
legend(["Koopman" "Kalman"])
title('Valve 1 ')
ax = gca;
ax.FontSize = 20;
set(ax, 'FontName', 'Times', 'FontSize',20,'FontWeight','normal','defaultAxesTickLabelInterpreter','latex')

ejes(4) = subplot(234)
plot(k,uplot(2,:),k,linear_control.uplot(2,:),'LineWidth',2);grid on
xlabel('Sampling Instant')
legend(["Koopman" "Kalman"])
title('Valve 2')
ax = gca;
ax.FontSize = 20;
set(ax, 'FontName', 'Times', 'FontSize',20,'FontWeight','normal','defaultAxesTickLabelInterpreter','latex')

ejes(5) = subplot(235)
plot(k,uplot(3,:),k,linear_control.uplot(3,:),'LineWidth',2);grid on
xlabel('Sampling Instant')
legend(["Koopman" "Kalman"])
title('Valve 3 ')
ax = gca;
ax.FontSize = 20;
set(ax, 'FontName', 'Times', 'FontSize',20,'FontWeight','normal','defaultAxesTickLabelInterpreter','latex')
grid on
linkaxes(ejes,'x');

function [h1, h2] = twoTank_v2(u1, u2, u3, t1, t2, h11, h21)
[~, Hsim] = ode45(@(t,h) tanksModel(t, h, u1, u2, u3), [t1 t2], [h11; h21]);
h1 = Hsim(end,1); h2 = Hsim(end,2) ;% Levels measurements
end

function dh = tanksModel(t, h, u1, u2, u3)

a2 = 4.38e-5;  A1 = 0.04; a3 = 4.601e-5; alpha1 = 0; alpha2 = 0;
A2 = 0.04;      q0 = 6.667e-5;  g = 9.81;

if t>800 && t<=4000
    alpha1 = 0e-6;
    alpha2 = 0e-5;
elseif t>1000 && t<=1500
    alpha1 = 0;
    alpha2 = 0;
elseif t>1500 && t<1900
    alpha1 = 0;
    alpha2 = 0;
else t>=3000;
    alpha1 = 2e-5;
    alpha2 = 0;
end

h1 = h(1);                                  % Tank 1 level
h2 = h(2);                                  % Tank 2 level
dh1 = 1/A1 * (u1*q0 - u2*a2*sqrt(2*g*h1)) - alpha1*sqrt(2*g*h1)/A1;% SE h1
dh2 = 1/A2 * (u2*a2*sqrt(2*g*h1) - u3*a3*sqrt(2*g*h2) - alpha2*sqrt(2*g*h2)); % SE for h2
dh = real([dh1; dh2]);
end

function theta = build_theta(X,U,n)

X =  X';
U = U';
theta = X;
k = 1;
for i = 0:n
    for j = 0: n
        if (j==1 && i == 0) || (j==0 && i == 1)
            k = k;
        else
            theta(:,k) = real(X(:,1).^(i).*X(:,2).^(j));
            k = k + 1;
        end
    end
end

theta = theta';
end

