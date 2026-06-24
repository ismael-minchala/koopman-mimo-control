%% - ***************** SIMULACION MPC TANKS MODEL *********************- %%
%% - ********* AND FAULT ESTIMATION BASED ON KALMAN FILTER ************- %%
%% - **************** JUAN FRANCISCO DURAN SIGUENZA *******************- %%
clear; clc; close all
%% - ********************** SYSTEM PARAMETERS *************************- %%
% a2 = 4.38e-5;  A1 = 0.04; a3 = 4.601e-5; alpha1 = 0; alpha2 = 0;
% A2 = 0.04;      q0 = 6.667e-5;  g = 9.81;

a2 = 4.17e-5;  A1 = 0.04; a3 = 1.8082e-5; alpha1 = 0; alpha2 = 0;
A2 = 0.04;      q0 = 3.37102e-5;  g = 9.81;

he1 = 0.05; he2 = 0.05;
ue1 = 0.9; ue2 = 1; ue3 = 1;
alphae1 = 0; alphae2 = 0;

%% - ******************** CONTROLLER PARAMETERS ***********************- %%
Delta_t = 1;                       %Tiempo de muestreo
% Np = 30;                           %Horizonte de predicción
% Nc = 5;                            %Horizonte de control
% lambda1 = 4;                      %Factor de penalización DeltaU1
% lambda2 = 4;                       %Factor de penalización DeltaU2
% lambda3 = 4;                       %Factor de penalización DeltaU2


Np = 30;                           %Horizonte de predicción
Nc = 5;                            %Horizonte de control
lambda1 = 1;                      %Factor de penalización DeltaU1
lambda2 = 1;                       %Factor de penalización DeltaU2
lambda3 = 1;                       %Factor de penalización DeltaU2


sp1 = 0.05*ones(1,Np);
sp2 = 0.05*ones(1,Np);
N_sim = 9000;                        %Tiempo de simulación
%% - ************************ RESTRICCIONES ***************************- %%
deltaumin1 = -0.1;                %DeltaU1 mínima
deltaumax1 = 0.1;                 %DeltaU1 máxima
deltaumin2 = -1;                %DeltaU2 mínima
deltaumax2 = 1;                 %DeltaU2 máxima
deltaumin3 = -1;                %DeltaU3 mínima
deltaumax3 = 1;                 %DeltaU3 máxima

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
Acont = [-(sqrt(2)*g*(a2*ue2+alphae1))/(2*sqrt(he1*g)*A1) 0 (-sqrt(2*g*he1))/A1 0;(sqrt(2)*g*(a2*ue2+alphae1))/(2*sqrt(he1*g)*A1) -(sqrt(2)*g*(a3*ue3+alphae2))/(2*sqrt(he2*g)*A2) 0 (-sqrt(2*g*he2))/A2;zeros(2,4)];
Bcont = [q0/A1 -(sqrt(2*g*he1)*a2)/A1 0;0 (sqrt(2*g*he1)*a2)/A2 -(sqrt(2*g*he2)*a2)/A2;zeros(2,3)];
Ccont = [1 0 0 0;0 1 0 0];
Dcont = zeros(2,3);

[Ad_cont,Bd_cont,Cd_cont,Dd_cont] = c2dm(Acont,Bcont,Ccont,Dcont,Delta_t);
%% %Dimensiones del proceso
n1=size(Ad_cont,1);           %número de estados
q=size(Cd_cont,1);            %número de salidas
e=size(Bd_cont,2);            %número de entradas

xe(:,1:2) = [0 0;0 0;0 0;0 0];
xe(:,1:2) = xe(:,1:2) - [he1;he2;0;0];
sigma = 3.6e-4;

Qf = sigma^3*eye(2);
Q = sigma*eye(2);
Q = [Q zeros(2,2);zeros(2,2) Qf];
Rkf = 0.5*eye(2);
P = 1e-10.*eye(n1);
Hkf = Cd_cont;

%% %Modelo aumentado (integrador embebido)
[Am,Bm,Cm,Dm]=mod_aumentado_mimo(Ad_cont,Bd_cont,Cd_cont,Dd_cont,n1,q);
%% %Cálculo de las matrices F y Phi
[F,Phi]=mat_f_phi(Am,Bm,Cm,Np,Nc,q,e);
%% %Adecuaciones de R y Puntos de ajuste para el caso multivariable
BarRs=[sp1;sp2];
BarRs=reshape(BarRs,q*Np,1);            %Vector de consignas
Rbase=eye(e*Nc,e*Nc);
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
    gamma1_max(k*e-1*e+1:1*k*e,1)=res_deltau_max(:,1);
end
for k=1:puntos
    gamma1_min(k*e-1*e+1:1*k*e,1)=res_deltau_min(:,1);
end
gamma1=[gamma1_max;gamma1_min];

for j=1:Np
    maty(j*q-q+1:j*q,1)=resy(1:q);
end
for j=1:Np
    maty(j*q-q+Np*q+1:j*q+Np*q,1)=resy(q+1:dim_y);
end

%% Cálculo de la matriz H
H=(Phi'*Phi+R);

%% Condiciones iniciales
xm=zeros(n1,1);                     %Planta a controlar
Xf=zeros(n1+q,1);                   %Modelo aumentado
u=zeros(e,1);                       %u(k-1)=0
y=zeros(q,1);
ures=zeros(e,1);

%ALGORITMO ITERATIVO
%*******************************************************************
t = 0:Delta_t:(N_sim*Delta_t);
h1(1:2) = [0;0]; h2(1:2) = [0;0];
N = length(t);
Du(:,1:2) = zeros(e,2) - [ue1;ue2;ue3];
y(:,1:2) = zeros(2,2);
f(:,1:2) = zeros(2,2);
for kkk=2:N-1
    tic
    Xf = [[h1(kkk);h2(kkk);f(:,kkk)] - ([h1(kkk-1);h2(kkk-1);f(:,kkk-1)]);[h1(kkk);h2(kkk)]];
    y(:,kkk) =[h1(kkk);h2(kkk)] ;
    for i=1:puntos
        gamma2(i*e-e+1:i*e,:)=-u+resu(1:e);
        gamma2(i*e-e+e*puntos+1:i*e+e*puntos,:)=u+resu(e+1:2*e);
    end
    compy=[-F*Xf;F*Xf];
    gamma3=maty+compy;
    %FORMULACION DE GAMMA TOTAL
    %********************************************************************
    gamma=[gamma1;gamma2;gamma3];
    %********************************************************************
    u=calc_accion_control(e,Nc,H,Phi,BarRs,F,Xf,u,M,gamma);
    
    Du(:,kkk) = u - [ue1;ue2;ue3];

    if u(2) < 0 
        u(2) = 0;
    end

    if u(3) < 0 
        u(3) = 0;
    end
    
    uplot(:,kkk) = u;

    P = Ad_cont*P*Ad_cont' + Q;
    K = P*Hkf'*(Hkf*P*Hkf'+Rkf)^(-1);
    xe(:,kkk+1) = Ad_cont*xe(:,kkk) + Bd_cont*Du(:,kkk) + K*(y(:,kkk) - (Hkf*xe(:,kkk) + [he1;he2]));
    ye(:,kkk) = Hkf*xe(:,kkk) + [he1;he2];
    P = (eye(n1)-K*Hkf)*P;
    f(:,kkk+1) = [zeros(2,2) eye(2)]*xe(:,kkk);
    r(:,kkk) = y(:,kkk) - (Hkf*xe(:,kkk) + [he1;he2]);
    %******************************************************************
    [h1(kkk+1), h2(kkk+1)] = twoTank_v2(u(1), u(2), u(3), t(kkk), t(kkk+1), h1(kkk), h2(kkk));
    h1(kkk+1) =  h1(kkk+1) + sigma*randn(1,1);
    h2(kkk+1) =  h2(kkk+1) + sigma*randn(1,1);
    aux(kkk) = toc;
end

%% Gráficas de respuesta
k=0:(N_sim-1);
figure(3)
subplot(231)
plot(k,h1(1:end-1),'LineWidth',2);grid on
xlabel('Sampling Instant')
ylabel('Height [m]')
title("Tank 1")
ax = gca;
ax.FontSize = 20;
set(ax, 'FontName', 'Times', 'FontSize',20,'FontWeight','normal','defaultAxesTickLabelInterpreter','latex')

subplot(232)
plot(k,h2(1:end-1),'LineWidth',2);grid on
xlabel('Sampling Instant')
ylabel('Height [m]')
title("Tank 2")
ax = gca;
ax.FontSize = 20;
set(ax, 'FontName', 'Times', 'FontSize',20,'FontWeight','normal','defaultAxesTickLabelInterpreter','latex')

subplot(233)
stairs(k,uplot(1,:),'LineWidth',2);grid on
xlabel('Sampling Instant')
ylabel('Aperture Coefficient')
title('Valve 1 ')
ax = gca;
ax.FontSize = 20;
set(ax, 'FontName', 'Times', 'FontSize',20,'FontWeight','normal','defaultAxesTickLabelInterpreter','latex')

subplot(234)
stairs(k,uplot(2,:),'LineWidth',2);grid on
xlabel('Sampling Instant')
ylabel('Aperture Coefficient')
title('Valve 2')
ax = gca;
ax.FontSize = 20;
set(ax, 'FontName', 'Times', 'FontSize',20,'FontWeight','normal','defaultAxesTickLabelInterpreter','latex')

subplot(235)
stairs(k,uplot(3,:),'LineWidth',2);grid on
xlabel('Sampling Instant')
ylabel('Aperture Coefficient')
title('Valve 3')
ax = gca;
ax.FontSize = 20;
set(ax, 'FontName', 'Times', 'FontSize',20,'FontWeight','normal','defaultAxesTickLabelInterpreter','latex')

subplot(236)
stairs(k,f(:,1:end-1)','LineWidth',2);grid on
xlabel('Sampling Instant')
ylabel('Coefficient')
title('Leak Faults')
legend('\alpha1','\alpha2')
ax = gca;
ax.FontSize = 20;
set(ax, 'FontName', 'Times', 'FontSize',20,'FontWeight','normal','defaultAxesTickLabelInterpreter','latex')

figure(4)
stairs(k,r','LineWidth',2);grid on
xlabel('Sampling Instant')
ylabel('Error [m]')
title('Generated Residual')
legend('\alpha1','\alpha2')
ax = gca;
ax.FontSize = 20;
set(ax, 'FontName', 'Times', 'FontSize',20,'FontWeight','normal','defaultAxesTickLabelInterpreter','latex')


function [h1, h2] = twoTank_v2(u1, u2, u3, t1, t2, h11, h21)
[~, Hsim] = ode45(@(t,h) tanksModel(t, h, u1, u2, u3), [t1 t2], [h11; h21]);
h1 = Hsim(end,1); h2 = Hsim(end,2) ;% Levels measurements
end

function dh = tanksModel(t, h, u1, u2, u3)


a2 = 4.38e-5;  A1 = 0.04; a3 = 4.601e-5; alpha1 = 0; alpha2 = 0;
A2 = 0.04;      q0 = 6.667e-5;  g = 9.81;

if t>1500 && t<=2000
    alpha1 = 1.5e-5;
    alpha2 = 0;
elseif t>2000 && t<=4000
    alpha1 = 0;
    alpha2 = 2.2e-5;
elseif t>4000 && t<6000
    alpha1 = 0;
    alpha2 = 0e-5;
elseif t>=6000
    alpha1 = 1e-5;
    alpha2 = 2e-5;
end


h1 = h(1);                                  % Tank 1 level
h2 = h(2);                                  % Tank 2 level
dh1 = 1/A1 * (u1*q0 - u2*a2*sqrt(2*g*h1)) - alpha1*sqrt(2*g*h1)/A1;% SE h1
dh2 = 1/A2 * (u2*a2*sqrt(2*g*h1) - u3*a3*sqrt(2*g*h2) - alpha2*sqrt(2*g*h2)); % SE for h2
dh = real([dh1; dh2]);
end

