%% - ***************** SIMULACION MPC TANKS MODEL *********************- %%
%% - ********* AND FAULT ESTIMATION BASED ON KALMAN FILTER ************- %%
%% - **************** JUAN FRANCISCO DURAN SIGUENZA *******************- %%
clear; close all; clc
%% - ********************** SYSTEM PARAMETERS *************************- %%
a2 = 4.38e-5;  A1 = 0.04; a3 = 4.601e-5; alpha1 = 0; alpha2 = 0;
A2 = 0.04;      q0 = 6.667e-5;  g = 9.81;

he1 = 0.0879; he2 = 0.1153;
ue1 = 0.2686; ue2 = 0.1947; ue3 = 0.32;
alphae1 = 0; alphae2 = 0;


sigma = 0;
%% - ******************** CONTROLLER PARAMETERS ***********************- %%
Delta_t = 1;                       %Tiempo de muestreo
Np = 30;                           %Horizonte de predicción
Nc = 5;                            %Horizonte de control
lambda1 = 4;                      %Factor de penalización DeltaU1
lambda2 = 2;                       %Factor de penalización DeltaU2
lambda3 = 2;                       %Factor de penalización DeltaU2

sp1 = 0.0879*ones(1,Np);
sp2 = 0.11*ones(1,Np);
N_sim = 7000;                        %Tiempo de simulación
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
%% - ******************** MODELO DE PREDICCION ************************- %%
Acont = [-(sqrt(2)*g*(a2*ue2+alphae1))/(2*sqrt(he1*g)*A1) 0;(sqrt(2)*g*(a2*ue2+alphae1))/(2*sqrt(he1*g)*A1) -(sqrt(2)*g*(a3*ue3+alphae2))/(2*sqrt(he2*g)*A2)];
Bcont = [q0/A1 -(sqrt(2*g*he1)*a2)/A1 0;0 (sqrt(2*g*he1)*a2)/A2 -(sqrt(2*g*he2)*a2)/A2];
Ccont = [1 0;0 1];
Dcont = zeros(2,3);
Bf = [-34.8839 0;-0.1336 -33.1908];

[Ad_cont,Bd_cont,Cd_cont,Dd_cont] = c2dm(Acont,Bcont,Ccont,Dcont,Delta_t);
%% - ****************** DIMENSIONES DEL PROCESO ***********************- %%
n1 = size(Ad_cont,1);           %número de estados
q = size(Cd_cont,1);            %número de salidas
e = size(Bd_cont,2);            %número de entradas
%% - ************ MODELO AUMENTADO - iNTEGRADOR EMBEBIDO **************- %%
[Am,Bm,Cm,Dm]=mod_aumentado_mimo(Ad_cont,Bd_cont,Cd_cont,Dd_cont,n1,q);
%% - **************** CALCULO DE LAS MATRICES F Y PHI *****************- %%
[F,Phi]=mat_f_phi(Am,Bm,Cm,Np,Nc,q,e);
%% -**Adecuaciones de R y Puntos de ajuste para el caso multivariable**- %%
BarRs=[sp1;sp2];
BarRs=reshape(BarRs,q*Np,1);            %Vector de consignas
Rbase=eye(e*Nc,e*Nc);
lambda_aux=[lambda1;lambda2;lambda3];
lambda=repmat(lambda_aux,Nc,1);
R=lambda.*Rbase;
%% - ********* Preparación de las matrices de restricciones ***********- %%
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
rex=tril(ones(puntos*e,Nc*e),0);
M2_aux=repmat(ro2,puntos,Nc);
M2_max=M2_aux.*rex;
M2_min=-M2_max;
M2=[M2_max;M2_min];

M3=[Phi;-Phi];

M=[M1;M2;M3];
[d1,~]=size(M);

for t=1:puntos
    gamma1_max(t*e-1*e+1:1*t*e,1)=res_deltau_max(:,1);
end
for t=1:puntos
    gamma1_min(t*e-1*e+1:1*t*e,1)=res_deltau_min(:,1);
end
gamma1=[gamma1_max;gamma1_min];

for j=1:Np
    maty(j*q-q+1:j*q,1)=resy(1:q);
end
for j=1:Np
    maty(j*q-q+Np*q+1:j*q+Np*q,1)=resy(q+1:dim_y);
end
%% - ********************* Cálculo de la matriz H *********************- %%
H=(Phi'*Phi+R);
%% - ********************** MATRICES DE PARIDAD ***********************- %%
s = 2;
Gx = Construct_GX(Ad_cont,Cd_cont,s);
Gu = Construct_GU(Ad_cont,Bd_cont,Cd_cont,s,q,e);
Gf = Construct_GU(Ad_cont,Bf,Cd_cont,s,q,size(Bf,2));
[U,~,~] = svd(Gx);
U1 = U(:,1:2);
U2 = U(:,3:4);
Gxo = U2';
S = eye(2);
Rps = S*Gxo;
save('AJUSTAR_RUTA'  % Ruta original de Windows — ajustar al directorio de datos,'Gx',"Gu","Rps")
%% - ********************** Condiciones iniciales *********************- %%
xm=zeros(n1,1);                     %Planta a controlar
Xf=zeros(n1+q,1);                   %Modelo aumentado
u=zeros(e,1);                       %u(k-1)=0
y=zeros(q,1);
ures=zeros(e,1);
%% - ************************* SIMULATION LOOP ************************- %%
t = 0:Delta_t:(N_sim*Delta_t);
h1(1:2) = [0;0]; h2(1:2) = [0;0];
y(:,1:2) = zeros(2,2);
N = length(t);

for kkk=2:N-1
    Xf = [[h1(kkk);h2(kkk)] - ([h1(kkk-1);h2(kkk-1)]);[h1(kkk);h2(kkk)]];
    y(:,kkk) = [h1(kkk);h2(kkk)];
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
    u=calc_accion_control(e,Nc,H,Phi,BarRs,F,Xf,u,M,gamma);
    uplot(:,kkk) = u;
    %******************************************************************
    if kkk > s
        r(:,kkk) = Rps*(reshape(y(:,kkk-s:kkk) - [he1;he2],[],1) - Gu*reshape(uplot(:,kkk-s:kkk) - [ue1;ue2;ue3],[],1));
    else
        r(:,kkk) = [0;0];
    end
    f(:,kkk) = pinv(Rps*Gf)*r(:,kkk);
    %******************************************************************
    [h1(kkk+1), h2(kkk+1)] = twoTank_v2(u(1), u(2), u(3), t(kkk), t(kkk+1), h1(kkk), h2(kkk));
    h1(kkk+1) =  h1(kkk+1) + sigma*randn(1,1);
    h2(kkk+1) =  h2(kkk+1) + sigma*randn(1,1);

end
save('AJUSTAR_RUTA'  % Ruta original de Windows — ajustar al directorio de datos,'Gx',"Gu","Rps",'h1',"h2");

%% - ************************* RESPONSE GRAFICS ***********************- %%
k=0:(N_sim-1);
figure(1)
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
stairs(k,f','LineWidth',2);grid on
xlabel('Sampling Instant')
ylabel('[m^2]')
title('Leakage Section')
legend('\alpha1','\alpha2')
ax = gca;
ax.FontSize = 20;
set(ax, 'FontName', 'Times', 'FontSize',20,'FontWeight','normal','defaultAxesTickLabelInterpreter','latex')

figure(2)

stairs(k,r','LineWidth',2);grid on
xlabel('Sampling Instant')
ylabel('Coefficient')
title('Residual Generator')
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

if t>800 && t<=2000
    alpha1 = 0;
    alpha2 = 0e-6;
elseif t>2000 && t<=4000
    alpha1 = 0;
    alpha2 = 0e-5;
elseif t>4000 && t<6000
    alpha1 = 0;
    alpha2 = 0.5e-5;
elseif t>=6000
    alpha1 = 0e-5;
    alpha2 = 0e-5;
end

h1 = h(1);                                  % Tank 1 level
h2 = h(2);                                  % Tank 2 level
dh1 = 1/A1 * (u1*q0 - u2*a2*sqrt(2*g*h1)) - alpha1*sqrt(2*g*h1)/A1;% SE h1
dh2 = 1/A2 * (u2*a2*sqrt(2*g*h1) - u3*a3*sqrt(2*g*h2) - alpha2*sqrt(2*g*h2)); % SE for h2
dh = real([dh1; dh2]);
end

function Gx = Construct_GX(A,C,s)

[nC,nA] = size(C);
Gx = zeros(nC*(s+1),nA);

for i = 0 : s
    Gx(i*nC+1:(i+1)*nC,:) = C*A^i;
end

end

function Gu = Construct_GU(Am,Bm,Cm,s,q,e)
h = Cm;
for i=1:s-1
    h(q*i+1:q*i+q,:) = h(q*i+1-q:q*i,:)*Am;
end
v=h*Bm;
Gu = zeros(q*(s+1),e*(s+1));
Gu(q+1:end,1:e)=v;
for i=1:s-1
    Gu(:,e*i+1:e*i+e)=[zeros(i*q+2,e);v(1:q*(s-i),:)];
end

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

