%% - ********************** -PENDULO DE FURUTA- ************************ %%
%% - ******************** -Juan Francisco Duran S- ********************* %%
clear; close all; clc 
%% - ********************* -Declarar Variables- ************************ %%
mr = 0.095; mp = 0.024; Lr = 0.085; Lp = 0.129; Jr = mr*Lr^2/3; br = 2.0282e-4;
bp = 4.7743e-14; Jp = mp*Lp^2/3; g = 9.81;
Rm = 8.4; km = 0.042; kt = 0.042;
l = Lp/2;
r = 0.085;
Bt = 7.7016;
%% - **************** -ESTIMATION OF UNKNOW PARAMETERS- **************** %%
pendulum = load("pendulum_time_series.mat");
t = pendulum.ans(1,:);
theta = pendulum.ans(2,:);
alpha = pendulum.ans(3,:);
voltage = pendulum.ans(4,:);

%% - ********************* -MODELO LINEALIZADO- ************************ %%
sigma1 = J1*m1*l0^2 + m0*m1*r0^2*r1^2 + J1*m0*r0^2 + J0*m1*r1^2 + J0*J1;
sigma2 = m1*l0^2 + m0*r0^2 + J0;
sigma3 = D0*R+Ke*Km;
A = [0 0 1 0;
     0 0 0 1; 
     0 -(g*l0*m1^2*r1^2)/(sigma1) -((m1*r1^2 + J1)*sigma3)/(R*sigma1) D1*l0*m1*r1/sigma1;
     0 (g*m1*r1*sigma2)/sigma1 (l0*m1*r1*sigma3)/(R*sigma1) -D1*sigma2/sigma1];
sigma1 = R*(m1*(J1*l0^2+J0*r1^2+m0*r0^2*r1^2)) + J1*(m0*r0^2+J0);
B = [0;0;(Km*(m1*r1^2+J1))/sigma1; -(Km*l0*m1*r1)/sigma1];
C = [1 0 0 0;0 1 0 0];
figure 
sys = ss(A,B,C,0);
margin(sys(2,1))
%% - ************************* -Control LAW- *************************** %%
Ts = 1e-5;
Mp = 0.1;
%% - ************** -DETERMINACION DE POLOS DESEADOS- ****************** %%
zeta = sqrt(log(Mp)^2/(log(Mp)^2+pi^2));
sigma = 4/Ts;                                           %Criterio 2%
wn = sigma/zeta;
poli_carac = [1 2*zeta*wn +wn^2];
pol_des = [roots(poli_carac)' -10 -10 -10 -10];
%% - ************** -DETERMINACION DE CONTROLABILIDAD- ***************** %%
%% - ********************** -Y OBSERVABILIDAD- ************************* %%
Controlab = ctrb(A,B);
if rank(Controlab) == size(A,1)
    disp("El sistema es controlable")
else
    disp("El sistema NO es controlable")
end

Observab = obsv(A,C);
if rank(Observab) == size(A,1)
    disp("El sistema es observable")
else
    disp("El sistema NO es observable")
end
%% - *************** -AÑADIR EL ESTADO DE INTEGRACION- ***************** %%
Ai = [A zeros(2,4)';C*A eye(2,2)];
Bi = [B;C*B];
Ci = [zeros(2,4) eye(2,2)];
Di  =zeros(2,1);
K = acker(Ai,Bi,pol_des);
%% - ******************* -OBSERVADOR DE ESTADOS- *********************** %%
pol_observ = pol_des(1:size(A,1))*10;
KE = acker(A',C(1,:)',pol_observ);

