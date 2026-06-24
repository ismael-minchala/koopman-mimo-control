clear; clc; close all
%% -**************************** Load Data *******************************-
load('AJUSTAR_RUTA'  % Ruta original de Windows — ajustar al directorio de datos);
% load('AJUSTAR_RUTA'  % Ruta original de Windows — ajustar al directorio de datos);
%% -********************* Variables of Interest **************************-
n = size(A,1);
m = size(B,2);
p = size(C,1);
deltaT = 0.001;
mr = 0.095; mp = 0.024; Lr = 0.085; Lp = 0.129; Jr = mr*Lr^2/3; br = 2.0282e-4;
bp = 4.7743e-14; Jp = mp*Lp^2/3; g = 9.81;
Rm = 8.4; km = 0.042; kt = 0.042; Jp_cm = mp*Lp^2/12;
l = Lp/2; r = 0.085;
Jt = Jr*Jp - mp^2*l^2*r^2;
Bt = 7.7016;
%% -*********************** Balanced System ******************************-
sys = ss(A,B,C,0,deltaT);
[sysb,G] = balreal(sys);
sysr = sysb; %xelim(sysb,6);
% sysr = sys;
%% -***************** Controlabillity - Controller ***********************-
Ahat = sysr.A;
Bhat = sysr.B;
Chat = sysr.C;
Q = 100*eye(6);
R = 1;
K = dlqr(Ahat,Bhat,Q,R);
figure(2)
plot(eig(Ahat),'*');
hold on
plot(eig(Ahat - Bhat*K),'*','Marker','square');
legend(["Original poles" "Moved poles"])
%% -******************** Observabillity - Observer ***********************-
% Q = full(C'*C);  % Penalización en el error de estimación
Q = 0.001*eye(6);
R = 1*eye(size(C,1));  % Penalización en la ganancia
L = dlqr(Ahat', Chat', Q, R)'; % Ganancia óptima del observador

figure(3)
plot(eig(Ahat - L*Chat),'*');
hold on
plot(eig(Ahat - Bhat*K),'*','Marker','square');
legend(["Observer poles" "Controller poles"])
