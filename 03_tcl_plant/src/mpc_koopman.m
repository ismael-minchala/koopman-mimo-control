%% - ***************** SIMULACION MPC TANKS MODEL *********************- %%
%% - ********* AND FAULT ESTIMATION BASED ON KALMAN FILTER ************- %%
%% - **************** JUAN FRANCISCO DURAN SIGUENZA *******************- %%
clear; close all; clc
rng(1234);
Model = load('AJUSTAR_RUTA'  % Ruta original de Windows — ajustar al directorio de datos);
load('parameters_controller.mat');

set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultTextFontName', 'Times New Roman');
set(groot, 'defaultLegendFontName', 'Times New Roman');

set(groot, 'defaultAxesFontSize', 18);  % or any other desired size
set(groot, 'defaultTextFontSize', 18);
set(groot, 'defaultLegendFontSize', 18);

set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

%% - ******************** CONTROLLER PARAMETERS ***********************- %%
Delta_t = 2;                       %Tiempo de muestreo
Np = xbest(1); %26;                           %Horizonte de predicción
Nc = xbest(2);  %4;                            %Horizonte de control
lambda1 = xbest(3);                       %Factor de penalización DeltaU1
lambda2 = xbest(4);                       %Factor de penalización DeltaU2

N_sim = 33200*3;
%yrr = [0.4*ones(1,N_sim);0.4*ones(1,N_sim)];
%yrr = [0.2*sin(2*pi*[1:N_sim] / N_sim);0.2*sin(2*pi*[1:N_sim] / N_sim)]+0.5;
%yrr = [[0.15*ones(1,floor(N_sim/2)) 0.05*ones(1,floor(N_sim/2))];[0.2*ones(1,floor(N_sim/2)) 0.1*ones(1,floor(N_sim/2))]];
% yrr = [[0.8*ones(1,floor(N_sim/2)) 0.8*ones(1,floor(N_sim/2))];[0.8*ones(1,floor(N_sim/2)) 0.5*ones(1,floor(N_sim/2))]];
% t = 0:N_sim-1;
% 
% % Frecuencia para una onda por el total de la simulación
% f = 1/(N_sim/2); 
% 
% % Señal diente de sierra entre 0 y 1
% signal = (sawtooth(2*pi*f*t, 1) + 1) / 2;
% 
% % Escalar y desplazar para que sea similar a tu ejemplo (amplitud 0.2 y offset 0.5)
% yrr = [0.2*signal + 0.5; 0.2*signal + 0.5];

min_duracion = 400;
valor_min = 0;
valor_max = 100;

% Inicializamos la matriz de referencias
yrr = zeros(2, N_sim);

for fila = 1:2
    idx = 1;
    while idx <= N_sim
        % Generamos un valor aleatorio entre 0 y 85
        valor = round(valor_min + (valor_max - valor_min) * rand())/100;
        % Elegimos una duración aleatoria >= min_duracion
        duracion = randi([min_duracion, 2*min_duracion]);
        % Nos aseguramos de no exceder el límite de muestras
        duracion = min(duracion, N_sim - idx + 1);
        % Llenamos el bloque
        yrr(fila, idx:idx+duracion-1) = valor;
        idx = idx + duracion;
    end
end

sp1 = yrr(1,1)*ones(1,Np);
sp2 = yrr(2,1)*ones(1,Np);
%% - ************************ RESTRICCIONES ***************************- %%
deltaumin1 = -.1;                %DeltaU1 mínima
deltaumax1 = .1;                 %DeltaU1 máxima
deltaumin2 = -.1;                %DeltaU2 mínima
deltaumax2 = .1;                 %DeltaU2 máxima

umin1 = 0;                      %U1 mínimo
umax1 = 1;                      %U1 máximo
umin2 = 0;                       %U2 mínimo
umax2 = 1;                        %U2 mínimo

ymin1 = 0;                      %Y1 mínimo
ymax1 = 1;                       %Y1 máximo
ymin2 = 0;                      %Y2 mínimo
ymax2 = 1;                       %Y2 máximo

puntos = Nc;
%% %Modelo de predicción.
A = Model.A;
B = Model.B;
C = Model.C;
D = 0;

% Supongamos que A, B, C, D son complejos
A_r = real(A);  A_i = imag(A);
B_r = real(B);  B_i = imag(B);
C_r = real(C);  C_i = imag(C);
D_r = real(D);  % D_i no se usa si salida es real

% Construcción de sistema real equivalente
Ad_cont = [A_r, -A_i; A_i, A_r];
Bd_cont = [B_r; B_i];
Cd_cont = [C_r, -C_i];
Dd_cont = 0;
liftFun = @(x)(liftFun_function(x,Model.phi_vec));

%% Dimensiones del proceso
n1 = size(Ad_cont,1);           %número de estados
q = size(Cd_cont,1);            %número de salidas
e = size(Bd_cont,2);            %número de entradas
%% %Modelo aumentado (integrador embebido)
[Am,Bm,Cm,Dm] = mod_aumentado_mimo(Ad_cont,Bd_cont,Cd_cont,Dd_cont,n1,q);
%% Observer
C = full(Cd_cont);
Q = eye(40); diag(xbest(5:44));
R = diag(xbest(45:end));  % Penalización en la ganancia
L = dlqr(Ad_cont', Cd_cont', Q, R)'; % Ganancia óptima del observador
%% %Cálculo de las matrices F y Phi
[F,Phi] = mat_f_phi(Am,Bm,Cm,Np,Nc,q,e);
%% %Adecuaciones de R y Puntos de ajuste para el caso multivariable
BarRs = [sp1;sp2];
BarRs=reshape(BarRs,q*Np,1);            %Vector de consignas
Rbase=eye(e*Nc,e*Nc);
lambda_aux=[lambda1;lambda2];
lambda=repmat(lambda_aux,Nc,1);
R=lambda.*Rbase;
%% % Preparación de las matrices de restricciones
res_deltau_max=[deltaumax1;deltaumax2];
res_deltau_min=[-deltaumin1;-deltaumin2];
resu=[umax1;umax2;-umin1;-umin2];
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
Xf = zeros(42,1);
yreal(:,1) = 17*ones(q,1);
Xpred(:,1) = zeros(40,1); %[real(liftFun(normalize_TS(yreal(:,1),Model.Tmax,Model.Tmin)));imag(liftFun(normalize_TS(yreal(:,1),Model.Tmax,Model.Tmin)))];
yhat(:,1) = C*Xpred(:,1);
u = zeros(e,1);                       %u(k-1)=0
ures = zeros(e,1);
ypred_norm(:,1) = zeros(2,1);
%ALGORITMO ITERATIVO
%*******************************************************************
t = 0:Delta_t:(N_sim-1)*Delta_t;
N = length(t);
for kkk = 2:N
    tic
    sp1 = yrr(1,kkk)*ones(1,Np);
    sp2 = yrr(2,kkk)*ones(1,Np);
    BarRs = [sp1;sp2];
    BarRs=reshape(BarRs,q*Np,1);

    for i=1:puntos
        gamma2(i*e-e+1:i*e,:)= -u+resu(1:e);
        gamma2(i*e-e+e*puntos+1:i*e+e*puntos,:)=u+resu(e+1:2*e);
    end

    compy=[-F*(Xf);F*(Xf)];
    gamma3=maty+compy;
    %FORMULACION DE GAMMA TOTAL
    %********************************************************************
    gamma=[gamma1;gamma2;gamma3];
    %********************************************************************
    u = calc_accion_control(e,Nc,H,Phi,BarRs,F,Xf,u,M,gamma);
    u = max(min(u, 1), 0);
    uplot(:,kkk) = u;
    u_apply = desnormalizeRange(u,Model.Umin,Model.Umax);

    sol = ode45(@(t,x)heat(t,x,u_apply(1),u_apply(2)),[t(kkk-1) t(kkk)],yreal(:,kkk-1));
    yreal(:,kkk) = sol.y(:,end);
    yreal_norm(:,kkk) = normalize_TS(yreal(:,kkk),Model.Tmax,Model.Tmin);
    
    Xpred(:,kkk) = Ad_cont*Xpred(:,kkk-1) + Bd_cont*[u(1);u(2)] + L*(yreal_norm(:,kkk) - yhat(:,kkk-1));
    yhat(:,kkk) = Cd_cont*Xpred(:,kkk);

    Xf = [Xpred(:,kkk) - Xpred(:,kkk-1);yreal_norm(:,kkk)];   
    aux(kkk) = toc;
end
yhat = desnormalizeRange(yhat,Model.Tmin,Model.Tmax);
%% - ********************* -Response of Control- *********************** %%
t = 0:Delta_t:(N_sim-1)*Delta_t;
yrrd = desnormalizeRange(yrr,Model.Tmin,Model.Tmax);
figure(1)
ejes(1) = subplot(221);
plot(t,yreal(1,:),t,yhat(1,:),t,yrrd(1,:),'--','LineWidth',2);grid on
xlabel('Time [s]')
ylabel('Temperature [C]')
legend(["Heatsink 1" "Estimated" "Reference"])
ax = gca;
ax.FontSize = 20;
set(ax, 'FontName', 'Times', 'FontSize',20,'FontWeight','normal','defaultAxesTickLabelInterpreter','latex')

ejes(2) = subplot(222);
plot(t,yreal(2,:),t,yhat(2,:),t,yrrd(2,:),'--','LineWidth',2);grid on
xlabel('Time [s]')
ylabel('Temperature [C]')
legend(["Heatsink 2" "Estimated" "Reference"])
ax = gca;
ax.FontSize = 20;
set(ax, 'FontName', 'Times', 'FontSize',20,'FontWeight','normal','defaultAxesTickLabelInterpreter','latex')

ejes(3) = subplot(223);
stairs(t,uplot(1,:),'LineWidth',2);grid on
xlabel('Time [s]')
legend('Heater 1 ')
title('Control Singal')
ax = gca;
ax.FontSize = 20;
set(ax, 'FontName', 'Times', 'FontSize',20,'FontWeight','normal','defaultAxesTickLabelInterpreter','latex')

ejes(4) = subplot(224);
stairs(t,uplot(2,:),'LineWidth',2);grid on
xlabel('Time [s]')
legend('Heater 2')
title('Control Singal')
ax = gca;
ax.FontSize = 20;
set(ax, 'FontName', 'Times', 'FontSize',20,'FontWeight','normal','defaultAxesTickLabelInterpreter','latex')
linkaxes(ejes,'x')
%% - ******************* -Accuracy of the Control- ********************* %%
serie_ref = yrrd;
serie_out = yreal;
Ts = 2;
[m, N] = size(serie_ref);
assert(m == 2, 'Solo soporta referencias 2xN');

% 1. Detectar indices de cambio (cuando cambia cualquier estado)
dif = any(diff(serie_ref,1,2) ~= 0, 1);
idx_cambios = [1 find(dif) + 1 N];

% 2. Calcular ITAE por bloque
r1 = [];
r2 = [];
itae = [];

for k = 1:length(idx_cambios)-1
    idx_ini = idx_cambios(k);
    idx_fin = idx_cambios(k+1)-1;

    % referencia constante en este bloque
    r = serie_ref(:,idx_ini);

    % salida del sistema en este bloque
    y = serie_out(:,idx_ini:idx_fin);

    % error absoluto
    e = vecnorm(y - r,2,1);

    % tiempo relativo
    t = (0:length(e)-1) * Ts;

    % ITAE
    ITAE = sum(t .* abs(e)) * Ts;

    r1(end+1) = r(1);
    r2(end+1) = r(2);
    itae(end+1) = ITAE;
end

resultados.r1 = r1;
resultados.r2 = r2;
resultados.itae = normalize(itae,2,'range') * 100;

xq = linspace(min(resultados.r1), max(resultados.r1), 100);
yq = linspace(min(resultados.r2), max(resultados.r2), 100);
[Xq, Yq] = meshgrid(xq, yq);

% Interpolación (natural o linear)
Zq = griddata(resultados.r1, resultados.r2, resultados.itae, Xq, Yq, 'natural');

Zq_filled = Zq;
Zq_filled(isnan(Zq)) = max(resultados.itae);  % Rellenar con valor alto (peor error)

% Graficar con colores suaves
figure(3);
subplot(121)
contourf(Xq, Yq, Zq_filled, 30, 'LineColor','none');
colormap(jet);   % Puedes cambiar a 'jet', 'hot', etc.
% xlabel({'Temperature in Heatsink 1 [$^\circ$C]'},{'(a)'}, 'Interpreter', 'latex');
xlabel(sprintf('Temperature in Heatsink 1 [$^\\circ$C]\n(a)'), 'Interpreter', 'latex');
ylabel('Temperature in Heatsink 2 [$^\circ$C]', 'Interpreter', 'latex');
title('Normalized ITAE');
c = colorbar;
c.Label.String = 'ITAE [%]';   % Texto de la leyenda
c.Label.Interpreter = "latex";
c.Label.FontSize = 18;
c.Label.FontName = 'Times New Roman';
%% - ********************* -Steady State Error- ************************ %%
serie_ref = yrrd;
serie_out = yreal;
[m, N] = size(serie_ref);
assert(m == 2, 'Solo soporta referencias 2xN');

Nwindow = 50;  % número de muestras para estado estable

r1_vals = [];
r2_vals = [];
ese_vals = [];

% --- Detectar todos los cambios de referencia ---
cambios1 = find([1 diff(serie_ref(1,:)) ~= 0]);
cambios2 = find([1 diff(serie_ref(2,:)) ~= 0]);
todos_cambios = unique([cambios1, cambios2, N]);  % incluir último índice
todos_cambios = sort(todos_cambios);

% --- Calcular ESE por bloque ---
for k = 1:length(todos_cambios)-1
    idx_ini = todos_cambios(k);
    idx_fin = todos_cambios(k+1)-1;

    % Promediar las últimas Nwindow muestras del bloque
    idx_ss = max(idx_ini, idx_fin-Nwindow+1):idx_fin;
    y_final = mean(serie_out(:, idx_ss), 2);

    % Referencias del bloque (tomamos las iniciales)
    r_block = serie_ref(:, idx_ini);

    % Error en estado estable (norma 2 normalizada)
    e_norm = norm((r_block - y_final)./r_block);

    % Guardar valores
    r1_vals(end+1) = r_block(1);
    r2_vals(end+1) = r_block(2);
    ese_vals(end+1) = e_norm;
end
eliminated_idx = find(ese_vals>1);

ese_vals(eliminated_idx) = [];
r1_vals(eliminated_idx) = [];
r2_vals(eliminated_idx) = [];
% --- Normalizar errores para colores ---
ese_norm = ese_vals*100;

% --- Crear malla 2D ---
xq = linspace(min(r1_vals), max(r1_vals), 100);
yq = linspace(min(r2_vals), max(r2_vals), 100);
[Xq, Yq] = meshgrid(xq, yq);

% --- Interpolación ---
Zq = griddata(r1_vals, r2_vals, ese_norm, Xq, Yq, 'natural');
Zq_filled = Zq;
Zq_filled(isnan(Zq)) = max(ese_norm);  % rellenar huecos

% --- Graficar heatmap ---
subplot(122)
contourf(Xq, Yq, Zq_filled, 30, 'LineColor','none');
colormap(jet);
colorbar;
c = colorbar;
xlabel(sprintf('Temperature in Heatsink 1 [$^\\circ$C]\n(b)'), 'Interpreter', 'latex');
ylabel('Temperature in Heatsink 2 [$^\circ$C]', 'Interpreter', 'latex');
title('Normalized Steady State Error');
c.Label.String = 'SSE [%]';   % Texto de la leyenda
c.Label.Interpreter = "latex";
c.Label.FontSize = 18;
c.Label.FontName = 'Times New Roman';

function X = desnormalizeRange(Xnorm, Xmin, Xmax)
    X = Xnorm .* (Xmax - Xmin) + Xmin;
end

function [Y_norm] = normalize_TS(X,Xmax,Xmin)

Y_norm = (X - Xmin) ./ (Xmax - Xmin);  % normalización manual

end

function results = ss_error_by_steps(y, r, tol_rel, tol_abs_min, consec)
% y, r: vectores (Nx1)
% tol_rel: tolerancia relativa (ej. 0.02 para 2%)
% tol_abs_min: tolerancia absoluta mínima (ej. 0.01 unidades)
% consec: número de muestras consecutivas dentro de la banda para considerar asentado

if nargin < 3, tol_rel = 0.02; end
if nargin < 4, tol_abs_min = 1e-3; end
if nargin < 5, consec = 50; end

N = numel(y);
% detectar bordes de cambio de referencia
edges = [1; find(diff(r)~=0)+1; N+1];

results = struct('r_value', {}, 'start', {}, 'stop', {}, 'settle_index', {}, 'e_ss', {}, 'e_ss_abs', {}, 'e_ss_pct', {});
cnt = 0;
for k = 1:(numel(edges)-1)
    start_idx = edges(k);
    stop_idx = edges(k+1)-1;
    r_val = r(start_idx);
    seg_y = y(start_idx:stop_idx);
    seg_len = numel(seg_y);

    % banda de tolerancia (abs)
    tol = max(tol_abs_min, tol_rel * abs(r_val));

    % buscar primer índice dentro de la banda que se mantenga 'consec' muestras
    settled = false;
    for j = 1:(seg_len - consec + 1)
        window = seg_y(j:(j+consec-1));
        if all(abs(window - r_val) <= tol)
            settle_idx = start_idx + j - 1; % índice relativo a y total
            settled = true;
            break;
        end
    end
    if ~settled
        % si no hay asentamiento detectado, usar la última fracción del segmento
        frac = 0.2; % usar último 20%
        settle_idx = start_idx + max(0, round((1-frac)*seg_len));
    end

    % calcular e_ss usando la última parte del segmento (desde settle_idx hasta stop_idx)
    idx_ss = settle_idx:stop_idx;
    e_ss = mean(r(idx_ss) - y(idx_ss));
    e_ss_abs = mean(abs(r(idx_ss) - y(idx_ss)));
    if r_val ~= 0
        e_ss_pct = 100 * e_ss_abs / r_val;
    else
        e_ss_pct = NaN;
    end

    cnt = cnt + 1;
    results(cnt).r_value = r_val;
    results(cnt).start = start_idx;
    results(cnt).stop = stop_idx;
    results(cnt).settle_index = settle_idx;
    results(cnt).e_ss = e_ss;
    results(cnt).e_ss_abs = e_ss_abs;
    results(cnt).e_ss_pct = e_ss_pct;
end
end
