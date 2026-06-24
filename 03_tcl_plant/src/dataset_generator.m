
clear;
% addpath('AJUSTAR_RUTA'  % Ruta original de Windows — ajustar al directorio de datos)
load('AJUSTAR_RUTA'  % Ruta original de Windows — ajustar al directorio de datos);
addpath('AJUSTAR_RUTA'  % Ruta original de Windows — ajustar al directorio de datos);

tclab;

% Parámetros
Tsim = 36000;
deltaT = 1;
Nsim = Tsim / deltaT;
ventana = 20;              % pasos para ventana de análisis
tiempo_min_cada_modo = 600; % [s] tiempo mínimo en cada modo

% Inicialización
h1(0); h2(0);
Temp(:,1) = [T1C(); T2C()];
U(:,1) = [round(100*rand()); round(100*rand())];  % entradas iniciales
modo = "calentando";
t_modo_actual = 1;
liftFun = @(x)(liftFun_function(x,phi_vec));
x_pred(:,1) = liftFun(normalize_TS(Temp(:,1),Tmax,Tmin));

for i = 2:Nsim
    % Medición
    Temp(:,i) = [T1C(); T2C()];
    difT(:,i) = Temp(:,i) - Temp(:,i-1);

    pasos_en_modo = i - t_modo_actual;

    switch modo
        case "calentando"
            h1(U(1,i-1)); h2(U(2,i-1));
            if i > ventana && (pasos_en_modo >= tiempo_min_cada_modo || ...
               (all(abs(mean(difT(:, i-ventana+1:i), 2)) < 0.01) && ...
                all(std(difT(:, i-ventana+1:i), 0, 2) < 0.01)))
                U(:,i) = [0; 0];
                h1(0); h2(0);
                modo = "enfriando";
                t_modo_actual = i;
                disp("-> Cambio a modo ENFRIANDO");
            else
                U(:,i) = U(:,i-1);
            end

        case "enfriando"
            h1(0); h2(0);
            U(:,i) = [0; 0];
            if i > ventana && (pasos_en_modo >= tiempo_min_cada_modo || ...
               (all(abs(mean(difT(:, i-ventana+1:i), 2)) < 0.01) && ...
                all(std(difT(:, i-ventana+1:i), 0, 2) < 0.01)))
                u1 = round(100*rand());
                u2 = round(100*rand());
                U(:,i) = [u1; u2];
                modo = "calentando";
                t_modo_actual = i;
                disp("-> Cambio a modo CALENTANDO");
            end
    end
    U_kop(:,i) = normalize_TS(U(:,i),Umax,Umin);
    x_pred(:,i) = A*x_pred(:,i-1) + B*U_kop(:,i);
    ypred(:,i) = desnormalizeRange(real(C*x_pred(:,i-1)),Tmin,Tmax);
    % Graficar cada 20 s
    if mod(i,20) == 0
        figure(1); clf;
        subplot(321); plot(1:i,Temp(1,:),1:i,ypred(1,:)); ylabel('T1 (°C)'); title('Temperatura 1');
        subplot(322); plot(1:i,Temp(2,:),1:i,ypred(2,:)); ylabel('T2 (°C)'); title('Temperatura 2');
        subplot(323); plot(U(1,:)); ylabel('u1 (%)'); title('Calentador 1');
        subplot(324); plot(U(2,:)); ylabel('u2 (%)'); title('Calentador 2');
        subplot(325); plot(difT(1,:)); ylabel('dT1'); title('Derivada T1');
        subplot(326); plot(difT(2,:)); ylabel('dT2'); title('Derivada T2');
        drawnow;
    elseif mod(i,3600) == 0
        % save("Data1_10h.mat","U","Temp","difT")
    end

    pause(deltaT);
end

h1(0); h2(0);
disp("Simulación completada.");

function [Y_norm] = normalize_TS(X,Xmax,Xmin)

Y_norm = (X - Xmin) ./ (Xmax - Xmin);  % normalización manual

end

function X = desnormalizeRange(Xnorm, Xmin, Xmax)
    X = Xnorm .* (Xmax - Xmin) + Xmin;
end