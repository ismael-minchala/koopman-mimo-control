%% - **************** -TCLab Temperature Estimation- ******************* %%
%% - ******************** -Juan Francisco Duran S- ********************* %%
clear; close all; clc
rng(1234)
%% - ********************* -Declarar Variables- ************************ %%
n = 2; deltaT = 1; m = 2;
addpath(fullfile('..','utils'));
load(fullfile('..','datos','Data4_12h_filter.mat'));
set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultTextFontName', 'Times New Roman');
set(groot, 'defaultLegendFontName', 'Times New Roman');

set(groot, 'defaultAxesFontSize', 18);  % or any other desired size
set(groot, 'defaultTextFontSize', 18);
set(groot, 'defaultLegendFontSize', 18);

set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
%% - ******************** -Unforced Dynamics gen- ********************** %%

[U,~,~] = normalize_TS(U);
Umin = [0;0];
Umax = [100;100];
[Temp,~,~] = normalize_TS(Temp);
 
Tmax = [100;100];
Tmin = [0;0];
ind = 1;
n_traj = 40; 
j = 1;
for i = 1 :  (n_traj)
    Xunforced(ind:ind+n-1,:) = Temp(:,j*500:(j+1)*500);
    ind = ind + n; j = j + 2;
end
Xunforced = Xunforced(:,50:450);
t = 0:1:length(Xunforced)-1;
%% - **************** -Unforced Dynamics Arrange Data- ***************** %%
N_efun = 20;

trajLen = size(Xunforced,2);
Ntraj_free =  size(Xunforced,1)/n - 0;
Traj = cell(1,Ntraj_free);
index = 1;
for j = 1:Ntraj_free
    textwaitbar(j, Ntraj_free, "Collecting data without control\n")
    xx = Xunforced(index:index+n-1,1);
    for i = 2: trajLen
        xx = [xx Xunforced(index:index+n-1,i)];
    end
    index = index + n;
    Traj{j} = xx;
end
Traj_test = Traj{end};
Traj(end) = [];
for i = 1:n
    F_vec{i} = [];
    for j = 1:numel(Traj)
        F_vec{i} = [F_vec{i};Traj{j}(i,:).'];
    end
end
%% - ******************** -Unforced Dynamics DMD- ********************** %%
X_DMD = []; Y_DMD = [];
for i = 1: numel(Traj)
    X_DMD = [X_DMD, Traj{i}(:,1:end-1)];
    Y_DMD = [Y_DMD, Traj{i}(:,2:end)];
end
A_DMD = Y_DMD*pinv(X_DMD);
% DMD eigenvalues
lam_dt = eig(A_DMD);

% Generate lattice of from DMD eigenvalues
deg = 20;
pows = monpowers(numel(lam_dt),deg);
lam_dt_lattice = [];
for i = 1:size(pows,1)
    lam_dt_lattice = [ lam_dt_lattice ;  prod(lam_dt.'.^pows(i,:)) ];
end
%lam_dt_lattice = delete_non_conjugate_pairs(lam_dt_lattice);
LAM_OPT{1} = lam_dt_lattice(1:N_efun/2);
LAM_OPT{2} = lam_dt_lattice(1:N_efun/2);

%% - ************************ -Forced Response- ************************ %%
n = 2;
forced_response = Temp(:,1:43000);
Uf  = smoothdata(U(:,1:43000),2,"gaussian");
Xforced = reshape(forced_response,n,[],40);
Xtest = Xforced(:,1:350,end);
Xforced = Xforced(:,1:350,1:39);
Ubig = reshape(Uf,m,[],40);
Utest = Ubig(:,1:350,end);
Ubig = Ubig(:,1:350,1:39);
ttest = 0:deltaT:length(Utest)-1;

Nsim = size(Xforced,2) - 1;
Ntraj_control =  size(Xforced,3);
X_Data = zeros(n,Nsim+1,Ntraj_control);
U_Data = zeros(m,Nsim,Ntraj_control);
X_Data(:,1,:) = Xforced(:,1,:);

for i = 1 : Nsim
    U_Data(:,i,:) = Ubig(:,i,:);
    X_Data(:,i+1,:) = Xforced(:,i+1,:);
end

% Convert to cell array
Xtraj = cell(1,Ntraj_free); Utraj = cell(1,Ntraj_control);
for i = 1:Ntraj_control
    Xtraj{i} = X_Data(:,:,i);
    Utraj{i} = U_Data(:,:,i);
end
%% - ******************** -Hyperparameters selection- ****************** %%
dim_lambda1 = 1;
dim_lambda2 = 1;
dim_lambda3 = 1;

vars = [
    optimizableVariable('lambda1',[0 1])
    optimizableVariable('lambda2',[0 1])
    optimizableVariable('lambda3',[0 1])
    ];

cost_fn = @(tbl) constrained_nonrecurrent_sets( ...
    F_vec, Ntraj_free, Ntraj_control, trajLen, ...
    Xtraj, Utraj, Xtest, Utest, ...
    [tbl.lambda1, tbl.lambda2, tbl.lambda3]);

results = bayesopt(cost_fn, vars,'AcquisitionFunctionName','expected-improvement-plus', ...
                       'MaxObjectiveEvaluations',50, ...
                       'Verbose',1, ...
                       'UseParallel',true);

% Mejor solución encontrada
bestParams = bestPoint(results);
xbest = [bestParams.lambda1, bestParams.lambda2, bestParams.lambda3, 2*bestParams.z];
fbest = results.MinObjective;

%% - ****************************** -Results- ************************** %%

% === Datos de entrada ===
p_all = eig(A);                      % Todos los polos
p_all = p_all(:);                    % Vector columna

p_ref = lam_dt(:).';                % Polos de referencia (vector fila)

% === Cálculo de distancia e intensidad ===
distancias = min(abs(p_all - p_ref), [], 2);
intensidad = 1 ./ (1 + distancias);  % Mayor intensidad = más cercano

% === Detección automática de región densa ===
D = abs(p_all - p_all.');
D(D == 0) = Inf;                     % ignorar distancia a sí mismo

epsilon = 0.05;                      % umbral de vecindad
densidad = sum(D < epsilon, 2);     % número de vecinos cercanos

[~, idx_max] = max(densidad);       % índice del polo más densamente rodeado
centro_zoom = p_all(idx_max);       % centro del zoom

delta = 0.05;                        % tamaño de la ventana de zoom
a = real(centro_zoom) - delta;
b = real(centro_zoom) + delta;
c = imag(centro_zoom) - delta;
d = imag(centro_zoom) + delta;

% === Gráfica principal ===
figure(6);
h_polos = scatter(real(p_all), imag(p_all), 50, intensidad, 'filled');
colormap('turbo'); colorbar;
xlabel('Real part'); ylabel('Imaginary Part');
title('Pole and Zero Diagram');
axis equal; grid on; hold on;

% Circunferencia unitaria
theta = linspace(0, 2*pi, 500);
plot(cos(theta), sin(theta), 'k--', 'LineWidth', 1.5, 'HandleVisibility','off');

% Polos de referencia
h_ref = plot(real(p_ref), imag(p_ref), 'kx', 'LineWidth', 2, 'MarkerSize', 10);
% === Cuadro de zoom en gráfica principal ===
rectangle('Position', [a, c, 2*delta, 2*delta], ...
          'EdgeColor', 'g', 'LineStyle', '--', 'LineWidth', 2);

% === Inset (axes para zoom) ===
axesZoom = axes('Position', [0.43, 0.22, 0.35, 0.17]);
scatter(axesZoom, real(p_all), imag(p_all), 50, intensidad, 'filled');
colormap(axesZoom, 'turbo');
axis(axesZoom, 'manual');
set(axesZoom, 'XLim', [a b], 'YLim', [c d]);
box(axesZoom, 'on');

% Polos de referencia también en el zoom
hold(axesZoom, 'on');
plot(axesZoom, real(p_ref), imag(p_ref), 'yx', 'LineWidth', 2, 'MarkerSize', 10);
grid on
legend([h_polos, h_ref], {'Koopman Poles', 'Linear Poles'}, 'Location', 'best');

function J = constrained_nonrecurrent_sets(F_vec,Ntraj_free,Ntraj_control,trajLen,Xtraj,Utraj,Xtest,Utest,parameters)
% Initial lambda
lam0 = 0.8 - (0.8-.79)*rand(20*2,1); 

% Optimization function
optimize_lam_fun_1 = @(x)(eigOptim_grad(x,F_vec{1},Ntraj_free-1,trajLen,parameters(1)));
optimize_lam_fun_2 = @(x)(eigOptim_grad(x,F_vec{2},Ntraj_free-1,trajLen,parameters(2)));

% Options

options = optimoptions('fmincon','Display','off','MaxFunEvals',1e6, 'SpecifyObjectiveGradient',true);

lam_opt_1 = fmincon(optimize_lam_fun_1, lam0, [], [], [], [], ...
    -1*ones(numel(lam0),1), 1*ones(numel(lam0),1), ...
    @unit_circle_constraint, options);
lam_opt_2 = fmincon(optimize_lam_fun_2, lam0, [], [], [], [], ...
    -1*ones(numel(lam0),1), 1*ones(numel(lam0),1), ...
    @unit_circle_constraint, options);
% % Get back minimizers

LAM_OPT{1} = lam_opt_1(1:20/2) + 1i*lam_opt_1((20/2)+1:end);
LAM_OPT{2} = lam_opt_2(1:20/2) + 1i*lam_opt_2((20/2)+1:end);
%% - *************** -Unforced Dynamics lifting function *************** %%

lam_powers_traj = cell(n,numel(LAM_OPT{1}));
for i = 1:n
    for j = 1:numel(LAM_OPT{i})
        lam_powers_traj{i,j} = bsxfun(@power,LAM_OPT{i}(j),0:trajLen-1);
    end
end

% Powers of eigenvalues along the trajectory
lam_powers_traj = cell(n,numel(LAM_OPT{1}));
for i = 1:n
    for j = 1:numel(LAM_OPT{i})
        lam_powers_traj{i,j} = bsxfun(@power,LAM_OPT{i}(j),0:trajLen-1);
    end
end


%Build linear operators mapping g0 to phi
Lbig = cell(1,n);
Lbig = [];
L = cell(n,numel(LAM_OPT{i}));
for i = 1:n
    Lbig{i} = [];
    for j = 1:numel(LAM_OPT{i})
        L{i,j} = bdiag(lam_powers_traj{i,j}.',Ntraj-1);
        Lbig{i} = [Lbig{i}, L{i,j}];
    end
end

for i = 1:n
    g00{i} = (Lbig{i}'*Lbig{i}) \ (Lbig{i}'*F_vec{i}); % Seems to be numerically more stable
    for j = 1:numel(LAM_OPT{i})
        g0{i,j} = g00{i}((j-1)*(Ntraj_free-1)+1:j*(Ntraj-1)); % Initial values for each eigenfunction separately
    end
end

X = [];
for j = 1:numel(Traj)
    X = [X Traj{j}];
end

%% Values of eigenfunctions

Val_phi = cell(n,numel(LAM_OPT{1}));
for i = 1:n
    for j = 1:numel(LAM_OPT{1})
        Val_phi{i,j} = (L{i,j}*g0{i,j}).';
    end
end


%% Interpolate

VAL_PHI_CONCAT = [];
for i = 1:size(Val_phi,1)
    for k = 1:size(Val_phi,2)
        phi{i,k} = scatteredInterpolant(X',Val_phi{i,k}');
        phi{i,k} = @(x)(phi{i,k}(x')');
        VAL_PHI_CONCAT = [VAL_PHI_CONCAT ; Val_phi{i,k} ];
    end
end

% Vectorize
ind = 1;
phi_vec = cell(1,N_efun);
for i = 1:n
    for j = 1:numel(LAM_OPT{i})
        phi_vec{ind} = phi{i,j};
        ind = ind + 1;
    end
end

%% Lifting function
liftFun = @(x)( liftFun_function(x,phi_vec) );
%% Matrices A and C
C = bdiag(ones(1,size(lam0,1)/2),n);  % Could be obtained through regression as well
A = diag(reshape(cell2mat(LAM_OPT),[],1)); % Put all eigenvalues on the diagonal

Obs =  obsvk(A,C,Nsim+1); % Not an efficient implementation
b = [];
nc = size(C,1);
Q = zeros(Ntraj*Nsim*nc,size(A,1)*m);
for q = 1 : Ntraj_control
    x0 = Xtraj{q}(:,1);
    Obsx0 = (Obs*liftFun(x0));
    for j = 1:Nsim
        b = [b ; Obsx0( j*nc + 1 : (j+1)*nc, : ) - Xtraj{q}(:,j+1)] ;
        tmp = 0;
        for k = 0 : j-1
            kprime = j - k -1;
            tmp = tmp + kron(Utraj{q}(:,k+1)',Obs(kprime*nc + 1 : (kprime+1)*nc,:));
        end
        Q((q-1)*Nsim*nc + (j-1)*nc + 1 : (q-1)*Nsim*nc + j*nc,:) = tmp;
    end
end
b = -b;
alpha = parameters(3);            % elegir mediante CV o L-curve
Breg = (Q'*Q + alpha*eye(size(Q,2))) \ (Q'*b);
B = reshape(Breg,20,2);

Npred = size(Xtest,2);
Xlift = liftFun(Xtest(:,1));
for i = 1 :Npred-1
    Xlift(:,i+1) = A*Xlift(:,i) + B*Utest(:,i);
end
Xpred = real(C*Xlift);
Xtrue = Xtest;
J = 100*norm(Xtrue'-Xpred')/norm(Xtrue);

figure(4)
lw = 2;
subplot(211)
plot([0:Npred-1]*deltaT,Xtrue(1,:),'b','linewidth',lw); hold on
plot([0:Npred-1]*deltaT,Xpred(1,:),'r','linewidth',lw);
legend('True','Predicted');
ylim([min(min(Xtrue(1,:)), min(Xpred(1,:)))-0.1,max(max(Xtrue(1,:)), max(Xpred(1,:)))+0.1])
xlabel('Sample');
ylabel('$[^\circ C]$');
title('Normalized Temperature Prediction in Heatsink 1')
subplot(212)
plot([0:Npred-1]*deltaT,Xtrue(2,:),'b','linewidth',lw); hold on
plot([0:Npred-1]*deltaT,Xpred(2,:),'r','linewidth',lw);
xlabel('Sample');
ylabel('[$^\circ$ C]');
title('Normalized Temperature Prediction in Heatsink 2')
grid on
end

function [c, ceq] = unit_circle_constraint(lams)
n = numel(lams)/2;  % Número de valores propios complejos

ceq = [];
for k = 1:2:n-1
    r1 = lams(k);
    i1 = lams(k+n);
    r2 = lams(k+1);
    i2 = lams((k+1)+n);
    ceq(end+1) = r1 - r2;
    ceq(end+1) = i1 + i2;
end

lams = lams(1:n) + 1i*lams(n+1:end);
c = abs(lams') - 1;
end


function [Y_norm,Xmin,Xmax] = normalize_TS(X)
Xmin = 0; %min(X, [], 2);  % mínimo por fila
Xmax = 100; %max(X, [], 2);  % máximo por fila

Y_norm = (X - Xmin) ./ (Xmax - Xmin);  % normalización manual

end

function X = desnormalizeRange(Xnorm, Xmin, Xmax)
    X = Xnorm .* (Xmax - Xmin) + Xmin;
end

