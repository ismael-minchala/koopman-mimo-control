%% - **************** -SINDy Frequency Deviation Prediction- *********** %%
%% - ************************ -Wind Turbine- ***************//********** %%
%% - ******************** -Juan Francisco Duran S- ********************* %%
clear; close all; clc
%% - ********************** -Fixed Parameters- ************************* %%
set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultTextFontName', 'Times New Roman');
set(groot, 'defaultLegendFontName', 'Times New Roman');

set(groot, 'defaultAxesFontSize', 18);  
set(groot, 'defaultTextFontSize', 18);
set(groot, 'defaultLegendFontSize', 18);

set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

%% - ************************* -Load Data- ***************************** %%
load('Datos_predictor.mat');
deltaT = Datos_predictor.time(2) - Datos_predictor.time(1); 

Data = reshape([Datos_predictor{1:150e3,2:end};Datos_predictor{240e3:320e3,2:end}],1,[]);
Data = 2 * (Data - min(Data)) / (max(Data) - min(Data)) - 1;


numSamples = length(Data)-1;

fprintf('Data Collection Done \n')

train = floor(numSamples*0.7);
Xnext = Data(1:train - 1000000)';
Xtest = Data(train+1:end)';
k = 0:1:length(Xnext)-1;
k2 = 0:1:length(Xtest)-1;
figure(1)

subplot(211)
plot(k,Xnext,'LineWidth',2)
title('Training Set','FontName', 'Times', 'FontSize',15,'FontWeight','normal')
xlabel('Sample')
ylabel('Frequency Deviation [Hz]')
grid on

subplot(212)
plot(k2,Xtest,'LineWidth',2)
title('Test Set','FontName', 'Times', 'FontSize',15,'FontWeight','normal')
xlabel('Sample')
ylabel('Frequency Deviation [Hz]')
grid on

%% - ************************ -BASIS FUNCTION- ************************* %%
X = Xnext(1:end-1);
Y = Xnext(2:end);
npol = 10;
n = size(Y,2);
%% - ************************ -SINDY ALGORITHM- ************************** %%
fprintf('Starting Lifting \n')
Theta = build_theta(X,npol);
Nlift = size(Theta,2);
lambda = 5e-4;%1e-2;
Xi = Theta\Y;
for numSamples=1:40     
    smallinds = (abs(Xi)<lambda); 
    Xi(smallinds) = 0; 
    for ind = 1:n 
        biginds = ~smallinds(:,ind);
        Xi(biginds,ind) = Theta(:,biginds)\Y(:,ind);
    end
end

fprintf( 'Regression residual %f \n', norm(Y - Theta*Xi,'fro') / norm(Y,'fro') );

%% - ************************* -Validation- **************************** %%
clear x
Nsim = 100e3;
Xprobe = Xtest(1:Nsim);
x(1) = Xprobe(1);
xlift = build_theta(x(1),npol);
for i = 2: Nsim
    x(i) = xlift*Xi;
    if mod(i,1) == 0 && i > 1
        xlift = build_theta(Xprobe(i),npol);
    end
end
error = norm(Xprobe' - x,'fro') / norm(Xprobe,'fro');
%% - ********************* -SINDY algorithm FIT- *********************** %%

figure(3)

plot(0:Nsim-1,Xprobe,'-b','LineWidth',2); hold on
plot(0:Nsim-1,x,'--r','LineWidth',2);
legend(["True" "Predictions"])
grid on

figure(4)

plot(error,'-b','LineWidth',2); hold on
legend(["Prediction error"])
grid on

error(isnan(error)) = [];
fprintf( 'Mean Error %f \n',mean(error)) ;

function theta = build_theta(X,n)
    Xs = X(:); 

    N = length(Xs);
    theta = zeros(N, n+1);

    theta(:,1) = 1;         % T0(x) = 1
    if n >= 1
        theta(:,2) = Xs;    % T1(x) = x
    end
    for j = 2:n
        % Tj(x) = 2*x*T_{j-1}(x) - T_{j-2}(x)
        theta(:,j+1) = 2 * Xs .* theta(:,j) - theta(:,j-1);
    end
end