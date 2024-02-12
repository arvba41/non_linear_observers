clear; clc; 

load("battery_sim_data/ocv_data.mat");
load("battery_sim_data/cyc_data.mat");

%% parameters
% models parameters
qocv = data.qbatt;
vocv = data.vbatt;
params.Ro = 1.5e-3; % internal resistance (estimate)
params.qn = 24;
% params.q0 = params.qn * cyc_data.Qinit * 1e-2;
params.q0 = params.qn * 0.7;

% measurements
q0_true = params.qn * cyc_data.Qinit;
tsim = cyc_data.tprc;
ibatt = cyc_data.ibatt_prc;
vbatt = cyc_data.vbatt_prc;
qbatt = cyc_data.qbatt_prc / 100 * params.qn;

ocv = SplineInterpol(qocv,vocv);

%% Battery model

% battery model (pure simulation with reasonable parameters)
fvocv = @(q) interp1(qocv,vocv,q,'spline'); % battery SOC-ocv function
Ro = params.Ro;                 % battery internal resistnace 
x_qb = @(x,u) - u / 3600;         % state equation
y_vb = @(x,u) fvocv(x) - u * Ro; % measurement equation

%% plotting ocv 
figure(100); clf; 

subplot(121); 
plot(qocv / params.qn * 1e2,vocv); box off;
ylabel('$v_b$ [V]','Interpreter','latex');
xlabel('$q_b$ [\%]','Interpreter','latex');

subplot(122); 
plot(qocv / params.qn * 1e2,(fvocv(qocv+0.01) - fvocv(qocv-0.01))/ 0.02 ); box off;
ylabel('$dv_b/dq$ [V/Ah]','Interpreter','latex');
xlabel('$q_b$ [\%]','Interpreter','latex');

hh = findall(gcf,'Type','Axes');
set(hh,'TickLabelInterpreter','latex');

% % %% EKF method
% % Ts = mean(diff(tsim)); % time step 
% % 
% % Guess = [10 0.1 10];
% % 
% % % model
% % modelEKF.f = @(x,u) x + Ts * x_qb(x,u);
% % % measurement 
% % modelEKF.h = @(x,u) y_vb(x,u);
% % % linearization f
% % modelEKF.fx = @(x) 1;
% % % linearization h
% % % modelEKF.hx = @(x,x_p) (fvocv(x) - fvocv(x_p)) ./ (x - x_p); % numeric differentiation on focv
% % modelEKF.hx = @(x) (fvocv(x + 0.01) - fvocv(x - 0.01)) ./ 0.02; % numeric differentiation on focv
% % 
% % % turing parameters          
% % modelEKF.Q = Guess(1); % model covriance matrix 
% % modelEKF.R = Guess(2); % measurement covriance matrix
% % init.P0 = Guess(3);  % Initial covariance
% % 
% % init.x0 = params.q0; % Initial battery capacity
% % 
% % % measuremetn data
% % mesh.y = vbatt;
% % mesh.u = ibatt;
% % 
% % q_EKF = EKF_batt(modelEKF,init,mesh);
% % 
% % % updating measuremet with the new estimate
% % vb_EKF = y_vb(q_EKF,ibatt);
% % 
% % % plots -- direct simulation results
% % figure(11); clf; set(gcf,'WindowStyle','docked');
% % 
% % tiledlayout(1,2)
% % 
% % ax(1) = nexttile();
% % plot(tsim / 3600,[qbatt; q_EKF] / params.qn * 1e2); box off;
% % ylabel('$q_b$ [\%]','Interpreter','latex');
% % xlabel('$t$ [h]','Interpreter','latex');
% % legend('tru','est','EKF');
% % legend('Interpreter','latex');
% % 
% % ax(2) = nexttile();
% % plot(tsim / 3600,[vbatt; vb_EKF]); box off;
% % ylabel('$v_b$ [V]','Interpreter','latex');
% % xlabel('$t$ [h]','Interpreter','latex');
% % 
% % linkaxes(ax,'x')
% % 
% % h = findall(gcf,'Type','Axes');
% % set(h,'TickLabelInterpreter','latex');
% % 
% % %% performing optimizaton to get the Q, R and P0 values
% % % measurment signals 
% % mesh.qbatt = qbatt;
% % 
% % % initial exit flag
% % exitflag = 0;
% % 
% % % initial guess
% % X0 = Guess; 
% % 
% % % otimization optionss
% % options = optimoptions('lsqnonlin','Display',...
% %     'iter','FunctionTolerance',1e-3);
% % 
% % % normal EKF
% % tStart = tic;
% % while ~exitflag
% %     [x,resnorm,residual,exitflag,output] = ...
% %         lsqnonlin(@(x) EKF_batt_opt(x,mesh,init,modelEKF), ...
% %         X0, zeros(size(X0)) + eps, [], options);
% %     X0 = x; % homotopy :)
% % end
% % tELF(2) = toc(tStart);
% % 
% % modelEKF.Q = x(1);
% % modelEKF.R = x(2);
% % init.P0 = x(3);
% % 
% % xhat_opt_EKF = EKF_batt(modelEKF,init,mesh);
% % 
% % %% plots -- optimal results
% % figure(12); clf; set(gcf,'WindowStyle','docked');
% % 
% % plot(tsim / 3600,[qbatt; xhat_opt_EKF] / params.qn * 1e2); box off;
% % ylabel('$q_b$ [\%]','Interpreter','latex');
% % xlabel('$t$ [h]','Interpreter','latex');
% % legend('tru','est','EKF');
% % legend('Interpreter','latex');
% % 
% % h = findall(gcf,'Type','Axes');
% % set(h,'TickLabelInterpreter','latex');
% % 
% % 
%% Coulumb counting

% simulated battery voltage iwith noisy data as input
N = length(tsim);
h = mean(diff(tsim)); % sample time
% initial value
qb_sim = zeros(size(tsim));
vb_sim = zeros(size(tsim));
qb_sim(1) = params.q0; 
vb_sim(1) = y_vb(qb_sim(1),ibatt(1));
% euler forward simulation 
for ii = 2:N
    qb_sim(ii)  = qb_sim(ii-1) + h * x_qb(qb_sim(ii-1),ibatt(ii)); 
    vb_sim(ii)  = y_vb(qb_sim(ii),ibatt(ii));
end

%% plots -- direct simulation results
figure(10); clf; set(gcf,'WindowStyle','docked');

tiledlayout(1,2)

ax(1) = nexttile();
plot(tsim / 3600,[qbatt; qb_sim]); box off;
ylabel('$q_b$ [Ah]','Interpreter','latex');
xlabel('$t$ [h]','Interpreter','latex');
legend('tru','est');
legend('Interpreter','latex');

ax(2) = nexttile();
plot(tsim / 3600,[vbatt; vb_sim]); box off;
ylabel('$v_b$ [V]','Interpreter','latex');
xlabel('$t$ [h]','Interpreter','latex');

linkaxes(ax,'x')

hh = findall(gcf,'Type','Axes');
set(hh,'TickLabelInterpreter','latex');

%% EKF method
Ts = mean(diff(tsim)); % time step 

Guess = [0.1 1 1];

% model
modelEKF.f = @(x,u) x - Ts * u / 3600;
% measurement 
modelEKF.h = @(x,u) ocv.f(x) - u * Ro;
% linearization f
modelEKF.fx = @(x) 1;
% linearization h
% modelEKF.hx = @(x,x_p) (fvocv(x) - fvocv(x_p)) ./ (x - x_p); % numeric differentiation on focv
modelEKF.hx = @(x) ocv.df(x); % numeric differentiation on focv

% turing parameters          
modelEKF.Q = Guess(1); % model covriance matrix 
modelEKF.R = Guess(2); % measurement covriance matrix
init.P0 = Guess(3);  % Initial covariance

init.x0 = params.q0; % Initial battery capacity

% measuremetn data
mesh.y = vbatt;
mesh.u = ibatt;

q_EKF = EKF_batt(modelEKF,init,mesh);

% updating measuremet with the new estimate
vb_EKF = y_vb(q_EKF,ibatt);

% plots -- direct simulation results
figure(11); clf; set(gcf,'WindowStyle','docked');

tiledlayout(1,2)

ax(1) = nexttile();
plot(tsim / 3600,[qbatt; qb_sim; q_EKF] / params.qn * 1e2); box off;
ylabel('$q_b$ [\%]','Interpreter','latex');
xlabel('$t$ [h]','Interpreter','latex');
legend('tru','est','EKF');
legend('Interpreter','latex');

ax(2) = nexttile();
plot(tsim / 3600,[vbatt; vb_sim; vb_EKF]); box off;
ylabel('$v_b$ [V]','Interpreter','latex');
xlabel('$t$ [h]','Interpreter','latex');

linkaxes(ax,'x')

h = findall(gcf,'Type','Axes');
set(h,'TickLabelInterpreter','latex');
%% performing an optimizaton 

% measurment signals 
mesh.qbatt = qbatt;

% initial exit flag
exitflag = 0;

% initial guess
X0 = Guess; 

% otimization optionss
options = optimoptions('lsqnonlin','Display',...
    'iter','FunctionTolerance',1e-4);

% normal EKF
tStart = tic;
while ~exitflag
    [x,resnorm,residual,exitflag,output] = ...
        lsqnonlin(@(x) EKF_batt_opt(x,mesh,init,modelEKF), ...
        X0, zeros(size(X0)) + eps, [], options);
    X0 = x; % homotopy :)
end
tELF(2) = toc(tStart);

modelEKF.Q = x(1);
modelEKF.R = x(2);
init.P0 = x(3);

xhat_opt_EKF = EKF_batt(modelEKF,init,mesh);

%% plots -- optimal results
figure(12); clf; set(gcf,'WindowStyle','docked');

plot(tsim / 3600,[qbatt; qb_sim; xhat_opt_EKF] / params.qn * 1e2); box off;
ylabel('$q_b$ [\%]','Interpreter','latex');
xlabel('$t$ [h]','Interpreter','latex');
legend('tru','est','EKF');
legend('Interpreter','latex');

h = findall(gcf,'Type','Axes');
set(h,'TickLabelInterpreter','latex');

%% FUNCTION: EKF
function xhat = EKF_batt(model,init,mesh)

  M = length(mesh.y);
  P = init.P0;
  n = numel(init.x0);
  
  xhat = zeros(n,M); 
  
  % initialize first two samples
  xhat(1) = init.x0;
  % xhat(2) = model.f(xhat(1), mesh.u(1));

  for ii = 1:M-1
    % Measurement update
    Ht = model.hx(xhat(ii));
    Kt = P * Ht' / (Ht * P * Ht' + model.R);
    xhat(ii) = xhat(ii) + Kt * (mesh.y(ii) - model.h(xhat(ii), mesh.u(ii)));
    P = P - Kt * Ht * P;

    % Time update
    xhat(ii+1) = model.f(xhat(ii), mesh.u(ii));
    Ft = model.fx(xhat(ii));
    P = Ft * P * Ft' + model.Q; %%%% assuming that G is 1
  end
end

%% FUNCTION: EKF - optimization

function residual = EKF_batt_opt(x,mesh,init,model)
    
  model.Q = x(1);
  model.R = x(2);
  init.P0 = x(3);
  
  xhat = EKF_batt(model,init,mesh);

  residual = mesh.qbatt - xhat;
end

