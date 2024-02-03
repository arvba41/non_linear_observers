clear; clc; 

load("battery_sim_data/ocv_data.mat");
load("battery_sim_data/cyc_data.mat");

%% Coulumb counting

% models parameters
qocv = data.qbatt;
vocv = data.vbatt;
params.Ro = 1.5e-3; % internal resistance (estimate)
params.qn = 24;
params.q0 = params.qn * 0.7;

% measurements
q0_true = params.qn * cyc_data.Qinit;
tsim = cyc_data.tprc;
ibatt = cyc_data.ibatt_prc;
vbatt = cyc_data.vbatt_prc;
qbatt = cyc_data.qbatt_prc / 100 * params.qn;

% battery model (pure simulation with reasonable parameters)
fvocv = @(q) interp1(qocv,vocv,q,'spline','extrap'); % battery SOC-ocv function
Ro = params.Ro;                 % battery internal resistnace 
x_qb = @(x,u) - u / 3600;         % state equation
y_vb = @(x,u) fvocv(x) - u * Ro; % measurement equation

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

h = findall(gcf,'Type','Axes');
set(h,'TickLabelInterpreter','latex');

%% EKF method
Ts = mean(diff(tsim)); % time step 

Guess = [0.001 5 50];

% model
modelEKF.f = @(x,u) x + Ts * (- u / 3600);
% measurement 
modelEKF.h = @(x,u) fvocv(x) - u * Ro;
% linearization f
modelEKF.fx = @(x) 1;
% linearization h
modelEKF.hx = @(x,x_p) (fvocv(x) - fvocv(x_p)) / (x - x_p); % numeric differentiation on focv

% turing parameters          
modelEKF.Q = diag(Guess(1)); % model covriance matrix 
modelEKF.R = diag(Guess(2)); % measurement covriance matrix
init.P0 = diag(Guess(3));  % Initial covariance

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
%% FUNCTION: EKF
function xhat = EKF_batt(model,init,mesh)

  M = length(mesh.y);
  P = init.P0;
  n = numel(init.x0);
  
  xhat = zeros(n,M); 
  
  % initialize first two samples
  xhat(1) = init.x0;
  xhat(2) = model.f(xhat(1), mesh.u(1));

  for ii = 2:M-1
    % Measurement update
    Ht = model.hx(xhat(ii),xhat(ii-1));
    Kt = P * Ht' / (Ht * P * Ht' + model.R);
    xhat(ii) = xhat(ii) + Kt * (mesh.y(ii) - model.h(xhat(ii), mesh.u(ii)));
    
    % Time update
    xhat(ii+1) = model.f(xhat(ii), mesh.u(ii));
    Ft = model.fx(xhat(ii));
    P = Ft * P * Ft' + model.Q; %%%% assuming that G is 1
  end
end


