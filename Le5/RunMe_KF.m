clear all; clc; saveplots = false; add_casadi_path = false;
% Ensuring casadi is added to the path 
if add_casadi_path 
    % add your casadi path here
    casadi_path = '~/PhDProject/matlab_functions/casadi-3.6.4-linux64-matlab2018b/';
    addpath(casadi_path);
end
import casadi.*

%% Simulating the model

N = 80; % number of samples

% model 
model.A = [0.99 0.2;-0.1 0.3];
model.B = [0; 1];
model.C = [1 -3];
model.R = 0.01;
model.Q = model.B * 1 * model.B';

% initial values
x0 = [0;0]; 
data.x = zeros(length(x0),N);
data.y = zeros(1,N);
data.x(:,1) = model.A * x0 + model.B * sqrt(1) * rand;
data.y(:,1) = model.C * data.x(:,1) + sqrt(model.R) * rand;
for ii = 2:N
    data.x(:,ii) = model.A * data.x(:,ii-1) + model.B * rand;
    data.y(:,ii) = model.C * data.x(:,ii) + sqrt(model.R) * rand;
end

figure(10); clf; set(gcf,"WindowStyle",'docked');
tiledlayout(1,2);

nexttile();
plot(1:N,data.x(1,:)); box off; hold on
ylabel('$x_1$','Interpreter','latex');
xlabel('$k$','Interpreter','latex');
% legend('meas','true'); legend('Interpreter','latex');

nexttile();
plot(1:N,data.x(2,:)); box off; hold on
ylabel('$x_2$','Interpreter','latex');
xlabel('$k$','Interpreter','latex');

h = findall(gcf,'Type','Axes');
set(h,'TickLabelInterpreter','latex');

% saving plots
textwidth = 14;
golden_ratio = (1 + sqrt(5));
textheight = textwidth / golden_ratio;
figsize = [textwidth, textheight];

% Set size and no crop
set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', figsize);
set(gcf, 'PaperUnits', 'normalized', 'PaperPosition', [0, 0, 1, 1]);

if saveplots
    print -dpdf ../doc/figures/ex5_MHE_sim.pdf
end
%% Estimation using KF

% initial guesses
init.x0 = x0;
init.P0 = diag([1, 1]);

xhatKF = KF(model, init, data); 

% plot KF estimate with simulation data 
figure(11); clf; set(gcf,"WindowStyle",'docked');
tiledlayout(1,2);

nexttile();
plot(1:N,data.x(1,:)); box off; hold on
plot(1:N,xhatKF(1,:)); box off; hold on
legend('sim','KF'); legend('Interpreter','latex');
ylabel('$x_1$','Interpreter','latex');
xlabel('$N$','Interpreter','latex');

nexttile();
plot(1:N,data.x(2,:)); box off; hold on
plot(1:N,xhatKF(2,:)); box off; hold on
legend('sim','KF'); legend('Interpreter','latex');
ylabel('$x_1$','Interpreter','latex');
xlabel('$N$','Interpreter','latex');

h = findall(gcf,'Type','Axes');
set(h,'TickLabelInterpreter','latex');

% saving plots
textwidth = 14.9;
golden_ratio = (1 + sqrt(5));
textheight = textwidth / golden_ratio;
figsize = [textwidth, textheight];

% Set size and no crop
set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', figsize);
set(gcf, 'PaperUnits', 'normalized', 'PaperPosition', [0, 0, 1, 1]);

if saveplots
    print -dpdf ../doc/figures/ex5_MHE_KF.pdf
end

%% Estimation using MHE 
% initialization
xhat_MHE = zeros(length(x0),N);
% xhat_MHE(:,1) = 0;
xhat_MHE(:,1) = data.x(:,1);

T = 2; % minimum time for MHE to begin

while T <= N % run untill final time 
    xhat_MHE(:,T) = MHE(model,init,data,T); % estimate stat using MHE 
    T = T + 1; % increase time 
end

%% plot MHE estimate with simulation data 
figure(12); clf; set(gcf,"WindowStyle",'docked');
tiledlayout(1,2);

nexttile();
plot(1:N,data.x(1,:)); box off; hold on
plot(1:N,xhatKF(1,:)); box off; hold on
plot(1:N,xhat_MHE(1,:)); box off; hold on
ylabel('$x_1$','Interpreter','latex');
xlabel('$N$','Interpreter','latex');

nexttile();
plot(1:N,data.x(2,:)); box off; hold on
plot(1:N,xhatKF(2,:)); box off; hold on
plot(1:N,xhat_MHE(2,:)); box off; hold on
L = legend('sim','KF','MHE'); 
L = legend('Interpreter','latex','NumColumns',3);
L.Layout.Tile = 'north';
ylabel('$x_2$','Interpreter','latex');
xlabel('$N$','Interpreter','latex');

h = findall(gcf,'Type','Axes');
set(h,'TickLabelInterpreter','latex');

% saving plots
textwidth = 14.9;
golden_ratio = (1 + sqrt(5));
textheight = textwidth / golden_ratio;
figsize = [textwidth, textheight];

% Set size and no crop
set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', figsize);
set(gcf, 'PaperUnits', 'normalized', 'PaperPosition', [0, 0, 1, 1]);

if saveplots
    print -dpdf ../doc/figures/ex5_MHE.pdf
    % print -dpdf ../doc/figures/ex5_MHE_10T.pdf
    % print -dpdf ../doc/figures/ex5_MHE_20T.pdf
end

%% MSE
MSE.KF = sum((data.x - xhatKF).^2, 2);
MSE.MHE = sum((data.x - xhat_MHE).^2, 2);

%% FUNCTION: KF
function xhat = KF(model,init,data)

  M = length(data.y);
  P = init.P0;
  n = numel(init.x0);
  
  xhat = zeros(n,M); xhat(:,1) = init.x0;
  
  for ii = 1:M-1
    % Measurement update
    Ht = model.C;
    Kt = P * Ht' / (Ht * P * Ht' + model.R);
    xhat(:,ii) = xhat(:,ii) + Kt * (data.y(ii) - model.C * xhat(:, ii));
    P = P - Kt * Ht * P;

    % Time update
    xhat(:,ii+1) = model.A * xhat(:, ii);
    Ft = model.A;
    P = Ft * P * Ft' + model.Q; %%%% assuming that G is 1
  end

end
%% FUNCTION: MHE % full-optimization (very time consuming)
function xhat = MHE(model,init,data,T)

% fetching data model information
A = model.A;
B = model.B;
C = model.C;
Q = model.Q;
R = model.R;

% fetching output 
y = data.y;

% optimization using casADi
opti = casadi.Opti();

x0_hat  = opti.variable(2, 1);  % initial state
x0      = init.x0;              % True initial state

P0      = init.P0;  % initial covariance matrix

w = opti.variable(1, T-1);      % noise variables

X = opti.variable(2, T);        % state vector
opti.subject_to(X(:, 1) == x0_hat);  % set initial state

vk = opti.variable(1, T-1);     % output error vector (y - yk)
% covariance matrices
R_inv = diag(ones(1, T-1)) / R; % measurement deviation penalty
Q_inv = diag(ones(1, T-1)) / 1; % states noise deviation penalty

% simulating the states (multiple-shooting)
for ii = 1:1:T-1
    opti.subject_to(X(:,ii+1) == A * X(:, ii) + B * w(ii)); 
    vk(ii) = y(ii) - C*X(:, ii);
end

opti.subject_to(0 <= w);    % noise to be grater than 0

% Cost funciont
J = vk * R_inv * vk' + w * Q_inv * w' + ...
    (x0_hat - x0)' / P0 * (x0_hat - x0);

opti.minimize(J);       % Onjective

opti.solver('ipopt');  
sol = opti.solve(); 

% fetching the solutions from thre OCP
    w_op = sol.value(w);
    x_temp = sol.value(x0_hat);
% simulating the system
for ii = 1:1:T-1
    x_temp = A*x_temp +B*w_op(ii);
end

xhat = x_temp; % taking final value of x at time T
end