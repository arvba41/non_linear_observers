clear all; clc;
% Ensure to add casadi to the path 
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
data.x(:,1) = model.A * x0 + sqrt(1)*rand(2,1);
data.y(:,1) = model.C * data.x(:,1) + sqrt(model.R)*randn(1);
for ii = 2:N
    data.x(:,ii) = model.A * data.x(:,ii-1) + model.B*rand;
    data.y(:,ii) = model.C * data.x(:,ii) + sqrt(model.R)*randn(1);
end

figure(10); clf; set(gcf,"WindowStyle",'docked');
tiledlayout(1,2);

nexttile();
plot(1:N,data.x(1,:)); box off; hold on
ylabel('$x_1$','Interpreter','latex');
xlabel('$k$','Interpreter','latex');
legend('mesh','true'); legend('Interpreter','latex');

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

% print -dpdf ../doc/figures/ex4_c_sim.pdf
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

% print -dpdf ../doc/figures/ex4_c_EKF.pdf

%% Estimation using MHE 
% initialization
xhat_MHE = zeros(length(x0),N);
xhat_MHE(:,1) = x0;

T = 10; % minimum time for MHE to begin

while T <= N % run untill final time 
    xhat_MHE(:,T) = MHE(model,init,data,T); % estimate stat using MHE 
    T = T + 1; % incriment time 
end

%% % plot KF estimate with simulation data 
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
legend('sim','KF','MHE'); 
legend('Interpreter','latex','NumColumns',3);
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

% print -dpdf ../doc/figures/ex4_c_EKF.pdf


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

vk = opti.variable(1, T-1); % output error vector
% covariance matrices
Q1 = diag(ones(1, T-1)) / R;    % measurement
Q2 = diag(ones(1, T-1));        % states

% simulating the states (multiple-shooting)
for ii = 1:1:T-1
    opti.subject_to(X(:,ii+1) == A * X(:, ii) + B * w(ii)); 
    vk(ii) = y(ii) - C*X(:, ii);
end

opti.subject_to(0 <= w);    % noise to be grater than 1

% Cost funciont
J = vk * Q1 * vk' + w * Q2 * w' + ...
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