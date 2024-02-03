clear; clc; 

%% Simulate the model

params.N = 50;
params.x0 = 0.1;

params.w_var = 10;
params.v_var = 1;

% simulation without noise
noise = false;
[~, data_tru] = sim_nlm(params,noise);

% plot simulation data without noise
% figure(10); clf; set(gcf,"WindowStyle",'docked');
% tiledlayout(1,2);
% 
% nexttile();
% plot(1:params.N,data_tru.x); box off;
% ylabel('$x$','Interpreter','latex');
% xlabel('$k$','Interpreter','latex');
% 
% nexttile();
% plot(1:params.N,data_tru.y); box off;
% ylabel('$y$','Interpreter','latex');
% xlabel('$k$','Interpreter','latex');

% simulation with noise
noise = true;
[model, data] = sim_nlm(params,noise);

% plot simulation data with noise
figure(10); clf; set(gcf,"WindowStyle",'docked');
tiledlayout(1,2);

nexttile();
plot(1:params.N,data.x); box off; hold on
scatter(1:params.N,data_tru.x,'kx');
ylabel('$x$','Interpreter','latex');
xlabel('$k$','Interpreter','latex');
ylim([-40 40])
legend('mesh','true'); legend('Interpreter','latex');

nexttile();
plot(1:params.N,data.y); box off; hold on
scatter(1:params.N,data_tru.y,'kx');
ylabel('$y$','Interpreter','latex');
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

%% EKF estimation

% EKF_inputs
modelEKF.f = @(x,k) 0.5 * x + 25 * x / (1 + x^2) + 8 * cos(1.2 * ...
    (k - 1));
modelEKF.h = @(x,v) x^2 / 20;
modelEKF.fx = @(x,k) 0.5 + 25 / (x^2 + 1) - (50 * x^2) / (x^2 + 1)^2;
modelEKF.hx = @(x,v) x / 10;
modelEKF.Q = params.w_var;
modelEKF.R = params.v_var;

% initial guesses
init.x0 = params.x0;
init.P0 = 1;

xhatEKF = EKF(modelEKF, init, data); 

% plot EKF estimate with simulation data with noise
figure(11); clf; set(gcf,"WindowStyle",'docked');

plot(1:params.N,xhatEKF); box off; hold on
scatter(1:params.N,data_tru.x,'kx');
ylabel('$x$','Interpreter','latex');
xlabel('$k$','Interpreter','latex');
ylim([-40 40])
legend('EKF','true'); legend('Interpreter','latex');

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
%% EKF estimate -- LSQ optimization

lsqoptions = optimoptions("lsqnonlin","Display","iter-detailed");
x = lsqnonlin(@(x) EKF_optim(x,params,data,data_tru), 2, 0, [], lsqoptions);

init.P0 = x;
xhatEKF_opt = EKF(modelEKF, init, data); 

% plot EKF estimate with simulation data with noise
figure(12); clf; set(gcf,"WindowStyle",'docked');

plot(1:params.N,xhatEKF); box off; hold on
plot(1:params.N,xhatEKF_opt,'--');
scatter(1:params.N,data_tru.x,'kx');
ylabel('$x$','Interpreter','latex');
xlabel('$k$','Interpreter','latex');
legend('EKF','EKF-opt','true'); legend('Interpreter','latex');
ylim([-40 40])

h = findall(gcf,'Type','Axes');
set(h,'TickLabelInterpreter','latex');

% saving plots
golden_ratio = (1 + sqrt(5));
textheight = textwidth / golden_ratio;
figsize = [textwidth, textheight];

% Set size and no crop
set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', figsize);
set(gcf, 'PaperUnits', 'normalized', 'PaperPosition', [0, 0, 1, 1]);

% print -dpdf ../doc/figures/ex4_c_EKF_opt.pdf
%% Particle filter 

N = 500; % particle filter number of samples

xp = params.x0 + sqrt(params.w_var)*randn(1,N); % initial 500 samples

x = params.x0; % initial guess
y_mes = x^2 / 20 + sqrt(params.v_var) * randn; % measurement equation with noise
% state estimates
xhatPF_mean = x; 
xhatPF_max = x;


% For visualization purposes
[fxk(1,:),xk(1,:)] = ksdensity(xp);

for ii = 2:params.N
% time update
    % state equation
    x = 0.5 * x + 25 * x / (1 + x^2) + 8 * cos(1.2 * (ii - 1)) + ...
        sqrt(params.w_var) * randn;
    % measurement  equation
    y = x^2 / 20 + params.v_var * randn;
    % update xp
    xp_nxt = 0.5 * xp + 25 * xp ./ (1 + xp.^2) + 8 * cos(1.2 * (ii - 1)) + ...
            sqrt(params.w_var) * randn(1,N); 
% measurement update
    % update y
    y_nxt = xp.^2 / 20 + sqrt(params.v_var)*randn; 
    % residual generation
    e_res = y - y_nxt;
    % updating weights
    w = normpdf(e_res, 0, params.v_var);
    % normalizing weights
    w = w / sum(w);
    % resampling
    idx = resample(w);
    xp = xp_nxt(idx);
% Estimate
    xhatPF_mean(ii) = mean(xp);
    w = ones(size(w)) / sum(w); % re-distributiong equal weights
    
    [fxk(ii,:),xk(ii,:)] = ksdensity(xp);
    % fxr = ksdensity(xp);
    % [~,idx] = max(fxr);
    % visualization of the particles
    % if ii == round(params.N / 2) || ii == params.N
    %     nexttile();
    %         stem(xp,w); box off
    %         ylabel('$w$','Interpreter','latex');
    %         xlabel('$x$','Interpreter','latex');
    % 
    %     nexttile();
    %         [fxk,xk] = ksdensity(xp);
    %         plot(xk,fxk); box off
    %         ylabel('$f(x)$','Interpreter','latex');
    %         xlabel('$x$','Interpreter','latex');
    % end

    % xhatPF_max(ii) = xp(idx);
    % xhatPF_max(ii) = max(xp);
end

% plot EKF estimate with simulation data with noise
figure(13); clf; set(gcf,"WindowStyle",'docked');

plot(1:params.N,xhatPF_mean); box off; hold on
% plot(1:params.N,xhatPF_max,'--');
scatter(1:params.N,data_tru.x,'kx');
ylabel('$x$','Interpreter','latex');
xlabel('$k$','Interpreter','latex');
ylim([-40 40])
legend('PF','true'); legend('Interpreter','latex');

h = findall(gcf,'Type','Axes');
set(h,'TickLabelInterpreter','latex');

% saving plots
golden_ratio = (1 + sqrt(5));
textheight = textwidth / golden_ratio;
figsize = [textwidth, textheight];

% Set size and no crop
set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', figsize);
set(gcf, 'PaperUnits', 'normalized', 'PaperPosition', [0, 0, 1, 1]);

% print -dpdf ../doc/figures/ex4_c_PF.pdf

%% Visualizing the particle filter
figure(100); set(gcf,'WindowStyle','docked'); clf;
    for ii = 1:params.N
        p = plot3(ii*ones(1,100),xk(ii,:),fxk(ii,:)); hold on; grid on;
        p.Color = [1 0 0] + ([-1 1 1] / params.N * (ii - 1));
    end
    zlabel('$f(x)$','Interpreter','latex');
    ylabel('$x$','Interpreter','latex');
    xlabel('$N$','Interpreter','latex');

    findall(gcf,'Type','Axes');
    h = set(gca,'TickLabelInterpreter','latex');
    
%% Quantitative comparison of the methods 

SQE_X_EKF = sum((data_tru.x - xhatEKF).^2);
SQE_X_PF = sum((data_tru.x - xhatPF_mean).^2);

figure(101); clf; set(gcf,'WindowStyle','docked');
Y = [SQE_X_EKF, SQE_X_PF] / 1e3;
bar(Y)
text(1:length(Y),Y,num2str(round(Y)'),'vert','bottom','horiz','center','Interpreter','latex'); 
box off

ylabel('SQE [$\times10^3$]','Interpreter','latex');
xticklabels({'EKF','PF'});

    findall(gcf,'Type','Axes');
    h = set(gca,'TickLabelInterpreter','latex');
%% %%%%
% FUNCTIONS
%%% %%%%
%% simulating the equation
function [model, data] = sim_nlm(params,noise)

N = params.N; % sample time 

% Model equations
model.f = @(x,k,w) 0.5 * x + 25 * x / (1 + x^2) + 8 * cos(1.2 * ...
    (k - 1)) + w;
% Measurement equations
model.h = @(x,v) x^2 / 20 + v;

% adding noise
if noise
    wn = randn(1,N) * sqrt(params.w_var);
    vn = randn(1,N) * sqrt(params.v_var);
else
    wn = zeros(1,N);
    vn = zeros(1,N);
end

% initializaiton
data.x(:,1) = 0.01;
data.y(:,1) = model.h(data.x(:,1),vn(1));

% simulation 
for ii = 2:N
    data.x(:,ii) = model.f(data.x(:,ii-1),ii,wn(ii)); % update state
    data.y(:,ii) = model.h(data.x(:,ii),vn(ii)); % update mesh
end

end


%% EKF 
function xhat = EKF(model,init,data)

  M = length(data.y);
  P = init.P0;
  n = numel(init.x0);
  
  xhat = zeros(n,M); xhat(:,1) = init.x0;
  
  for ii = 1:M-1
    % Measurement update
    Ht = model.hx(xhat(:,ii));
    Kt = P * Ht' / (Ht * P * Ht' + model.R);
    xhat(:,ii) = xhat(:,ii) + Kt * (data.y(:, ii) - model.h(xhat(:, ii)));

    % Time update
    xhat(:,ii+1) = model.f(xhat(:, ii), ii + 1);
    Ft = model.fx(xhat(:, ii));
    P = Ft * P * Ft' + model.Q; %%%% assuming that G is 1
  end
end

%% EKF-optim residual generation
function residual = EKF_optim(x,params,data,data_tru)
    % EKF_inputs
    modelEKF.f = @(x,k) 0.5 * x + 25 * x / (1 + x^2) + 8 * cos(1.2 * ...
        (k - 1));
    modelEKF.h = @(x,v) x^2 / 20;
    modelEKF.fx = @(x,k) 0.5 + 25 / (x^2 + 1) - (50 * x^2) / (x^2 + 1)^2;
    modelEKF.hx = @(x,v) x / 10;
    modelEKF.Q = params.w_var;
   modelEKF.R = params.v_var;
    
    % initial guesses
    init.x0 = params.x0;
    init.P0 = x;
    
    xhatEKF = EKF(modelEKF, init, data); 
    
    residual = xhatEKF - data_tru.x;
end 