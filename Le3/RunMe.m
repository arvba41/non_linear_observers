clear; clc; close all; 

plots = false; % 
saveplots = false;

%% Simulate robot movement
% parameters 
params.Ts = 0.1;
params.V = 3;

[model, data, params] = sim_robot(params,plots);

%% EKF estimations 
% intelligent guess for covariance and intial covance matrices

Guess = [[0.05 0.05 0.05],...
         [0.08 0.2],...
         [0.01 0.1 0.01]];

% model
modelEKF.f = @(x,u) [x(1) + params.Ts * params.V * cos(x(3));...
                     x(2) + params.Ts * params.V * sin(x(3));...
                     x(3) + params.Ts * u];
% measurement 
modelEKF.h = @(x,u) [sqrt(x(1)^2 + x(2)^2);...
                     atan(x(2) / x(1))];
% linearization f
modelEKF.fx = @(x,u) [1 0 -params.Ts * params.V * sin(x(3));...
                      0 1 params.Ts * params.V * cos(x(3));...
                      0 0 1];
% linearization h                  
modelEKF.hx = @(x,u) [x(1) / sqrt(x(1)^2 + x(2)^2), ...
                            x(2) / sqrt(x(1)^2 + x(2)^2), ...
                            0;...
                      -x(2) / (x(1)^2 * (x(2)^2 / x(1)^2 + 1)), ...
                            1 / (x(1) * (x(2)^2 / x(1)^2 + 1)), ...
                            0];

% turing parameters 
modelEKF.Q = diag(Guess(1:3));
modelEKF.R = diag(Guess(4:5));
init.P0 = diag(Guess(6:8)); % Initial covariance

init.x0 = [10; 5; 0]; % Initial stafte estimate 

% normal EKF
tStart = tic; 
xhatEKF = EKF(modelEKF,init,data); 
tELF(1) = toc(tStart);

% EKF with SR-alg
tStart = tic; 
xhatEKF_SR = EKF_sr_new(modelEKF,init,data); 
tELF_SR(1) = toc(tStart);

%% performing an optimization

% initial exit flag
exitflag = 0;

% initial guess
x0 = Guess; 

% otimization optionss
options = optimoptions('lsqnonlin','Display',...
    'iter','FunctionTolerance',1e-6);

% normal EKF
tStart = tic;
while ~exitflag
    [x,resnorm,residual,exitflag,output] = ...
        lsqnonlin(@(x) EKF_res(x,data,init,modelEKF), ...
        x0, zeros(size(x0)) + eps, [], options);
    x0 = x; % homotopy :)
end
tELF(2) = toc(tStart);

%% running EKF with new values
modelEKF.Q = diag(x(1:3));
modelEKF.R = diag(x(4:5));
init.P0 = diag([x(6) x(7) x(8)]); % Initial covariance

xhatEKF_opt = EKF(modelEKF,init,data);

disp('-------------------------------------------------------------------');

%% resetting exitflag
exitflag = 0;
% resetting initial guess
x0 = Guess; 

% EKF with SR-alg
tStart = tic;
while ~exitflag
    [x_SR,resnorm,residual,exitflag,output] = ...
        lsqnonlin(@(x) EKF_SR_res(x,data,init,modelEKF), ...
        x0, zeros(size(x0)), [], options);
    x0 = x_SR; % homotopy :)
end
tELF_SR(2) = toc(tStart);

%% running EKF with new values
modelEKF.Q = diag(x_SR(1:3));
modelEKF.R = diag(x_SR(4:5));
init.P0 = diag([x_SR(6) x_SR(7) x_SR(8)]); % Initial covariance

xhatEKF_SR_opt = EKF_sr(modelEKF,init,data);

%% plots 
plot_path(20,data,xhatEKF,xhatEKF_opt,saveplots,'EKF_Le3');
plot_path(30,data,xhatEKF_SR,xhatEKF_SR_opt,saveplots,'EKF_Le3_SR');
plot_error(21,data,xhatEKF,xhatEKF_opt,saveplots,'EKF_Le3_err');
plot_error(31,data,xhatEKF_SR,xhatEKF_SR_opt,saveplots,'EKF_Le3_SR_err');

%% important functions

% residual for EKF optimization 
function residual = EKF_res(x,data,init,modelEKF)

    %%% running EKF 
    % tuning parameters
    modelEKF.Q = diag([x(1) x(2) x(3)]);
    modelEKF.R = diag([x(4) x(5)]);
    
    % P0 initialization
    init.P0 = diag([x(6) x(7) x(8)]); % Initial covariance
    
    Xhat = EKF(modelEKF,init,data);

    Xhat_sum = sum(Xhat,1);

    X = sum(data.x,1);
    residual = X - Xhat_sum;
end

% residual for EKF with SR-alg optimization 
function residual = EKF_SR_res(x,data,init,modelEKF)

    %%% running EKF 
    % tuning parameters
    modelEKF.Q = diag([x(1) x(2) x(3)]);
    modelEKF.R = diag([x(4) x(5)]);
    
    % P0 initialization
    init.P0 = diag([x(6) x(7) x(8)]); % Initial covariance
    
    Xhat = EKF_sr(modelEKF,init,data);

    Xhat_sum = sum(Xhat,1);
    
    X = sum(data.x,1);
    residual = X - Xhat_sum;
end

