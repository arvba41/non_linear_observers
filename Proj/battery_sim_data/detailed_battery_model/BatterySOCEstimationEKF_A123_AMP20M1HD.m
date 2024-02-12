clear; clc; 

%% loading datasheet data

datasheet_data = Battery_A123_AMP20M1HD_data();

% extracting the battery data at 25 degC 

vocv_ds = datasheet_data.OCV_vec(:,4); 
qsoc_ds = datasheet_data.SOC_vec;
q0_ds   = datasheet_data.AH;
R0_ds   = datasheet_data.R0_vec(:,4);

params.ocv  = SplineInterpol(qsoc_ds / 100 * q0_ds, vocv_ds);
params.R0   = SplineInterpol(qsoc_ds / 100 * q0_ds, R0_ds);
params.R0_v = median(R0_ds); %% See the figure 10 and it will make sense
params.q0   = q0_ds;

% plotting OCV curves
qsweep = linspace(0,1,1000) * q0_ds;

figure(10); clf; set(gcf,'WindowStyle','docked');

subplot(2,2,1:2); 
plot(qsweep / q0_ds * 100, params.ocv.f(qsweep)); hold on;
plot(qsoc_ds, vocv_ds, 'x'); box off;
ylabel('$v_b$ [V]','Interpreter','latex');
xlabel('SOC [\%]','Interpreter','latex');

subplot(223); 
plot(qsweep / q0_ds * 100, params.ocv.df(qsweep)); box off;
ylabel('$v_b$ [V]','Interpreter','latex');
xlabel('SOC [\%]','Interpreter','latex');

subplot(224); 
plot(qsweep / q0_ds * 100, params.R0.f(qsweep) * 1e3); hold on;
plot(qsoc_ds, R0_ds * 1e3, 'x'); box off; 
yline(params.R0_v * 1e3,'k--',{'Median'},'Interpreter','latex')
ylabel('$R_b$ [m$\Omega$]','Interpreter','latex');
xlabel('SOC [\%]','Interpreter','latex');

h = findall(gcf,'Type','Axes');
set(h,'TickLabelInterpreter','latex');

clear vocv_ds qsoc_ds qsweep q0_ds h R0_ds % clearing data
%% simulating the data in simulink
simdata = sim("sim_mld_batt_A123AMP20M1HD.slx");

%% plotting simulation results
tsim    = simdata.tout;
soc     = simdata.Real_SOC;
ibatt   = simdata.Current;
vbatt   = simdata.Voltage;

figure(11); clf; set(gcf,'WindowStyle','docked');

subplot(311); 
plot(tsim / 3600, soc * 100); box off;
ylabel('SOC [\%]','Interpreter','latex');
xlabel('$t$ [h]','Interpreter','latex');

subplot(312); 
plot(tsim / 3600, ibatt); box off;
ylabel('$i_b$ [A]','Interpreter','latex');
xlabel('$t$ [h]','Interpreter','latex');

subplot(313); 
plot(tsim / 3600, vbatt); box off;
ylabel('$v_b$ [V]','Interpreter','latex');
xlabel('$t$ [h]','Interpreter','latex');

h = findall(gcf,'Type','Axes');
set(h,'TickLabelInterpreter','latex');

clear h ibatt soc vbatt tsim Ts initialSOC 

params.Ts = 1;
%% Battery model (simple R-model)

focv = @(q) params.ocv.f(q);   % battery ocv function of capacity
% Ro = @(q) params.R0.f(1);    % battery internal resistnace function of capacity

model.f = @(x,u) - u / 3600;                % state equation
model.h = @(x,u) focv(x) - u * params.R0_v; % measurement equation

init.q0 = 0.3 * params.q0;

%% Coulumb Counting (ClbCnt)

% measured data 
ibatt   = simdata.Current;
tsim    = simdata.tout;

N = length(tsim); 
h = mean(diff(tsim)); % sample time

% initializing vectors
qb_sim = zeros(size(tsim));
vb_sim = zeros(size(tsim));
qb_sim(1) = init.q0; 
vb_sim(1) = model.h(qb_sim(1),ibatt(1));

% euler forward simulation 
for ii = 2:N
    qb_sim(ii)  = qb_sim(ii-1) + h * model.f(qb_sim(ii-1),ibatt(ii)); 
    vb_sim(ii)  = model.h(qb_sim(ii),ibatt(ii));
end

figure(12); clf; set(gcf,'WindowStyle','docked');

subplot(311); 
plot(simdata.tout / 3600, simdata.Real_SOC * 100); hold on
plot(tsim / 3600, qb_sim / params.q0 * 100); box off;
ylabel('SOC [\%]','Interpreter','latex');
xlabel('$t$ [h]','Interpreter','latex');
legend('Real','Clb. Cnt.')
legend('Interpreter','latex','NumColumns',2);

subplot(312); 
plot(tsim / 3600, ibatt); box off;
ylabel('$i_b$ [A]','Interpreter','latex');
xlabel('$t$ [h]','Interpreter','latex');

subplot(313); 
plot(tsim / 3600, simdata.Voltage); hold on
plot(tsim / 3600, vb_sim); box off;
ylabel('$v_b$ [V]','Interpreter','latex');
xlabel('$t$ [h]','Interpreter','latex');

h = findall(gcf,'Type','Axes');
set(h,'TickLabelInterpreter','latex');

ClbCnt.tsim     = tsim;
ClbCnt.Current  = ibatt;
ClbCnt.SOC      = qb_sim / params.q0 * 100;
ClbCnt.Voltage  = vb_sim;

clear h vb_sim qb_sim tsim h N ibatt ii % cleaning data

%% EKF 

Ts = mean(diff(simdata.tout)); 
model.f  = @(x,u) x + Ts * ( -u / 3600); % making model descrete time

model.fx = @(x) 1;                  % state linearization 
model.hx = @(q) params.ocv.df(q);   % measurement linearization 

% tuning parameters
model.Q = 0.01;
model.R = 1;
init.P0 = 1;

% initial state 
init.x0 = init.q0;

% measurement 
mesh.y  = simdata.Voltage;
mesh.u  = simdata.Current;

% EKF filter
M = length(mesh.y);
P = init.P0;
n = numel(init.x0);

% initialize first sample
xhat = zeros(n,M); % estimated state
xhat(1) = init.x0;

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

qb_EKF = xhat;

% updating measuremet with the new estimate
vb_EKF = model.h(qb_EKF, simdata.Current');

figure(13); clf; set(gcf,'WindowStyle','docked');

subplot(311); 
plot(simdata.tout / 3600, simdata.Real_SOC * 100); hold on
plot(simdata.tout / 3600, qb_EKF / params.q0 * 100); box off;
ylabel('SOC [\%]','Interpreter','latex');
xlabel('$t$ [h]','Interpreter','latex');
legend('Real','EKF')
legend('Interpreter','latex','NumColumns',2);

subplot(312); 
plot(simdata.tout / 3600, simdata.Current); box off;
ylabel('$i_b$ [A]','Interpreter','latex');
xlabel('$t$ [h]','Interpreter','latex');

subplot(313); 
plot(simdata.tout / 3600, simdata.Voltage); hold on
plot(simdata.tout / 3600, vb_EKF); box off;
ylabel('$v_b$ [V]','Interpreter','latex');
xlabel('$t$ [h]','Interpreter','latex');

h = findall(gcf,'Type','Axes');
set(h,'TickLabelInterpreter','latex');

EKF.tsim     = simdata.tout;
EKF.Current  = simdata.Current';
EKF.SOC      = qb_EKF / params.q0 * 100;
EKF.Voltage  = vb_EKF;

clear h qb_EKF vb_EKF tsim h N ii P n Ht Ft Kt M Ts xhat % cleaning data

%% Plot -- Error

soc_error = sqrt((simdata.Real_SOC' * 100 - EKF.SOC).^2);

figure(14); clf; set(gcf,'WindowStyle','docked'); 
semilogy(simdata.tout / 3600, soc_error); box off; hold on;
yline(1,'k--',{'$\pm$ 1\%'},'Interpreter','latex','LabelHorizontalAlignment','left');
yline(median(soc_error),'k-',{['RMS: ' num2str(round(median(soc_error),1)) '\%']},'Interpreter','latex','LabelHorizontalAlignment','left');
title('\bf Root Square Error','Interpreter','latex');
ylabel('$\Delta$SOC [\%]','Interpreter','latex')
xlabel('$t$ [h]','Interpreter','latex')


h = findall(gcf,'Type','Axes');
set(h,'TickLabelInterpreter','latex');

clear h soc_error % cleaning data



