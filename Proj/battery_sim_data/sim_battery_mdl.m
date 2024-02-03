clear; clc; 

% load the simulation parameters
Qo = 24;
C_rate = 1/100;
I = Qo * C_rate;
cyc_time = Qo / I * 3600;

% simulating the simulink model 
simdata = sim("simple_battery_mdl_pulses.slx","StopTime",num2str(cyc_time*4));

% data collections
tsim    = simdata.tout;
ibatt   = simdata.ibatt;
vbatt   = simdata.vbatt;
SOC     = simdata.soc;

qbatt   = cumtrapz(tsim,-ibatt) / 3600 + Qo; % battery capacity (current integration)
%% plots -- simulation results
figure(10); clf; set(gcf,'WindowStyle','docked');

tiledlayout(3,2)

nexttile(1);
plot(tsim,ibatt*1e3); box off;
ylabel('$i_b$ [mA]','Interpreter','latex');
xlabel('$t$ [s]','Interpreter','latex');

nexttile(3);
plot(tsim,vbatt); box off;
ylabel('$v_b$ [V]','Interpreter','latex');
xlabel('$t$ [s]','Interpreter','latex');

nexttile(5);
plot(tsim,qbatt); box off;
ylabel('$q_b$ [Ah]','Interpreter','latex');
xlabel('$t$ [s]','Interpreter','latex');

nexttile(2,[3 1]);
plot(qbatt,vbatt); box off;
ylabel('$v_b$ [V]','Interpreter','latex');
xlabel('$q_b$ [Ah]','Interpreter','latex');

h = findall(gcf,'Type','Axes');
set(h,'TickLabelInterpreter','latex');
%% determining the modes
mode.idx_crg    = find(ibatt > 0); % charging index
mode.idx_dcrg   = find(ibatt < 0); % discharging index

% capacity calculation
qbatt_crg   = qbatt(mode.idx_crg);
qbatt_dcrg  = qbatt(mode.idx_dcrg);
% ocv calculation
vbatt_crg   = vbatt(mode.idx_crg);
vbatt_dcrg  = vbatt(mode.idx_dcrg);

% ocv curve interpolation
qbatt_int = linspace(0,Qo,1000);
vbatt_crg_int = interp1(qbatt_crg,vbatt_crg,qbatt_int,'spline');
vbatt_dcrg_int = interp1(qbatt_dcrg,vbatt_dcrg,qbatt_int,'spline');
%% plots -- simulation ocv curves
figure(11); clf; set(gcf,'WindowStyle','docked');

tiledlayout(1,2), 

nexttile();
plot(qbatt_crg,vbatt_crg); box off; hold on;
plot(qbatt_int,vbatt_crg_int,'--');
ylabel('$v_b$ [V]','Interpreter','latex');
xlabel('$q_b$ [Ah]','Interpreter','latex');
title('Charging OCV','Interpreter','latex');

nexttile();
plot(qbatt_dcrg,vbatt_dcrg); box off; hold on;
plot(qbatt_int,vbatt_dcrg_int,'--');
ylabel('$v_b$ [V]','Interpreter','latex');
xlabel('$q_b$ [Ah]','Interpreter','latex');
title('Disharging OCV','Interpreter','latex');

h = findall(gcf,'Type','Axes');
set(h,'TickLabelInterpreter','latex');
%% battery ocv curve export
data.qbatt = qbatt_crg;
data.vbatt = vbatt_crg;
save('ocv_data','data');
%% simulating test cycle
Qo = 24;
C_rate = 0.8;
I = Qo * C_rate;
cyc_time = Qo / I * 3600;
Qinit = 80;
tend = 6 * 3600;

% simulating the simulink model 
simdata = sim("simple_battery_mdl_cycle.slx");

% data collections
tsim    = simdata.tout;
ibatt   = simdata.ibatt;
vbatt   = simdata.vbatt;
SOC     = simdata.soc;
%% plots -- simulation cycle
figure(12); clf; set(gcf,'WindowStyle','docked');

tiledlayout(1,3)

nexttile();
plot(tsim / 3600,ibatt); box off;
ylabel('$i_b$ [A]','Interpreter','latex');
xlabel('$t$ [h]','Interpreter','latex');

nexttile();
plot(tsim,vbatt); box off;
ylabel('$v_b$ [V]','Interpreter','latex');
xlabel('$t$ [s]','Interpreter','latex');

nexttile();
plot(tsim,SOC); box off;
ylabel('$q_b/q_0$ [\%]','Interpreter','latex');
xlabel('$t$ [s]','Interpreter','latex');

h = findall(gcf,'Type','Axes');
set(h,'TickLabelInterpreter','latex');
%% exporting simlated data
cyc_data.tsim = tsim;
cyc_data.ibatt  = ibatt;
cyc_data.vbatt  = vbatt;
cyc_data.soc    = SOC;
cyc_data.Qinit  = Qinit;
save('cyc_data','cyc_data');
%% processing data with noise and also incrasing the sample time to 100s
tprc = tsim(1):5:tsim(end);

n_ibatt = randn(size(tprc)) * sqrt(0.5) + 0.5 * sin(2 * pi * 50 * tprc);
n_vbatt = randn(size(tprc)) * sqrt(0.5e-3);
ibatt_prc = interp1(tsim,ibatt,tprc,'nearest') + n_ibatt;
vbatt_prc = interp1(tsim,vbatt,tprc,'linear') + n_vbatt;
qbatt_prc = interp1(tsim,SOC,tprc,'linear');
%% plots -- processed data - simulation cycle 
figure(13); clf; set(gcf,'WindowStyle','docked');

tiledlayout(1,2)

ax(1) = nexttile();
plot(tprc,ibatt_prc); box off;
ylabel('$i_b$ [A]','Interpreter','latex');
xlabel('$t$ [s]','Interpreter','latex');

ax(2) = nexttile();
plot(tprc,vbatt_prc); box off;
ylabel('$v_b$ [V]','Interpreter','latex');
xlabel('$t$ [s]','Interpreter','latex');

linkaxes(ax,'x')

h = findall(gcf,'Type','Axes');
set(h,'TickLabelInterpreter','latex');
%% exporting processed simlated data
cyc_data.tprc = tprc;
cyc_data.ibatt_prc = ibatt_prc;
cyc_data.vbatt_prc = vbatt_prc;
cyc_data.qbatt_prc = qbatt_prc;
save('cyc_data','cyc_data')