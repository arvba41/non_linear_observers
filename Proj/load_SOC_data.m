clear; clc; 

%% loading data
% addpath('batt_online_data/SP1_25C_LC_OCV_11_5_2015/');

data = readmatrix("batt_online_data/SP1_25C_LC_OCV_11_5_2015/11_5_2015_low current OCV test_SP20-1.xlsx");

% rmpath('batt_online_data/SP1_25C_LC_OCV_11_5_2015/');
%% data segrigation

ibatt = data(:,4);
tmesh = data(:,1);
vbatt = data(:,3) / 1e3;

q0 = 2152.38;
qbatt = cumtrapz(tmesh,ibatt) / 3600 + q0;

%% plots

figure(10); clf; set(gcf,'WindowStyle','docked');

tiledlayout(3,2)

nexttile(1);
plot(tmesh,ibatt);
ylabel('$i_b$ [mA]','Interpreter','latex');
xlabel('$t$ [s]','Interpreter','latex');

nexttile(3);
plot(tmesh,vbatt);
ylabel('$v_b$ [V]','Interpreter','latex');
xlabel('$t$ [s]','Interpreter','latex');

nexttile(5);
plot(tmesh,qbatt);
ylabel('$q_b$ [Ah]','Interpreter','latex');
xlabel('$t$ [s]','Interpreter','latex');

nexttile(2,[3 1]);
plot(qbatt,vbatt);
ylabel('$v_b$ [V]','Interpreter','latex');
xlabel('$q_b$ [Ah]','Interpreter','latex');

h = findall(gcf,'Type','Axes');
set(h,'TickLabelInterpreter','latex');