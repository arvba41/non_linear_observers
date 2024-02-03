clear; 

%% system definition

% states 
A = [0 1; 0 0];
C = [1, 0];
system = @(t,y) A*y;

% simulation 
y0 = [0; 1];
% tspan = [0 6];
tspan = linspace(0,6,5000);
[t_true,y_true] = ode45(@(t,y) system(t,y),tspan,y0);
% adding noise
y_noise = y_true(:,1) + (rand(size(t_true)) - 0.5)*0.5;

% simulation of full-order observer
y0 = [1; 0];
y_val = @(t) interp1(t_true, y_true(:,1), t, 'spline');
[t_fo,y_fo] = ode45(@(t,y) full_order(y,y_val(t),A,C), tspan, y0);
y_temp = interp1(t_true, y_true, t_fo, 'spline');
e_fo = y_temp - y_fo;

% simulation of full-order observer with noise
y_valN = @(t) interp1(t_true, y_noise(:,1), t, 'spline');
[t_foN,y_foN] = ode45(@(t,y) full_order(y,y_valN(t),A,C), tspan, y0);
y_temp = interp1(t_true, y_true, t_foN, 'spline');
e_foN = y_temp - y_foN;

% simulation of reduced-order observer 
y0 = 0;
[t_ro,w] = ode45(@(t,y) reduced_order(y,y_val(t)),tspan,y0);
x2hat = w + 5 * y_val(t_ro);
y_temp = interp1(t_true, y_true(:,2), t_ro, 'spline');
e_ro = y_temp - x2hat;

% simulation of reduced-order observer (with noise)
[t_roN,w] = ode45(@(t,y) reduced_order(y,y_valN(t)),tspan,y0);
x2hatN = w + 5 * y_val(t_roN);
y_temp = interp1(t_true, y_true(:,2), t_roN, 'spline');
e_roN = y_temp - x2hatN;

%% plots
figure(1); clf; 

tiledlayout(2,2); 

ax(1) = nexttile; 
yyaxis left
plot(t_true,y_true(:,1)); box off;
hold on;
plot(t_fo,y_fo(:,1)); box off;
ylabel('$x_1$','Interpreter','latex')
yyaxis right
plot(t_true,y_true(:,2)); box off;
hold on;
plot(t_fo,y_fo(:,2)); box off;
ylabel('$x_2$','Interpreter','latex')
xlabel('$t$ [s]','Interpreter','latex')
title('\bf no noise','Interpreter','latex')

ax(2) = nexttile; 
yyaxis left
plot(t_true,y_true(:,1)); box off;
hold on;
plot(t_foN,y_foN(:,1)); box off;
ylabel('$x_1$','Interpreter','latex')
yyaxis right 
plot(t_true,y_true(:,2)); box off;
hold on;
plot(t_foN,y_foN(:,2)); box off;
ylabel('$x_2$','Interpreter','latex')
xlabel('$t$ [s]','Interpreter','latex')
title('\bf noise','Interpreter','latex')

ax(3) = nexttile; 
yyaxis left
plot(t_fo,e_fo(:,1)); box off;
ylabel('$\tilde x_1$','Interpreter','latex')
yline(0,'--')
yyaxis right
plot(t_fo,e_fo(:,2)); box off;
ylabel('$\tilde x_2$','Interpreter','latex')
xlabel('$t$ [s]','Interpreter','latex')
yline(0,'--')

ax(4) = nexttile; 
yyaxis left
plot(t_foN,e_foN(:,1)); box off;
ylabel('$\tilde x_1$','Interpreter','latex')
yline(0,'--')
yyaxis right
plot(t_foN,e_foN(:,2)); box off;
ylabel('$\tilde x_2$','Interpreter','latex')
xlabel('$t$ [s]','Interpreter','latex')
yline(0,'--')

% saving plots
textwidth = 14;
golden_ratio = (1 + sqrt(5)) / 2;
textheight = textwidth / golden_ratio;
figsize = [textwidth, textheight];

% Get latex font in ticks
h = findall(gcf,'Type','axes'); % An array if you have subplots
set(h, 'TickLabelInterpreter', 'latex')

% Set size and no crop
set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', figsize);
set(gcf, 'PaperUnits', 'normalized', 'PaperPosition', [0, 0, 1, 1]);

% print -dpdf ../doc/figures/ex6_p1.pdf
%% plots

figure(2); clf; 

tiledlayout(2,2); 

ax(1) = nexttile; 
plot(t_true,y_true(:,2)); box off;
hold on;
plot(t_ro,x2hat); box off;
ylabel('$x_2$','Interpreter','latex')
xlabel('$t$ [s]','Interpreter','latex')
title('\bf no noise','Interpreter','latex')

ax(2) = nexttile; 
plot(t_true,y_true(:,2)); box off;
hold on;
plot(t_roN,x2hatN); box off;
ylabel('$x_2$','Interpreter','latex')
xlabel('$t$ [s]','Interpreter','latex')
title('\bf noise','Interpreter','latex')

ax(3) = nexttile; 
plot(t_ro,e_ro); box off;
ylabel('$\tilde x_2$','Interpreter','latex')
xlabel('$t$ [s]','Interpreter','latex')
yline(0,'--');

ax(4) = nexttile; 
plot(t_roN,e_roN); box off;
ylabel('$\tilde x_2$','Interpreter','latex')
xlabel('$t$ [s]','Interpreter','latex')
yline(0,'--');

% saving plots
textwidth = 14;
golden_ratio = (1 + sqrt(5)) / 2;
textheight = textwidth / golden_ratio;
figsize = [textwidth, textheight];

% Get latex font in ticks
h = findall(gcf,'Type','axes'); % An array if you have subplots
set(h, 'TickLabelInterpreter', 'latex')

% Set size and no crop
set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', figsize);
set(gcf, 'PaperUnits', 'normalized', 'PaperPosition', [0, 0, 1, 1]);

% print -dpdf ../doc/figures/ex6_p2.pdf

%% functions 
function xdot = full_order(x,y,A,C)

    K = [2*5; 5^2]; % feedback gain 
    
    xdot = A * x + K*(y - C * x);
    
end

function wdot = reduced_order(w,y)

    L = 5; % feedback gain 
    
    wdot = -L * w - L^2 * y;
    
end

% clear; clc; 
% 
% %% system definition
% 
% x20 = 5;
% % states
% x = @(t) [x20 * t;...
%           x20 ];
% % output
% yo = @(t) x20 * t; 
% 
% % output with noise
% yoN = @(t) x20 * t + (rand() - 0.5)*5; % noise with zero average 
% 
% %% full-order observer
% 
% tspan = [0 5];
% xh10 = [0; 0];
% [t1, x1] = ode45(@(t,y) full_order(y, yo(t)), tspan, xh10);
% 
% %% full-order with noise
% 
% [t1N, x1N] = ode45(@(t,y) full_order(y, yoN(t)), tspan, xh10);
% 
% %% reduced-order observer
% 
% mass = [1 0 0;...v
%         0 0 0;...
%         0 0 0]; % mass matrix
% ode_opts = odeset('Mass',mass); 
% xvec0 = [1; -3; 5];
% [t2, x2] = ode15s(@(t,y) r_obs(y, yo(t)), tspan, xvec0, ode_opts);
% 
% %% reduced-order observer with noise
% 
% xvec0 = [1; -3; 4];
% [t2N, x2N] = ode15s(@(t,y) r_obs(y, yoN(t)), tspan, xvec0, ode_opts);
% 
% %% plots
% 
% figure(1); clf; 
% 
% % true states
% xv1 = @(t) x20 * t; 
% xv2 = @(t) x20 + 0*t; 
% xt1 = xv1(t1);
% xt2 = xv2(t1);
% xt1N = xv1(t1N);
% xt2N = xv2(t1N);
% xtt1 = xv1(t2);
% xtt2 = xv2(t2);
% xtt1N = xv1(t2N);
% xtt2N = xv2(t2N);
% 
% subplot(221)
% yyaxis left
% plot(t1,[xt1, x1(:,1)]); box off;
% ylabel('$x_1$','Interpreter','latex')
% yyaxis right
% plot(t1,[xt2, x1(:,2)]); box off;
% ylabel('$x_2$','Interpreter','latex')
% xlabel('$t$ [s]','Interpreter','latex')
% title('full','Interpreter','latex')
% 
% subplot(222)
% yyaxis left
% plot(t1N,[xt1N, x1N(:,1)]); box off;
% ylabel('$x_1$','Interpreter','latex')
% yyaxis right 
% plot(t1N,[xt2N, x1N(:,2)]); box off;
% ylabel('$x_2$','Interpreter','latex')
% xlabel('$t$ [s]','Interpreter','latex')
% title('full (noise)','Interpreter','latex')
% 
% subplot(223)
% yyaxis left
% plot(t2,[xtt1, x2(:,2)]); box off;
% ylabel('$x_1$','Interpreter','latex')
% yyaxis right
% plot(t2,[xtt2, x2(:,3)]); box off;
% ylabel('$x_2$','Interpreter','latex')
% xlabel('$t$ [s]','Interpreter','latex')
% title('reduced','Interpreter','latex')
% legend('$x$','$\hat x$')
% legend('Interpreter','latex');
% 
% subplot(224)
% yyaxis left
% plot(t2N,[xtt1N, x2N(:,2)]); box off;
% ylabel('$x_1$','Interpreter','latex')
% yyaxis right
% plot(t2N,[xtt2N, x2N(:,3)]); box off;
% ylabel('$x_2$','Interpreter','latex')
% xlabel('$t$ [s]','Interpreter','latex')
% title('reduced (noise)','Interpreter','latex')
% 
% % saving plots
% textwidth = 14.9;
% golden_ratio = (1 + sqrt(5)) / 2;
% textheight = textwidth / golden_ratio;
% figsize = [textwidth, textheight];
% 
% % Get latex font in ticks
% h = findall(gcf,'Type','axes'); % An array if you have subplots
% set(h, 'TickLabelInterpreter', 'latex')
% 
% % Set size and no crop
% set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', figsize);
% set(gcf, 'PaperUnits', 'normalized', 'PaperPosition', [0, 0, 1, 1]);
% 
% % print -dpdf ../doc/figures/ex6.pdf
% 
% %% function
% 
% function xdot = full_order(x,y) 
% 
% %%% initialization 
%     xdot = zeros(2,1);
%    
% %%% state equations
%     A = [0, 1;...
%          0, 0];
%     
%     K = 5;
%     
%     xdot = A * x + K * (y - x(1));
%     
% end
% 
% function xdot = r_obs(x,y) 
% 
% %%% initialization 
%     % x(1) = w
%     % x(2) = x1
%     % x(3) = x2
% 
%     xdot = zeros(3,1);
%     % xdot(1) = dot w
%     % xdot(2) = 0 = y - x(2)
%     % xdot(3) = 0 = x(3) - (x(1) + Ky) 
%     
% %%% state equations
%     
%     K = 5; %% same K as before
%     
%     xdot(1) = - K * x(1) - K * K * y;
%     xdot(2) = y - x(2);
%     xdot(3) = x(1) + K * y - x(3);
%     
% end
% 
