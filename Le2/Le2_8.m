clear; clc; 

%% system simultation 
u = 1; % input
fx = @(x,u) [x(1)^2 * x(2)^2 + u;...
             1 - x(2) * x(1)^3 - x(2) - x(2) / x(1) * u]; % system

tspan = [0 5];
y0 = [0.1; -1];
[t,y] = ode45(@(t,y) fx(y,u),tspan, y0);

%% observer
mesh = @(tt) interp1(t,y,tt,'spline');
A = [-1, 1; 0 0]; C = [1 0]; B = [0; 1];
fy = @(tt)[0; interp1(t,y(:,2),tt,'spline')];
K = [0.5; 1];

fz = @(t,z,u) A * z + B * u + fy(t) + K * (mesh(t)' - C * z);
[ts,Z] = ode45(@(t,y) fz(t,y,u),tspan, y0);
X_hat = [Z(:,2), Z(:,1)./Z(:,2)];

%% plots 
figure(1); clf; 
tiledlayout(2,2);

ax(1) = nexttile; 
plot(t,y(:,1),ts,X_hat(:,1)); box off;
ylabel('$x_2$','Interpreter','latex');
xlabel('$t$','Interpreter','latex');

ax(1) = nexttile; 
plot(t,y(:,2),ts,X_hat(:,2)); box off;
ylabel('$x_2$','Interpreter','latex');
xlabel('$t$','Interpreter','latex');

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
