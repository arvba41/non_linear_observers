clear; 

omega = 1;
v0 = 10;

% system true states
x = @(t) 10*sin(omega*t);
v = @(t) 10*cos(omega*t);

A = [0 1; -omega^2 0];
C = [0 1];

y0 = [5; 0];
tspan = [0 10];
[t,y] = ode45(@(t,y) observ(y,v(t),C,A,omega),tspan, y0);

% y0 = [1; 0];
% tspan = [0 10];
% [t1,y1] = ode45(@(t,y) observ(y,v(t),C,A,omega),tspan, y0);

%% plots

figure(1); clf; 

tiledlayout(3,2)

nexttile([2 1])
plot(t,x(t),'LineWidth',2); ylabel('$x(t)$ [m]','Interpreter','latex') 
xlabel('$t$ [s]','Interpreter','latex'); hold on; 
plot(t,y(:,1),'--','LineWidth',2);
% plot(t1,y1(:,1),':','LineWidth',2);
box off;
ylim([-12 12]);

nexttile([2 1])
plot(t,v(t),'LineWidth',2); ylabel('$v(t)$ [m/s]','Interpreter','latex'); 
xlabel('$t$ [s]','Interpreter','latex'); hold on; 
plot(t,y(:,2),'--','LineWidth',2);
% plot(t1,y1(:,2),':','LineWidth',2);
box off;
ylim([-12 12]);

nexttile()
plot(t,x(t) - y(:,1),'LineWidth',2); ylabel('$e_x(t)$ [m]','Interpreter','latex'); 
xlabel('$t$ [s]','Interpreter','latex'); hold on; 
% plot(t1,y1(:,2),':','LineWidth',2);
box off;

nexttile()
plot(t,v(t) - y(:,2),'LineWidth',2); ylabel('$e_v(t)$ [m/s]','Interpreter','latex'); 
xlabel('$t$ [s]','Interpreter','latex'); hold on; 
% plot(t1,y1(:,2),':','LineWidth',2);
box off;

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

% print -dpdf ../doc/figures/ex417.pdf

%% functions

function xdot = observ(x,y,C,A,omega)

K = [0; 2*omega];

xdot = A*x + K*(y - C*x);

end