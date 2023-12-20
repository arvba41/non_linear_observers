clear; clc; 

tvec = linspace(-0.02,0.02,100);
u = sin(2*pi*50*tvec);

x = cos(2*pi*50*tvec);

h = x.^2;

hx = 2 * x;

figure(1); clf; 
tiledlayout(2,1)

nexttile;
plot(tvec,x); box off;
ylabel('$x$','Interpreter','latex');
xlabel('$t$','Interpreter','latex');

nexttile;
plot(x,[h; hx]); box off;
ylabel('$y$','Interpreter','latex');
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
