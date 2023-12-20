
t = linspace(0,20,2000);
y = sawtooth(2* pi * 0.25 * t);
% y(y <= -2/1000) = 0;
y1 = y;
y(t > 4 + 10/1000) = NaN;
y1(t <= 4 + 10/1000) = NaN;

clf;
subplot(211)
p1 = plot(t,y); hold on;
p2 = plot(t,y1,'--');
p2.Color = p1.Color;
ylabel('$h(x)$','Interpreter','latex');
xlabel('$x$','Interpreter','latex');
box off; xlim([0 14]); ylim([0 1])

subplot(212)
p1 = plot(t,t); hold on;
p2 = plot(t,t+2,'--');
p3 = plot(t,t+4,'--');
p2.Color = p1.Color;
p3.Color = p1.Color;
ylabel('$x$','Interpreter','latex');
xlabel('$t$','Interpreter','latex');
box off; xlim([0 4])
yline([2 4],'-',{' ','exitation area'},'Interpreter','latex')

% saving plots
textwidth = 13;
golden_ratio = (1 + sqrt(5)) / 2;
textheight = textwidth / golden_ratio;
figsize = [textwidth, textheight];

% Get latex font in ticks
h = findall(gcf,'Type','axes'); % An array if you have subplots
set(h, 'TickLabelInterpreter', 'latex')

% Set size and no crop
set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', figsize);
set(gcf, 'PaperUnits', 'normalized', 'PaperPosition', [0, 0, 1, 1]);

print -dpdf ../doc/figures/ex2_3.pdf