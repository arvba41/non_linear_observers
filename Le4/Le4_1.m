clear; clc; 

%% section 1a
% done with sttistics and machine learning toolbox 

N = 500;
w = 1/N + zeros(1,N);

for ii = 1:N
    x(ii) = normrnd(0,1);
end

[fx,xk]=ksdensity(x);

% plots section 1s
figure(10); clf; set(gcf,'WindowStyle','docked');

tiledlayout(2,2)

nexttile();
plot(x); box off;
ylabel('$x$','Interpreter','latex');

nexttile();
stem(x,w); box off;
ylabel('$w$','Interpreter','latex');
xlabel('$x$','Interpreter','latex');

nexttile();
histogram(x); box off;
ylabel('$x$','Interpreter','latex');

nexttile();
plot(xk,fx); box off;
ylabel('$f(x)$','Interpreter','latex');
xlabel('$x$','Interpreter','latex');

h = findall(gcf,'Type','Axes');
set(h,'TickLabelInterpreter','latex');

% saving plots
textwidth = 14.9;
golden_ratio = (1 + sqrt(5)) * 0.5;
textheight = textwidth / golden_ratio;
figsize = [textwidth, textheight];

% Set size and no crop
set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', figsize);
set(gcf, 'PaperUnits', 'normalized', 'PaperPosition', [0, 0, 1, 1]);

% print -dpdf ../doc/figures/ex4_a.pdf

%% section 1b

y = sin(2 * pi * x.^2) .* abs(x).^(1/2);
[fy,yk]=ksdensity(y);

figure(11); clf; set(gcf,'WindowStyle','docked');

tiledlayout(1,2)

nexttile();
stem(y,w); box off;
ylabel('$w$','Interpreter','latex');
xlabel('$y$','Interpreter','latex');

nexttile();
plot(yk,fy); box off;
ylabel('$f(y)$','Interpreter','latex');
xlabel('$y$','Interpreter','latex');

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

% print -dpdf ../doc/figures/ex4_b.pdf