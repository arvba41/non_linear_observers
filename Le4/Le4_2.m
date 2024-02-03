clear; clc;
%% particle generation (initial)

N = 100; 
w = 1 / N * ones(1,N);
for ii = 1:N
    x(ii) = rand;
end

[fx,xk] = ksdensity(x);

% plotting the particle desity 
figure(10); clf; set(gcf,"WindowStyle","docked");
tiledlayout(1,2)

nexttile();
stem(x,w); box off;
ylabel('$w$','Interpreter','latex');
xlabel('$x$','Interpreter','latex');

nexttile();
plot(xk,fx); box off;
ylabel('$f(x)$','Interpreter','latex');
xlabel('$x$','Interpreter','latex');

h = findall(gcf,'Type','Axes');
set(h,'TickLabelInterpreter','latex');

%% mesurement update

sigma = 0.1;
mu = 0;
y = 0.7; % measuremenrt sample
for ii = 1:N
    wn(ii) = normpdf(y - x(ii)^2, mu, sigma);
    xn(ii) = rand; % new particles
end

wu = wn / sum(wn); % weight update (normalized)
[fxn,xnk] = ksdensity(xn);

% plotting the particles and a kernel density estimate of the distribution
% after measurement update
figure(11); clf; set(gcf,"WindowStyle","docked");
tiledlayout(1,2)

nexttile();
stem(xn,wu); box off;
ylabel('$w_n$','Interpreter','latex');
xlabel('$x_n$','Interpreter','latex');

nexttile();
plot(xnk,fxn); box off;
ylabel('$f(x_n)$','Interpreter','latex');
xlabel('$x_n$','Interpreter','latex');

h = findall(gcf,'Type','Axes');
set(h,'TickLabelInterpreter','latex');

%% number of effective particles
Neff = 1 / sum(wn.^2);

%% resampling
% the goal to re-sample use xn and then sample randomly around 
% a desired weight
idx = resample(wu); % desired weight index
xr = xn(idx);
wr = ones(1,N) * 1/N;
[fxr,xrk] = ksdensity(xr);

% plotting the particles and a kernel density estimate of the distribution
% after measurement update resampling
figure(12); clf; set(gcf,"WindowStyle","docked");
tiledlayout(1,2)

nexttile();
stem(xr,wr); box off;
ylabel('$w_r$','Interpreter','latex');
xlabel('$x_r$','Interpreter','latex');

nexttile();
plot(xrk,fxr); box off;
ylabel('$f(x_r)$','Interpreter','latex');
xlabel('$x_r$','Interpreter','latex');

h = findall(gcf,'Type','Axes');
set(h,'TickLabelInterpreter','latex');

%% Expextec value of x
x_hat = mean(xr); % mean value
[~,idx] = max(fxr); % take the maximum of the PDF
x_hat_new = xrk(idx);
