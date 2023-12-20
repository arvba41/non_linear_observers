clear
close all

%% Simulate robot movement
Ts = 0.1;
V = 3;

model.f = @(x,u,w) [x(1) + Ts*V*cos(x(3));x(2) + Ts*V*sin(x(3));x(3) + Ts*u];
model.h = @(x,u) [sqrt(x(1,:).^2 + x(2,:).^2);atan(x(2,:)./x(1,:))];

x0 = [10;10;0];
Tfinal = 20;

data.t = (0:Ts:Tfinal);
N = numel(data.t);
data.x = zeros(3,N); data.x(:,1) = x0;

data.u = zeros(1,N);
data.u(data.t>=5 & data.t < 10) = -1;
data.u(data.t>=10 & data.t < 12) = 0;
data.u(data.t>=12 & data.t < 15) = 1;
data.u(data.t>=15 & data.t < 17) = -1;
for t=1:N-1
  data.x(:,t+1) = model.f(data.x(:,t),data.u(t),0);
end
data.y0 = model.h(data.x, data.u);
data.y  = data.y0 + diag([0.6,0.03])*randn(2,N); 

%% Plot mobile robot path
figure(1)
plot( data.x(1,:)', data.x(2,:)' )
hold on
plot( data.x(1,1)', data.x(2,1)', 'bo', ...
  data.x(1,end)', data.x(2,end)', 'rx');
hold off
axis([0 max(data.x(1,:))*1.1 0 max(data.x(2,:))*1.1]);
axis square
xlabel('x');
ylabel('y');

%% Plot measurements
figure(2)
subplot( 211 )
plot( data.t', [data.y(1,:)', data.y0(1,:)'])
xlabel('t [s]');
ylabel('[m]')
title('Range measurement');
box off

subplot( 212 )
plot( data.t', 180/pi*[data.y(2,:)', data.y0(2,:)'])
xlabel('t [s]');
ylabel('[deg]');
title('Angle measurement');
box off

%% EKF estimations
% model
modelEKF.f = @(x,u) [x(1) + Ts * V * cos(x(3));...
                     x(2) + Ts * V * sin(x(3));...
                     x(3) + Ts * u];
% measurement 
modelEKF.h = @(x,u) [sqrt(x(1)^2 + x(2)^2);...
                     atan(x(2) / x(1))];
% linearization f
modelEKF.fx = @(x,u) [1 0 -Ts * V * sin(x(3));...
                      0 1 Ts * V * cos(x(3));...
                      0 0 1];
% linearization h                  
modelEKF.hx = @(x,u) [x(1) / sqrt(x(1)^2 + x(2)^2),                 x(2) / sqrt(x(1)^2 + x(2)^2),   0;...
                      x(2) / (x(1)^2 * (x(2)^2 / x(1)^2 + 1)),  1 / (x(2) * (x(2)^2 / x(1)^2 + 1)), 0];

% turing parameters 
% modelEKF.Q = diag([0.8 0.8 0.1]);
% modelEKF.R = diag([0.8,0.5]); % diag([0.6,0.03])
modelEKF.Q = diag([0 0 0.1]);
modelEKF.R = diag([100 5]) * 10;

init.x0 = [10; 5; 0]; % Initial stafte estimate 
init.P0 = diag([1000, 1000, 0]); % Initial covariance
tic; xhatEKF = EKF(modelEKF,init,data); toc

%% performing an optimization

params.Ts = Ts;
params.V = V;

x0 = [0 0 0.1 ...
    3 10 ...
    10, 10, 1];

% x0 = x; % homotopy :)
options = optimoptions('lsqnonlin','Display','iter','FunctionTolerance',1e-12);
[x,resnorm,residual,exitflag,output] = lsqnonlin(@(x) EKF_res(x,data,params), x0, zeros(size(x0)), [], options);

%% running EKF with new values
modelEKF.Q = diag(x(1:3));
modelEKF.R = diag(x(4:5));
init.P0 = diag([x(6) x(7) x(8)]); % Initial covariance

tic; xhatEKF_opt = EKF(modelEKF,init,data); toc

%% performing an optimization for all 

params.Ts = Ts;
params.V = V;

% x00 = [0 0 0.1, 0 0 0, 0 0 0.1 ...
%     1 0, 0 10 ...
%     1 0 0,0 1 0,0 0 1];

x00 = x_f; % homotopy :)
options = optimoptions('lsqnonlin','Display','iter','FunctionTolerance',1e-12);
[x_f,resnorm_f,residual_f,exitflag_f,output_f] = lsqnonlin(@(x) EKF_res_all(x,data,params), x00, zeros(size(x00)), [], options);

%% running EKF with new values
modelEKF.Q = [x_f(1:3); x_f(4:6); x_f(7:9)];
modelEKF.R = [x_f(10:11); x_f(12:13)];
init.P0 = [x_f(14:16); x_f(17:19); x_f(20:22)]; % Initial covariance

tic; xhatEKF_opt_f = EKF(modelEKF,init,data); toc

%% Plot path 
figure(3); clf; 

subplot(221);
plot( data.x(1,:)', data.x(2,:)', 'LineWidth', 1)
hold on
plot( data.x(1,1)', data.x(2,1)', 'bo', ...
  data.x(1,end)', data.x(2,end)', 'rx');
plot( (data.y(1,:).*cos(data.y(2,:)))', (data.y(1,:).*sin(data.y(2,:)))', 'k')
hold off
axis([0 max(data.x(1,:))*1.1 0 max(data.x(2,:))*1.1]);
axis square
xlabel('x');
ylabel('y');
title('Direct computation from measurements');


subplot(222);
plot( data.x(1,:)', data.x(2,:)', 'LineWidth', 1)
hold on
plot(xhatEKF(1,:)', xhatEKF(2,:)');
plot( data.x(1,1)', data.x(2,1)', 'bo', ...
  data.x(1,end)', data.x(2,end)', 'rx');
hold off
axis([0 max(data.x(1,:))*1.1 0 max(data.x(2,:))*1.1]);
axis square
xlabel('x');
ylabel('y');
title('EKF estimation')


subplot(223);
plot( data.x(1,:)', data.x(2,:)', 'LineWidth', 1)
hold on
plot(xhatEKF_opt(1,:)', xhatEKF_opt(2,:)');
plot( data.x(1,1)', data.x(2,1)', 'bo', ...
  data.x(1,end)', data.x(2,end)', 'rx');
hold off
axis([0 max(data.x(1,:))*1.1 0 max(data.x(2,:))*1.1]);
axis square
xlabel('x');
ylabel('y');
title('EKF estimation (opt)')

subplot(224);
plot( data.x(1,:)', data.x(2,:)', 'LineWidth', 1)
hold on
plot(xhatEKF_opt_f(1,:)', xhatEKF_opt_f(2,:)');
plot( data.x(1,1)', data.x(2,1)', 'bo', ...
  data.x(1,end)', data.x(2,end)', 'rx');
hold off
axis([0 max(data.x(1,:))*1.1 0 max(data.x(2,:))*1.1]);
axis square
xlabel('x');
ylabel('y');
title('EKF estimation (opt-full)')


%%
figure(20)

tiledlayout(2,3,'TileIndexing','columnmajor');

nexttile();
plot( data.t', [xhatEKF_opt(1,:)', xhatEKF(1,:)', data.x(1,:)'])
xlabel('t [s]');
ylabel('[m]')
title('x position estimation');
box off

nexttile();
plot( data.t', [xhatEKF_opt(1,:)'-data.x(1,:)', xhatEKF(1,:)'-data.x(1,:)'])
xlabel('t [s]');
ylabel('[m]')
title('Estimation error');
box off

nexttile();
plot( data.t', [xhatEKF_opt(2,:)', xhatEKF(2,:)', data.x(2,:)']); 
xlabel('t [s]');
ylabel('[m]')
title('y position estimation');
box off

nexttile();
plot( data.t', [xhatEKF_opt(2,:)'-data.x(2,:)', xhatEKF(2,:)'-data.x(2,:)'])
xlabel('t [s]');
ylabel('[m]');
title('Estimation error');
box off

nexttile();
plot( data.t', 180/pi*[xhatEKF_opt(3,:)', xhatEKF(3,:)', data.x(3,:)'])
xlabel('t [s]');
ylabel('[deg]')
title('Heading estimation');
box off

nexttile();
plot( data.t', 180/pi*([xhatEKF_opt(3,:)'-data.x(3,:)', xhatEKF(3,:)'-data.x(3,:)']))
xlabel('t [s]');
ylabel('[deg]');
title('Estimation error');
box off

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

newcolors = ["#1f77b4" "#ff7f0e" "#2ca02c" "#d62728"];
colororder(newcolors)

%% functions
%%%% optimization part 1
function Xhat = modelEKF_opti(x,data,params)

Ts = params.Ts;
V = params.V;

% model
modelEKF.f = @(x,u) [x(1) + Ts * V * cos(x(3));...
                     x(2) + Ts * V * sin(x(3));...
                     x(3) + Ts * u];
% measurement 
modelEKF.h = @(x,u) [sqrt(x(1)^2 + x(2)^2);...
                     atan(x(2) / x(1))];
% linearization f
modelEKF.fx = @(x,u) [1 0 -Ts * V * sin(x(3));...
                      0 1 Ts * V * cos(x(3));...
                      0 0 1];
% linearization h                  
modelEKF.hx = @(x,u) [x(1) / sqrt(x(1)^2 + x(2)^2),                 x(2) / sqrt(x(1)^2 + x(2)^2),   0;...
                      x(2) / (x(1)^2 * (x(2)^2 / x(1)^2 + 1)),  1 / (x(2) * (x(2)^2 / x(1)^2 + 1)), 0];

% turing parameters 
% modelEKF.Q = diag([0.8 0.8 0.1]);
% modelEKF.R = diag([0.8,0.5]); % diag([0.6,0.03])
modelEKF.Q = diag([x(1) x(2) x(3)]);
modelEKF.R = diag([x(4) x(5)]);

init.x0 = [10; 5; 0]; % Initial stafte estimate 
init.P0 = diag([x(6) x(7) x(8)]); % Initial covariance
Xhat = EKF(modelEKF,init,data);
end 

function residual = EKF_res(x,data,params)
Xhat = sum(modelEKF_opti(x,data,params),1);
X = sum(data.x,1);
residual = X - Xhat;
end 


%%%%% optimizaion with inter cross correlations
%%%% optimization part 1
function Xhat = modelEKF_opti_all(x,data,params)

Ts = params.Ts;
V = params.V;

% model
modelEKF.f = @(x,u) [x(1) + Ts * V * cos(x(3));...
                     x(2) + Ts * V * sin(x(3));...
                     x(3) + Ts * u];
% measurement 
modelEKF.h = @(x,u) [sqrt(x(1)^2 + x(2)^2);...
                     atan(x(2) / x(1))];
% linearization f
modelEKF.fx = @(x,u) [1 0 -Ts * V * sin(x(3));...
                      0 1 Ts * V * cos(x(3));...
                      0 0 1];
% linearization h                  
modelEKF.hx = @(x,u) [x(1) / sqrt(x(1)^2 + x(2)^2),                 x(2) / sqrt(x(1)^2 + x(2)^2),   0;...
                      x(2) / (x(1)^2 * (x(2)^2 / x(1)^2 + 1)),  1 / (x(2) * (x(2)^2 / x(1)^2 + 1)), 0];

% turing parameters 
% modelEKF.Q = diag([0.8 0.8 0.1]);
% modelEKF.R = diag([0.8,0.5]); % diag([0.6,0.03])
modelEKF.Q = [x(1:3); x(4:6); x(7:9)];
modelEKF.R = [x(10) x(11); x(12) x(13)];

init.x0 = [10; 5; 0]; % Initial stafte estimate 
init.P0 = [x(14:16);x(17:19);x(20:22)]; % Initial covariance
Xhat = EKF(modelEKF,init,data);
end 

function residual = EKF_res_all(x,data,params)
Xhat = sum(modelEKF_opti_all(x,data,params),1);
X = sum(data.x,1);
residual = X - Xhat;
end 