clear 
Ts = 0.1;
V = 3;
u = -1;

tvec = 0:Ts:20; % time vector
N = numel(tvec); % number of sample points

% input control sequence 
u = zeros(1,N);
u(tvec >= 5  & tvec < 10) = -1;
u(tvec >= 10 & tvec < 12) = 0;
u(tvec >= 12 & tvec < 15) = 1;
u(tvec >= 15 & tvec < 17) = -1;
% u = 1;

n=3;      % number of states
q = 0.01;   % std of process 
r = 1;  % std of measurement

Qs = diag([0 0 0]); % std matrix of process
Rs = diag([0.6, 0.03].^2); % std of measurement  

f = @(x,i) [x(1) + Ts*V*cos(x(3));...
        x(2) + Ts*V*sin(x(3));...
        x(3) + Ts*u(i)];  % nonlinear state equations

h = @(x,i) [sqrt(x(1).^2 + x(2).^2);...
        atan(x(2)./x(1))]; % measurement equation

s = [10; 10; 0];          % initial state
x = [10; 5 ; 0];          % initial state with noise
S = diag([0.01 0.01 0.01]);% initial square root of state covraiance

xV = zeros(n,N);          % estmate        % allocate memory
sV = zeros(n,N);          % actual
zV = zeros(2,N);

for kk = 1 : N
  z = h(s,kk) + [0.6 ; 0.03]*randn; % measurments
  sV(:,kk)= s;                      % save actual state
  zV(:,kk)  = z;                    % save measurment
  [x, S] = ukf(f,x,S,h,z,Qs,Rs,kk); % ekf 
  xV(:,kk) = x;                     % save estimate
  s = f(s,kk);                      % update process 
end
for kk=1:3                                 % plot results
  subplot(3,1,kk)
  plot(1:N, sV(kk,:), '-', 1:N, xV(kk,:), '--')
end


figure(3); clf; 
plot(sV(1,:),sV(2,:)); hold on; 
plot(xV(1,:),xV(2,:))