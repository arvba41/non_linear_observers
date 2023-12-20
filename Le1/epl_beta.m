clear; clc; 

%% system oberver with known beta

beta  = 1;
tspan = [0 10]; 
u = @(t) heaviside(t - 2);
x0 = 0;
[t1, x1] = ode45(@(t,y) lpf(y, beta, u(t)),tspan,x0);
y1 = x1;
%% plots

figure(1)
subplot(211);
plot(t1,y1)
ylabel('y'); xlabel('t');

subplot(212);
plot(t1,u(t1))
ylabel('u'); xlabel('t');

%% observer with known beta
x10 = 0;

yy = @(t) interp1(t1,y1,t);
x1 = @(t) interp1(t1,x1,t);

[t2, x2] = ode45(@(t,y) obs1(y, x1(t), beta, yy(t)), tspan, x0);
y2 = x2;

%% plots

figure(2)
plot(t1,y1,t2,y2)

%% functions
function xdot = lpf(x, beta, u)    
    xdot = - beta * x + u;
end

function xdot = obs1(x, xr, beta, y) 
    A = - beta;
    
    K = 0; 
    
    xdot =  A * xr + K * (y - x);
end