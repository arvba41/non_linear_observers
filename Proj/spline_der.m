clear
%% Generate some data and define spline interpolator
fun = @(x) 1 ./ (1 + exp(-x));  % sigmoid

x = linspace(-5, 5, 20);
y = fun(x);

figure(10)
plot(x, y, x, y, 'x')
box off

ocv = SplineInterpol(x, y);

%% Define spline interpolator

xi = linspace(-5, 5, 200);
figure(20)
plot(xi, ocv.f(xi))
hold on
plot(xi, ocv.df(xi))
hold off
box off



