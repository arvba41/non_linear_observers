%% Battery State-of-Charge Estimation
% This example shows how to estimate the battery state-of-charge (SOC) by 
% using a Kalman filter. The initial SOC of the battery is equal to 0.5. 
% The estimator uses an initial condition for the SOC equal to 0.8. The 
% battery keeps charging and discharging for 6 hours. The extended Kalman 
% filter estimator converges to the real value of the SOC in less than 10 
% minutes and then follows the real SOC value. To use a different Kalman 
% filter implementation, in the SOC Estimator (Kalman Filter) block, set 
% the Filter type parameter to the desired value.

% Copyright 2022 The MathWorks, Inc.

%% Model

open_system('BatterySOCEstimation')

set_param(find_system('BatterySOCEstimation','FindAll', 'on','type','annotation','Tag','ModelFeatures'),'Interpreter','off')

%% Simulation Results
%
% The plot below shows the real and estimated battery state-of-charge.
%


BatterySOCEstimationPlotSOC;

%% Results from Real-Time Simulation
%%
%
% This example has been tested on a Speedgoat Performance real-time target 
% machine with an Intel(R) 3.5 GHz i7 multi-core CPU. This model can run 
% in real time with a step size of 50 microseconds.