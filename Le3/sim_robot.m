%% Simulate robot movement
function [model, data, params] = sim_robot(params,plots)
    
    Ts = params.Ts; % sample time 
    V = params.V; % constant velocity
    
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

    %% plotting data 
    if plots
        % Plot mobile robot path
        figure(1);
        plot( data.x(1,:)', data.x(2,:)' )
        hold on
        plot( data.x(1,1)', data.x(2,1)', 'bo', ...
          data.x(1,end)', data.x(2,end)', 'rx');
        hold off
        axis([0 max(data.x(1,:))*1.1 0 max(data.x(2,:))*1.1]);
        axis square
        xlabel('x');
        ylabel('y');
        
        % Plot measurements
        figure(2);
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
    end

    params.Ts = Ts;
    params.V = V;
end