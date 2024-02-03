function xhat = EKF(model,init,data)

  M = length(data.y);
  P = init.P0;
  n = numel(init.x0);
  
  xhat = zeros(n,M); 
  xhat(:,1) = init.x0;
  
  for i=1:M-1
    % Measurement update
    Hi = model.dy(xhat(:, i));
    Ki = P*Hi'*inv(Hi*P*Hi' + model.R());
    P = P - Ki*Hi*P;
    xhat(:,i) = xhat(:, i) + Ki*(data.y(:, i) - model.y(xhat(:, i)));
    % Time update
    xhat(:,i+1) = model.x(xhat(:, i), i+1);
    Fi = model.dx(xhat(:, i));
    P = Fi*P*Fi' + model.Q();
  end
end