function xhat = EKF(model,init,data)

  M = length(data.y);
  P = init.P0;
  n = numel(init.x0);
  
  xhat = zeros(n,M); xhat(:,1) = init.x0;
  
  for ii = 1:M-1
    % Measurement update
    Ht = model.hx(xhat(:,ii));
    Kt = P * Ht' / (Ht * P * Ht' + model.R);
    xhat(:,ii) = xhat(:,ii) + Kt * (data.y(:, ii) - model.h(xhat(:, ii)));
    P = P - Kt * Ht * P;
    
    % Time update
    xhat(:,ii+1) = model.f(xhat(:, ii), data.u(ii));
    Ft = model.fx(xhat(:, ii));
    P = Ft * P * Ft' + model.Q; %%%% assuming that G is 1
  end
end
