function xhat = EKF_sr(model,init,data)

  M = length(data.y);
  P = init.P0;
  n = numel(init.x0);
  
  xhat = zeros(n,M); xhat(:,1) = init.x0;
  
  for ii = 1:M-1
    % Measurement update
    Ht = model.hx(xhat(:,ii));
    
    % creating A matrix (cell with matrices)
    A = cell(2, 2); 
    A{1,1} = sqrtm(model.R);
    A{1,2} = Ht * sqrtm(P);
    A{2,1} = zeros(3,2);
    A{2,2} = sqrtm(P);
    
    % QR Factorization
    A = cell2mat(A);
    [Q,~] = qr(A'); 
    
    A_theta = A * Q;
    X = A_theta(1:2, 1:2);
    Y = A_theta(3:5, 1:2);
    Z = A_theta(3:5, 3:5);
    
    P = Z^2;
    Kt = Y / X; % equivalent to Y * inv(X)

    xhat(:,ii) = xhat(:,ii) + Kt * (data.y(:, ii) - model.h(xhat(:, ii)));

    % Time update
    xhat(:,ii+1) = model.f(xhat(:, ii), data.u(ii));
    Ft = model.fx(xhat(:, ii));
    P = Ft * P * Ft' + model.Q; %%%% assuming that G is 1
  end
end
