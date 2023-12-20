clear; 

syms e p1 p2 p3

A = [-2 , -1;...
     e  , -1];
 
 AT = [-2, e;...
       -1, -1];
 
C = [1, 1];

sigma =[p1, p2; p2, p3];

eq = AT * sigma + sigma * A + C' * C == zeros(size(A));

eqq = [eq(1), eq(2), eq(4)];

S = solve(eqq,[p1, p2, p3]);

Sigma = [S.p1, S.p2; S.p2, S.p3];

% lyap(A,C' * C)

% eq = A' * sigma + sigma * A + C' * C == zeros(size(A))
% 
% sigma_o = solve(eq ,sigma)

% 
% int_val1 = exp(A'*t)
% int_val2 = C'
% int_val3 = C 
% int_val4 = exp(A*t)
% 
% int_val = int_val1 * int_val2 * int_val3 * int_val4
% 
% int_val_simp = simplify(int_val)
% 
% % exp(-4*t)*(exp(t*(epsilon + 2)) + 1)*(exp(t*(conj(epsilon) + 2)) + 1)
% % 
% % 2*exp(-t)*(exp(t*conj(epsilon)) + exp(-2*t))
% % 
% % 2*exp(-3*t)*(exp(t*(epsilon + 2)) + 1)
% % 
% % 4*exp(-2*t)
% 
% sigma_o = int(int_val, t)
% 
% sigma_o_simp = simplify(sigma_o)
% 
% rank(sigma_o_simp)
% 
% %% 
% sigma_o_simp_simp = [1/4 - 2/(epsilon - 2) - 1/(2*epsilon), 2/3 - 2/(epsilon - 1);...
%     2/3 - 2 /(epsilon - 1), 2]
% 
% rank(sigma_o_simp_simp)
% 
% %% implementation usinfg simulation to check how the the integration looks like
% % tspan = [0 10]; 
% % sigma0 = [0; 0; 0; 0];
% % 
% % [t,y] = ode45()
% % 
% % 
% % %% functions
% % function sigma_o = gramian_fun(t,epsilon)
% %     integrand
% %     sigma_o(1) = exp(-4*t)*(exp(t*(epsilon + 2)) + 1)*(exp(t*(conj(epsilon) + 2)) + 1);
% %     sigma_o(2) = 2*exp(-t)*(exp(t*conj(epsilon)) + exp(-2*t));
% %     sigma_o(3) = 2*exp(-3*t)*(exp(t*(epsilon + 2)) + 1);
% %     sigma_o(4) = 4*exp(-2*t);
% %     sigma
% %     exp(t*(conj(epsilon) - 2))/(conj(epsilon) - 2) - exp(-4*t)/4 + exp(t*(epsilon - 2))/(epsilon - 2) + exp((t*(epsilon^2 + abs(epsilon)^2))/epsilon)/(2*real(epsilon))
% %     (2*exp(t*(epsilon - 1)))/(epsilon - 1) - (2*exp(-3*t))/3
% %     (2*exp(t*conj(epsilon))*exp(-t))/(conj(epsilon) - 1) - (2*exp(-3*t))/3
% %     -2*exp(-2*t)
%  
% % end
% 
