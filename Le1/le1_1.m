
clear; 

syms u1 u2 u3 x1 x2 x3 v1 v2 v3 lambda

L = sym('L', [3 3]);

A = [u1, 1, 0;
     0, u1, 1;   
     0, 0, u2];
 
C = [0, 0, 1];

X = [x1; ...
     x2; ...
     x3];

V = [v1; ...
     v2; ...
     v3];
 
O = [C; C*A; C*A^2];

rO = rref(O);

eqns_Ons = rO * V

NS = [lambda*eye(3) - A;C];

rLA = rref(NS);

eqns_rLA = rLA * V

% null(NS)
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
