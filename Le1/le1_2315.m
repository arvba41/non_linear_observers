clear;

%% transfer function
syms a k 

% s = tf('s'); 

% H = (s + 1) / (s * (s + 3))
% G = k / (s + a)


A = [-(3+a) 1 0;...
     -(3*a+k) 0 1;...
     -k 0 0];
 
C = [1 0 0];

O = [C; C*A; C*A^2];

subs(O,[a k],[0 0])
rank(subs(O,[a k],[0 0]))