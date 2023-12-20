clear;

syms a k

A = [-(3+a) 1 0;...
     -(3*a+k) 0 1;...
     -k 0 0];
 
C = [1 0 0];

O = [C; C*A; C*A^2]