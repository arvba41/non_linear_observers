clear;

syms L mu c1

A = [L 1 0;...
     0 L 1;...
     0 0 L];

C = sym('C',[1 3]);

% C = [0 c1 c2];

O = [C; C*A; C*A^2]

rank(O)


Aa = [L 1 0;...
     0 L 0;...
     0 0 mu];

Cc = sym('C',[1 3]);

Oc = [Cc; Cc*Aa; Cc*Aa^2]