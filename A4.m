clear;
clc;
A = [0.1 0.2 0.3;
     0.4 0.5 0.6;
     0.7 0.8 0.9];
B = ones(3,1);
C = ones(1,3);
D = 0;
P = sdpvar(3);
W = sdpvar(1,3,'full');
mat1 = [P A*P+B*W;
       (A*P+B*W)' P];

F = [mat1 >= 0; P>=0];
optimize(F, -P);
W_is = value(W)
K_d = value(W)*inv(value(P))
