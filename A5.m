clear;
clc;
A = - eye(3);
B = ones(3,1);
C = ones(1,3);
D = 0;
P = sdpvar(3);
W = sdpvar(1,3,'full');
mat1 = [P A'*P+C'*W;
       (A'*P+C'*W)' P];

F = [mat1 >= 0; P>=0];
optimize(F, -P);
W_is = value(W)
L = inv(value(P))*value(W)'