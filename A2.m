clear;
clc;
A = - eye(3);
B = ones(3,1);
C = ones(1,3);
D = 0;
P = sdpvar(3);
gam = sdpvar(1);
mat = [A'*P*A-P A'*P*B C';
       (A'*P*B)' B'*P*B-gam D';
       C D -gam];
F = [mat <= 0; P>=0];
optimize(F, gam);
H_inf_norm = value(gam)