clear;
clc;
A = [0.1 0.2 0.3;
     0.4 0.5 0.6;
     0.7 0.8 0.9];
Bd1 = ones(3,1);
Bd2 = ones(3,1);
Cd1 = ones(1,3);
Dd12 = 1;
Dd11 = 0;
P = sdpvar(3);
Z = sdpvar(1);
Fd = sdpvar(1,3);
gam = sdpvar(1);
mat1 = [P A*P+Bd2*Fd Bd1 zeros(3,1);
       (A*P+Bd2*Fd)' P zeros(3,1) P*Cd1'-Fd'*Dd12';
       Bd1' zeros(1,3) gam*eye(1) Dd11';
       zeros(1,3) (P*Cd1'-Fd'*Dd12')' Dd11 gam*eye(1)];
F = [mat1 >= 0; P>=0];
optimize(F, gam);
Kd = value(Fd)*inv(value(P))
H2_norm = (value(gam))