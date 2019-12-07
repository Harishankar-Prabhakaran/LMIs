clear;
clc;
A = [0.1 0.2 0.3;
     0.4 0.5 0.6;
     0.7 0.8 0.9];
Bd1 = ones(3,1);
Bd2 = ones(3,1);
Cd1 = ones(1,3);
Dd12 = 1;
P = sdpvar(3);
Z = sdpvar(1);
Fd = sdpvar(1,3);
musqr = sdpvar(1);
mat1 = [P A*P+Bd2*Fd Bd1;
       (A*P+Bd2*Fd)' P zeros(3,1);
       Bd1' zeros(1,3) eye(1)];
mat2 = [Z Cd1*P+Dd12*Fd;
        (Cd1*P+Dd12*Fd)' P];
F = [mat1 >= 0; mat2>=0; P>=0; Z>=0 trace(Z)<=musqr];
optimize(F, musqr);
Kd = value(Fd)*inv(value(P))
H2_norm = sqrt(value(musqr))