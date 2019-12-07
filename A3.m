clear;
clc;
A = - eye(3);
B = ones(3,1);
C = ones(1,3);
D = 0;
P = sdpvar(3);
Z = sdpvar(1);
musqr = sdpvar(1);
mat1 = [P A*P B;
       (A*P)' P zeros(3,1);
       B' zeros(1,3) eye(1)];
mat2 = [Z C*P;
        (C*P)' P];
F = [mat1 >= 0; mat2>=0; P>=0; trace(Z)<=musqr];
optimize(F, musqr);
H2_norm = sqrt(value(musqr))