clear;
clc;
A = [0.1 0.2 0.3;
     0.4 0.5 0.6;
     0.7 0.8 0.9];
Bd11 = ones(3,1);
Bd12 = ones(3,1);
Bd2 = ones(3,1);
Cd11 = ones(1,3);
Cd12 = ones(1,3);
Dd1111 = zeros(3,1);
Dd1112 = ones(3,1);
Dd1121 = ones(1,1);
Dd1122 = ones(1,1);
Dd121 = 1;
Dd122 = 1;
%Hinf norm is gam
gam = 1;
P = sdpvar(3);
Z = sdpvar(1);
Fd = sdpvar(1,3,'full');
musqr = sdpvar(1);
mat1 = [P A*P-Bd2*Fd Bd11;
       (A*P-Bd2*Fd)' P zeros(3,1);
       Bd11' zeros(1,3) eye(1)];
   
mat2 = [P A*P-Bd2*Fd Bd12 zeros(3,1);
       (A*P-Bd2*Fd)' P zeros(3,1) P*Cd12'-Fd'*Dd122';
       Bd12' zeros(1,3) gam*eye(1) Dd1122';
       zeros(1,3) (P*Cd12'-Fd'*Dd122')' Dd1122 gam*eye(1)];
   
   
mat3 = [Z Cd11*P+Dd121*Fd;
        (Cd11*P+Dd121*Fd)' P];
F = [mat1 >= 0; mat2>=0; mat3>=0; P>=0; Z>=0 trace(Z)<=musqr];
optimize(F, musqr);
Kd = value(Fd)*inv(value(P))
Hinf_norm = gam
H2_norm = sqrt(value(musqr))
