clear;
clc;
A = - eye(5); %try this I to get +ve definite P
%B = [0 0; 0 0; 0 0; 20 0; 0 20]; 
%Q = eye(5);
P = sdpvar(5,5);

F = [(A'*P*A-P)<=0; P >= 0.0001];
optimize(F, -P);

d = eig(value(P));
fprintf('\n');
if d > 0
    disp(['THE POSITIVE DEFINITE MATRIX P']);
    value(P)
else
    disp(['No positive definite P exists']);
    disp(['System is unstable in open-loop']);
    disp(['Eigen values of A must be non negative']);
    Eigen_A = eig(A)
end