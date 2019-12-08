clc;
clear;

A = [-1 1 0 1 0 1; -1 -2 -1 0 0 1; 1 0 -2 -1 1 1; -1 1 -1 -2 0 0; -1 -1 1 1 -2 -1; 0 -1 0 0 -1 -3];
B = [0 -1 -1; 0 0 0; -1 1 1; -1 0 0; 0 0 1; -1 1 1];
C = [0 1 0 -1 -1 -1; 0 0 0 -1 0 0; 1 0 0 0 -1 0];
D = [0 0 0; 0 0 0; 0 0 0];
B1 = [B zeros(6,3)];
B2 = B;
C1 = [C; zeros(3,6)];
C2 = C;
D11 = [D zeros(3,3); zeros(3,3) zeros(3,3)];
D12 = [D; eye(3)];
D21 = [D eye(3)];
D22 = D;

An = sdpvar(6,6,'full');
Bn = sdpvar(6,3,'full');
Cn = sdpvar(3,6,'full');
Dn = sdpvar(3,3,'full');
X1 = sdpvar(6);
Y1 = sdpvar(6);
Z = sdpvar(6);
gam = sdpvar(1);

mat1 = [X1 eye(6) X1*A+Bn*C2 An X1*B1+Bn*D21 zeros(6,6);
        eye(6) Y1 A+B2*Dn*C2 A*Y1+B2*Cn B1+B2*Dn*D21 zeros(6,6);
       (X1*A+Bn*C2)' (A+B2*Dn*C2)' X1 eye(6) zeros(6,6) C1'+C2'*Dn'*D12';
       (An)' (A*Y1+B2*Cn)' eye(6) Y1 zeros(6,6) Y1*C1'+Cn'*D12';
       (X1*B1+Bn*D21)' (B1+B2*Dn*D21)' zeros(6,6) zeros(6,6) -gam*eye(6) D11'+D21'*Dn'*D12';
       zeros(6,6) zeros(6,6) (C1'+C2'*Dn'*D12')' (Y1*C1'+Cn'*D12')' (D11'+D21'*Dn'*D12')' -gam*eye(6)];
   

mat2 = [X1 eye(6);
        eye(6) Y1];
    
F = [mat1 >= 0; mat2>=0];
F = [F; X1>=0.0001; Y1>=0.0001; Z>=0.0001];
optimize(F, gam);
Hinf_norm = (value(gam))

[L,U] = lu(1-(value(X1)*value(Y1)));
X2 = L; Y2 = U';    

temp1 = [X2 X1*B2;
         zeros(3,6)   eye(3)];
temp2 = [An-X1*A*Y1 Bn;
         Cn Dn];
temp3 = [Y2' zeros(6,3);
         C2*Y1 eye(3)];
Mdk = inv(value(temp1))*value(temp2)*inv(value(temp3));
Adk =  Mdk(1:6,1:6);
Bdk =  Mdk(1:6,7:9);
Cdk =  Mdk(7:9,1:6);
Ddk =  Mdk(7:9,7:9);

%controller
Ddc = inv(eye(3)+Ddk*D22)*Ddk;
Cdc = (1-Ddc*D22)*Cdk;
Bdc = Bdk*(1-Ddc*D22);
Adc = Adk - Bdc*inv(eye(3)-D22*Ddc)*D22*Cdc;

display('the controller matrixes')
Adc
Bdc
Cdc
Ddc


