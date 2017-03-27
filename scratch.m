clear;clc;

d=sym('d', 'real');
f=sym('f','real');
v=sym('v','real');
b=sym('b','real'); 
b=d^2/2;
syms A B Qf r R Pn Pnm1 Pnm2 Pnm3 Pnm4 Pnm5
A=[1 d;0 1];
B=[b;d];
Qf=[f 0;0 0];
Q=[0 0;0 0];
R=r;

Pn=Qf;
Pnm1=simplify(Q+A'*Pn*A-A'*Pn*B*inv(R+B'*Pn*B)*B'*Pn*A);
Pnm2=simplify(Q+A'*Pnm1*A-A'*Pnm1*B*inv(R+B'*Pnm1*B)*B'*Pnm1*A);
Pnm3=simplify(Q+A'*Pnm2*A-A'*Pnm2*B*inv(R+B'*Pnm2*B)*B'*Pnm2*A);
Pnm4=simplify(Q+A'*Pnm3*A-A'*Pnm3*B*inv(R+B'*Pnm3*B)*B'*Pnm3*A);
Pnm5=simplify(Q+A'*Pnm4*A-A'*Pnm4*B*inv(R+B'*Pnm4*B)*B'*Pnm4*A);
Pnm6=simplify(Q+A'*Pnm5*A-A'*Pnm5*B*inv(R+B'*Pnm5*B)*B'*Pnm5*A);


simplify(-inv(R+B'*Pnm1*B)*B'*Pnm1*A)
simplify(-inv(R+B'*Pnm2*B)*B'*Pnm2*A)
simplify(-inv(R+B'*Pnm3*B)*B'*Pnm3*A)
simplify(-inv(R+B'*Pnm4*B)*B'*Pnm4*A)
simplify(-inv(R+B'*Pnm5*B)*B'*Pnm5*A)
simplify(-inv(R+B'*Pnm6*B)*B'*Pnm6*A)








