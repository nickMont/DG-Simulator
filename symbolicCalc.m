clear;clc;

dt=sym('dt','real');
syms R L A B Q u1 u2 x1 x2 v1 v2 x u up1 up2 ue1 ue2 ue3 up3
q1e=sym('q1e','real');
q1p=sym('q1p','real');
q2e=sym('q2e','real');
q2p=sym('q2p','real');
re=sym('re','real');
rp=sym('rp','real');

B=[eye(2)*dt^2/2; eye(2)*dt];
A=[eye(2) eye(2)*dt; zeros(2,2) eye(2)];
Rp=rp*eye(2);
Re=re*eye(2);
Qp=diag([q1p q1p q2p q2p]);
Qe=diag([-q1e -q1e q2e q2e]);
u=[u1;u2]; x=[x1;x2;v1;v2];
syms e ekp1 Ce Cp
e1=sym('e1','real');
e2=sym('e2','real');
e3=sym('e3','real');
e4=sym('e4','real');
umax=sym('umax','real');
normE=sym('normE','real');
e=[e1;e2;e3;e4];

inv(Rp+B'*Qp*B)
B'*Qp*A*e
ue1=-inv(Re+B'*Qe*B)*(B'*Qe*A*e - B'*Qe*B*umax*[e1;e2]/normE);
up1=inv(Rp+B'*Qp*B)*(B'*Qp*A*e + B'*Qp*B*umax*[e1;e2]/normE);

Ce=simplify((A*e+B*ue1-B*up1)'*Qe*(A*e+B*ue1-B*up1)+ue1'*Re*ue1)
Cp=simplify((A*e+B*ue1-B*up1)'*Qp*(A*e+B*ue1-B*up1)+up1'*Rp*up1)











