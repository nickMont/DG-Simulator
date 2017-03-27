clear;clc;

dt=sym('dt','real');
syms re rp R L A B Q q1p q1e q2p q2e u1 u2 x1 x2 v1 v2 x u up1 up2 ue1 ue2 ue3 up3
B=[eye(2)*dt^2/2; eye(2)*dt];
A=[eye(2) eye(2)*dt; zeros(2,2) eye(2)];
Rp=rp*eye(2);
Re=re*eye(2);
Qp=diag([q1p q1p q2p q2p]);
Qe=diag([-q1e -q1e q2e q2e]);
u=[u1;u2]; x=[x1;x2;v1;v2];

%lvl1
up1=inv(Rp+B'*Qp*B)*(-B'*Qp*B*u+B'*Qp*A*x);
ue1=simplify(inv(Re+B'*Qe*B)*( B'*Qe*B*up1+B'*Qe*A*x));
%lvl2
up2=inv(Rp+B'*Qp*B)*(-B'*Qp*B*ue1+B'*Qp*A*x);
ue2=simplify(inv(Re+B'*Qe*B)*( B'*Qe*B*up2+B'*Qe*A*x));
%lvl3
up3=inv(Rp+B'*Qp*B)*(-B'*Qp*B*ue2+B'*Qp*A*x);
ue3=simplify(inv(Re+B'*Qe*B)*( B'*Qe*B*up3+B'*Qe*A*x));

evalvars={'dt','q1p','q1e','q2p','q2e','rp','re'};
evalsubs={.25 , 1 , 2 , 1 , 0 , 1 , 2 };
diff_12=vpa(simplify(subs(ue2(1)-ue1(1),evalvars,evalsubs)),5)
diff_23=vpa(simplify(subs(ue3(1)-ue2(1),evalvars,evalsubs)),5)
















