function J = J_logDetY(u_step)
%log-det(Y) cost for a single step of u

u=u_step;

global PPmat_GG xest_GG ePV_GG r1_GG r0_GG QQ_GG;
PP=PPmat_GG;
xcurr=xest_GG;
ePVe=ePV_GG;
r1=r1_GG;
r0=r0_GG;
QQ=QQ_GG;

[nX,~,~,~,~,dt]=get_problemSpecs;
%RR=@(d) r0*diag([ones(nX/2,1)])+r1*[abs(d(1)) 0;0 abs(d(2))];
RR2=@(d) r0*diag([ones(nX/2,1)])+r1*eye(2)*norm(d);

oneDiagNx=diag(ones(nX/2,1));
A_tr=[oneDiagNx dt*oneDiagNx
    zeros(nX/2,nX/2) oneDiagNx];
B_tr=[.5*dt^2*oneDiagNx
    dt*oneDiagNx];
H=[1 0 0 0
    0 1 0 0]; %position sensor

xPV2=A_tr*xcurr+B_tr*u;
ePV2=A_tr*ePVe;

PPap=A_tr'*PP*A_tr+QQ;
RR_est=RR2(xPV2(1:2)-ePV2(1:2));

PPinv=inv(PPap)+H'*inv(RR_est)*H;

Y=PPinv;

J=-log(det(Y));


end

