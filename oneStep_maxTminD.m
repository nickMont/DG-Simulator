function [cineq,ceq] = oneStep_maxTminD(u)

global umax_GG xCurr_GG distMax_GG dt_GG u_GG nX_GG cd_GG
umax=umax_GG;
dt=dt_GG;
distMax=distMax_GG;
xCurr=xCurr_GG;
uPur=u_GG;
nX=nX_GG;
cd=cd_GG;

A_tr=[eye(nX/2) (1*dt-cd*dt^2/2)*eye(nX/2); zeros(nX/2,nX/2) (1-cd*dt)*eye(nX/2)];
B_tr=[dt^2/2*eye(nX/2); dt*eye(nX/2)];

x2=A_tr*xCurr+B_tr*u-B_tr*uPur;

ceq=[];
cineqT=norm(u)-umax;
cineqD=distMax-norm(x2(1:nX/2));
cineq=[cineqT
    cineqD];


end

