clear;clc;
%Deterministic DG main


n = 2; %time steps to rendezvous
nX=6;  %num x-states, here pos+vel
nU=3;
xPV0=[10;0;0 ; 0;0;0];
dtset=1; %update period, seconds
%set blanks for fdyn_eva or fdyn_pur for unknown pursuer/evader dynamics
set_problemSpecs(nX,nU,'fdyn_eva','fdyn_pur',xPV0,dtset);

%nX/2 since you only care that xFeva=xFpur, don't care about velocity
Arendezvous=[zeros(3,(nX+nU)*(n-1)) diag(ones(nX/2,1)) zeros(nX/2,nX/2) zeros(nU,nU) ];

%generation of xu0, the initial estimate for full x- and u- state history
%plan: u=0 until last time step
xu0=[];
v0=xPV0(1+nX/2:end);
xprev=xPV0(1:nX/2);
for i=1:n-1
    xnew=xprev+dtset*v0;
    vnew=v0;
    xu0=[xu0;xnew;vnew;zeros(nU,1)];
    xprev=xnew;
    vprev=vnew;
end
evaFinalPos=feval('fdyn_eva',n);
ufinal=1/dtset*(evaFinalPos-(xprev+dtset*v0));
xu0=[xu0;evaFinalPos;v0+ufinal*dtset;ufinal];
[xuhist,J]=fmincon('J_pursuer',xu0,[],[],Arendezvous,evaFinalPos,[],[],'dynamicsConstraint');

xuhist
J














