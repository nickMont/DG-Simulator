function [cineq, ceq] = dynamicsConstraint(xu_hist)


[nX,nU,~,~,xPV0,dt]=get_problemSpecs;
sizeT=length(xu_hist);
n=sizeT/(nX+nU);

g=9.81;

xhist=[];  %nX x n
uhist=[];  %nU x n

xutemp=reshape(xu_hist,[nX+nU,n]);
xhist=xutemp(1:nX,:);
uhist=xutemp(nX+1:end,:);

xPV=xPV0;

oneDiagNx=diag(ones(nX/2,1));
A_tr=[oneDiagNx dt*oneDiagNx
      zeros(nX/2,nX/2) oneDiagNx];
B_tr=[.5*dt^2*oneDiagNx
    dt*oneDiagNx];

for i=1:n
    xPVdyn=A_tr * xPV(:,i) + B_tr*uhist(:,i) ;% + [0;0;-g*dt^2/2 ; 0;0;-g*dt];
    xPV=[xPV xPVdyn];
end

trajectoryErrorFromTruth = xPV(1:nX,2:end)-xhist;
ceqDyn=reshape(trajectoryErrorFromTruth,[n*nX,1]);

if nX>=4
    cineqAvoidCircle=[];
    for i=1:n
        circCent=zeros(nX/2,1); circCent(1:2)=[5;5];
        cineqAvoidCircle=[cineqAvoidCircle; 2-norm(xhist(1:2,i)-circCent)];
    end
    cineq=cineqAvoidCircle;
else
    cineq=[];
end
    
ceq=ceqDyn;


end

