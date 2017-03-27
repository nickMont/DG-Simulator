function [xhist,xPVhist] = xPVsim(uhist,xPV0)
%uhist is n blocks of size nUx1 as a stacked vector (n*nU x 1)
%xPV0 is nX x 1  initial PV
%xhist is nX/2 x n  positions
%xPVhist is nX*n x 1  position+velocity

[nX,nU,~,~,~,dt]=get_problemSpecs;
xcurr=xPV0;

oneDiagNx=diag(ones(nX/2,1));
A_tr=[oneDiagNx dt*oneDiagNx
    zeros(nX/2,nX/2) oneDiagNx];
B_tr=[.5*dt^2*oneDiagNx
    dt*oneDiagNx];

n=length(uhist)/nU; %num timesteps
xPV=xcurr;
xhist=[];
xPVhist=[];

for i=1:n
    xPV=A_tr*xPV+B_tr*uhist((i-1)*nU+1:(i-1)*nU+nU);
    xhist=[xhist xPV(1:nX/2)];
    xPVhist=[xPVhist; xPV];
end

end

