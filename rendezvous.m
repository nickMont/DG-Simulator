function violation = rendezvous(uhist)

[nX,nU,~,~,~,dt]=get_problemSpecs;
global xest_GG ePV_GG
xcurr=xest_GG;
ePVe=ePV_GG;

oneDiagNx=diag(ones(nX/2,1));
A_tr=[oneDiagNx dt*oneDiagNx
    zeros(nX/2,nX/2) oneDiagNx];
B_tr=[.5*dt^2*oneDiagNx
    dt*oneDiagNx];


n=length(uhist)/nU; %num timesteps
xPV=xcurr;
ePV=ePVe;
for i=1:n
    ePV=A_tr*ePV;
    xPV=A_tr*xPV+B_tr*uhist((i-1)*nU+1:(i-1)*nU+nU);
end

violation=ePV(1:2)-xPV(1:2);

end

