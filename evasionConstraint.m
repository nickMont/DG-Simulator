function [cineq,ceq] = evasionConstraint(uhist)

global PPmat_GG xEva_GG r1_GG r0_GG QQ_GG n0_GG nmax_GG xPurEst_GG uPurGuess_GG;
PP=PPmat_GG;
yEva=xEva_GG;
xPur=xPurEst_GG;
r1=r1_GG;
r0=r0_GG;
QQ=QQ_GG;
nmax=nmax_GG;
n0=n0_GG;
uPurHist=uPurGuess_GG;

[nX,nU,~,~,~,dt]=get_problemSpecs;

%assumes nX/2 = nU
oneDiagNx=diag(ones(nX/2,1));
A_tr=[oneDiagNx dt*oneDiagNx
    zeros(nX/2,nX/2) oneDiagNx];
B_tr=[.5*dt^2*oneDiagNx
    dt*oneDiagNx];

safeRegion=2*ones(1+nmax-n0,1);  %minimum safe distance, can expand with Rxx
%safeRegion(1) corresponds to t=tcurrent

evasionCon=[];

yEvaHist=xPVsim(uhist,yEva);
xPurHist=xPVsim(uPurHist,xPur);

if nmax-n0+1>=2
    for i=2:nmax-n0+1
        evasionCon=[evasionCon; safeRegion(i)-norm(yEvaHist(:,i)-xPurHist(:,i))];
    end
end


cineq=evasionCon;

ceq=[];

end










