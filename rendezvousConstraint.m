function [cineq,ceq] = rendezvousConstraint(uhist)

global PPmat_GG xest_GG ePV_GG r1_GG r0_GG QQ_GG;
PP=PPmat_GG;
xcurr=xest_GG;
ePVe=ePV_GG;
r1=r1_GG;
r0=r0_GG;
QQ=QQ_GG;

[nX,nU,~,~,~,dt]=get_problemSpecs;
RR=@(d) r0*diag([ones(nX/2,1)])+r1*[abs(d(1)) 0;0 abs(d(2))];
%RR2=@(d) r0*diag([ones(nX/2,1)])+r1*eye(2)*norm(d);

% [xhist,xPVhist]=xPVsim(uhist,xcurr);
% n=length(uhist)/nU;
% circcent=[5;5];
% circeq=[];
% if n>=2
%     for i=2:n %simple projection stuff
%         %http://stackoverflow.com/questions/1073336/circle-line-segment-collision-detection-algorithm
%         AB=xhist(:,i-1)-xhist(:,i);
%         AC=xhist(:,i-1)-circcent;
%         a1=dot(AC,AB)/norm(AB);
%         a2=AC-a1;
%         circeq=-a2+2;
%     end
% end

cineq=[];

ceq=rendezvous(uhist);

end










