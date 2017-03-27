function J = J_PE_I_pur(u)
%pursuit-evasion with information

[nX,nU,~,~,~,dt]=get_problemSpecs; %#ok<ASGLU>
global xest_GG ePV_GG
xcurr=xest_GG;
ePVe=ePV_GG;

n=length(u)/nU; %num timesteps
ntot=n;

an=.1;

Rfull=.25*eye(nU*n);

Ju=u'*Rfull*u;
Jinfo=J_logDetY(u(1:nU));
J=Jinfo+Ju+an*ntot*dt;

end

