function pnmdot = rhodesV2Ode(t,PN)

rhodesV2Const;

global tM_GG MM_GG
tM=tM_GG;
MM=MM_GG;
M=zeros(nX,nX);
tmLen=length(tM);
for i=1:nX
    for j=1:nX
        dm=reshape(MM(i,j,:),[tmLen,1]);
        M(i,j)=interp1(tM,dm,t);
    end
end
        
pEvec=PN(1:nX^2);
nEvec=PN(nX^2+1:2*nX^2);

P=reshape(pEvec,[nX,nX]);
N=reshape(nEvec,[nX,nX]);

invV1=inv(V1);
invV2=inv(V2);
invR1=inv(R11);
invR2=inv(R22);

Pdot=-(P*F+F'*P-P*(G1*invR1*G1'+G2*invR2*G2')*P);
Ndot=-(N*F+F'*N+P*(G1*invR1*G1'+G2*invR2*G2')*P-(P+N)*G1*invR1*G1'*(P+N)+...
    -N*M*H2'*invV2*H2-H2'*invV2*H2*M*N);

Pdot=-Pdot;
Ndot=-Ndot;


%note the negative from the jacobian
pnmdot = [reshape(Pdot,[nX^2,1]); reshape(Ndot,[nX^2,1])];



end

