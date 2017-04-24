function mmdot = rhodesV2_M(t,mm)

global tPN_GG PP_GG NN_GG
tPN=tPN_GG;
PP=PP_GG;
NN=NN_GG;
rhodesV2Const;
M=reshape(mm,[nX,nX]);

P=zeros(nX,nX);
N=zeros(nX,nX);
tpnLen=length(tPN);
for i=1:nX
    for j=1:nX
        P(i,j)=interp1(tPN,reshape(PP(i,j,:),[tpnLen,1]),t);
        N(i,j)=interp1(tPN,reshape(NN(i,j,:),[tpnLen,1]),t);
    end
end

A=F-G1*inv(R11)*G1'*(P+N);

Mdot=A*M+M*A'-M*H2'*inv(V2)*H2*M;
mmdot=reshape(Mdot,[nX^2,1]);    

end

