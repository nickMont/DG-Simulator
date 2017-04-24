function pnmdot = rhodesOde(t,PNM)

rhodesConst;

pEvec=PNM(1:nX^2);
nEvec=PNM(nX^2+1:2*nX^2);
mEvec=PNM(2*nX^2+1:end);

P=reshape(pEvec,[nX,nX]);
N=reshape(nEvec,[nX,nX]);
M=reshape(mEvec,[nX,nX]);

invV1=inv(V1);
invV2=inv(V2);
invR1=inv(R11);
invR2=inv(R22);

A=F-G1*invR1*G1'*(P+N);

Pdot=-(P*F+F'*P-P*(G1*invR1*G1'+G2*invR2*G2')*P);
Ndot=-(N*F+F'*N+P*(G1*invR1*G1'+G2*invR2*G2')*P-(P+N)*G1*invR1*G1'*(P+N)+...
    -N*M*H2'*invV2*H2-H2'*invV2*H2*M*N);
Mdot=A*M+M*A'-M*H2'*invV2*H2*M;

Pdot=-Pdot;
Ndot=-Ndot;
Mdot= Mdot; %defined forwards in time

%note the negative from the jacobian
pnmdot = [reshape(Pdot,[nX^2,1]); reshape(Ndot,[nX^2,1]); reshape(Mdot,[nX^2,1])];



end

