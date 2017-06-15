function [x1,P1] = pointMassKFStep(x0,P0,Q0,u1,u2,z,H,R0,dt)
global cddGG
cdd=cddGG; %drag coefficient

nX=length(x0);

eyeHalfNX=eye(nX/2);
zerosHalfNX=zeros(nX/2,nX/2);
A_tr=[eyeHalfNX (1*dt-cdd*dt^2/2)*eyeHalfNX; zerosHalfNX (1-cdd*dt)*eyeHalfNX];
B_tr=[dt^2/2*eyeHalfNX; dt*eyeHalfNX];
Gammak=[dt^2/2*eyeHalfNX; dt*eyeHalfNX];

x1bar=A_tr*x0+B_tr*u2-B_tr*u1;
P1bar=A_tr*P0*A_tr'+Gammak*Q0*Gammak';
Sk=H*P1bar*H'+R0;
Kk=P1bar*H'*inv(Sk);
x1=x1bar+Kk*(z-H*x1bar);
P1=(eye(nX)-Kk*H)*P1bar;


end

