%Contains constant matrices for Kumar1980 solution
global nX_GG nU_GG
nX=nX_GG;
nU=nU_GG;

eyeHalf=eye(nX/2);
zeroHalf=zeros(nX/2,nX/2);
F=[zeroHalf eyeHalf;zeroHalf zeroHalf];
G1=[zeroHalf;eyeHalf];
G2=[zeroHalf;eyeHalf];

R21=zeroHalf;
R12=zeroHalf;
R11=.2*eyeHalf;
R22=.1*eyeHalf;

W=.1*eye(nX); %process noise covariance
H1=eye(nX);
H2=H1;
V1=.5*eye(nX); %measurement noise covariance
V2=.5*eye(nX);

P0=.5*eye(nX);
Q1f=10*eye(nX);
Q2f=5*eye(nX);

