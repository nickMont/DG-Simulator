%Contains constant matrices for Kumar1980 solution
global nX_GG
nX=nX_GG;

qPurF=10;
rPur=5;
qEvaF=10;
rEva=20;

eyeHalf=eye(nX/2);
zeroHalf=zeros(nX/2,nX/2);
F=[zeroHalf eyeHalf;zeroHalf zeroHalf];
G1=[zeroHalf;eyeHalf];
G2=[zeroHalf;eyeHalf];
umaxPur=1;
umaxEva=1;

R21=zeroHalf;
R12=zeroHalf;
R11=rPur*eyeHalf;
R22=rEva*eyeHalf;

W=.1*eye(nX); %process noise covariance
H1=eye(nX);
H2=H1;
V1=.05*eye(nX); %measurement noise covariance
V2=.05*eye(nX);

P0=.1*eye(nX);
Q1f=qPurF*[eyeHalf zeroHalf; zeroHalf zeroHalf];
Q2f=qEvaF*[eyeHalf zeroHalf; zeroHalf zeroHalf];

tint=1/.4; %estimated interval between updates in ODE for noise calc




