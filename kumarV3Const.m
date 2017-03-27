%Contains constant matrices for Kumar1980 solution when roles are swapped

kumarV2Const;

%sign changes caused by "roleswapping"
R21Temp = -R21;
R12Temp = -R12;
R11Temp = -R22;
R22Temp = -R11;
V2Temp = V2;
V1Temp = V1;
G1Temp = -G1;
G2Temp = -G2;
Q1fTemp = Q2f;
Q2fTemp = Q1f;

eyeHalf=eye(nX/2);
zeroHalf=zeros(nX/2,nX/2);
F=[zeroHalf eyeHalf;zeroHalf zeroHalf];
G1=G2Temp;
G2=G1Temp;
umaxPur=1;
umaxEva=1;

R21=R12Temp;
R12=R21Temp;
R11=R22Temp;
R22=R11Temp;

W=.1*eye(nX); %process noise covariance
H1=eye(nX);
H2=H1;
V1=V2Temp; %measurement noise covariance
V2=V1Temp;

P0=.1*eye(nX);
Q1f=Q2fTemp;
Q2f=Q1fTemp;

tint=1/.4; %estimated interval between updates in ODE for noise calc




