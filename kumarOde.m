function p1p2p3q1q2q3q4_dot = kumarOde(t,p1p2p3q1q2q3q4)
%note: P1,P2,P3 go forwards in time
%Q1,Q2,Q3,Q4 go backwards in time, hence the odd signage in derivatives

kumarConst
nn=length(p1p2p3q1q2q3q4);

p1Evec=p1p2p3q1q2q3q4(1:nX^2);
p2Evec=p1p2p3q1q2q3q4(nX^2+1:2*nX^2);
p3Evec=p1p2p3q1q2q3q4(2*nX^2+1:3*nX^2);
q1Evec=p1p2p3q1q2q3q4(3*nX^2+1:4*nX^2);
q2Evec=p1p2p3q1q2q3q4(4*nX^2+1:5*nX^2);
q3Evec=p1p2p3q1q2q3q4(5*nX^2+1:6*nX^2);
q4Evec=p1p2p3q1q2q3q4(6*nX^2+1:7*nX^2);
P1=reshape(p1Evec,[nX,nX]);
P2=reshape(p2Evec,[nX,nX]);
P3=reshape(p3Evec,[nX,nX]);
Q1=reshape(q1Evec,[nX,nX]);
Q2=reshape(q2Evec,[nX,nX]);
Q3=reshape(q3Evec,[nX,nX]);
Q4=reshape(q4Evec,[nX,nX]);


invV1=inv(V1);
invV2=inv(V2);
invR1=inv(R11);
invR2=inv(R22);

% %METHOD 1: Hardcoded jacobians (has stiffness issues)
% P1dot=F*P1+P1*F'-P1*H1'*invV1*H1*P1-P1*H2'*invV2*H2*P1+W;
% 
% P2dot=F*P2+P2*F'-P2*H2'*invV2*H2*P2+G1*invR1*G1'*(Q1+Q3)*(P2-P3)+(P2-P3')*(Q1+Q3)'*G1*invR1*G1'+W;
% 
% P3dot=F*P3+P3*F'-P1*H1'*invV1*H1*P3-P1*H2'*invV2*H2*P3-P1*(Q1+Q3)'*G1*invR1*G1'+P3*(Q1+Q3)'*G1*invR1*G1'+...
%     -P3*H2'*invV2*H2*P2+P1*H2'*invV2*H2*P2+W; %#ok<*MINV>
% 
% Q1dot=+Q1*F+F'*Q1+Q1*G1*invR1*G1'*Q1-Q1*G2*inv(R22)*G2'*Q3-H2'*invV2*H2*P2*Q3'+H2'*invV2*H2*P2*Q4;
% 
% Q3dot=-Q1dot+(Q1+Q3)*F+F*(Q1+Q3)+(Q1+Q3)'*G1*invR1*G1'*(Q1+Q3)-Q3*P2*H2'*invV2*H2-H2'*invV2*H2*P2*Q3';
% 
% Q2dot=+Q2*F+F'*Q2-Q2*G2*invR2*G2'*Q2+Q2*G1*invR1*G1'*Q1+Q1'*G1*invR1*G1'*Q2-Q1'*G1*invR1*R21*invR1*G1'*Q1;
% 
% Q4dot=-(-F'-Q1'*G1*invR1*G1'+Q2*G2*invR2*G2'+H2'*invV2*H2*P2)*Q4-Q4*(-F-G1*invR1*G1'*Q1+G2*invR2*G2'*Q2+P2*H2'*invV2*H2)+...
%     +Q3'*G1*invR1*G1'*Q3+Q2*G2*invR2*R12*invR2*G2'*Q2+Q3'*G2*invR2*G2'*Q2+Q2*G2*invR2*G2'*Q3;



% %METHOD 2: Negation of derivatives (inexact)
P1dot = F*P1+P1*F'-P1*H1'*invV1*H1*P1-P1*H2'*invV2*H2*P1+W;

P2dot = F*P2+P2*F'-P2*H2'*invV2*H2*P2-G1*invR1*G1'*(Q1+Q3)*(P2-P3)-(P2-P3')*(Q1+Q3)'*G1*invR1*G1'+W;

P3dot = F*P3+P3*F'-P1*H1'*invV1*H1*P3-P1*H2'*invV2*H2*P3+P1*(Q1+Q3)'*G1*invR1*G1'-P3*(Q1+Q3)'*G1*invR1*G1'+...
    -P3*H2'*invV2*H2*P2+P1*H2'*invV2*H2*P2+W; %#ok<*MINV>

Q1dot = -Q1*F-F'*Q1+Q1*G1*invR1*G1'*Q1-Q1*G2*inv(R22)*G2'*Q3+H2'*invV2*H2*P2*Q3'-H2'*invV2*H2*P2*Q4;

Q2dot = -Q2*F-F'*Q2-Q2*G2*invR2*G2'*Q2+Q2*G1*invR1*G1'*Q1+Q1'*G1*invR1*G1'*Q2-Q1'*G1*invR1*R21*invR1*G1'*Q1;

Q3dot = -Q1dot-(Q1+Q3)*F-F'*(Q1+Q3)+(Q1+Q3)'*G1*invR1*G1'*(Q1+Q3)+Q3*P2*H2'*invV2*H2+H2'*invV2*H2*P2*Q3';

Q4dot = (-F'+Q1'*G1*invR1*G1'-Q2*G2*invR2*G2'+H2'*invV2*H2*P2)*Q4+Q4*(-F+G1*invR1*G1'*Q1-G2*invR2*G2'*Q2+P2*H2'*invV2*H2)+...
    +Q3'*G1*invR1*G1'*Q3+Q2*G2*invR2*R12*invR2*G2'*Q2+Q3'*G2*invR2*G2'*Q2+Q2*G2*invR2*G2'*Q3;

Q1dot=-Q1dot;
Q2dot=-Q2dot;
Q3dot=-Q3dot;
Q4dot=-Q4dot;

p1dot=reshape(P1dot,[nX^2,1]);
p2dot=reshape(P2dot,[nX^2,1]);
p3dot=reshape(P3dot,[nX^2,1]);
q1dot=reshape(Q1dot,[nX^2,1]);
q2dot=reshape(Q2dot,[nX^2,1]);
q3dot=reshape(Q3dot,[nX^2,1]);
q4dot=reshape(Q4dot,[nX^2,1]);

p1p2p3q1q2q3q4_dot=[p1dot;p2dot;p3dot;q1dot;q2dot;q3dot;q4dot];


end

