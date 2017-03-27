function q1q2q3q4_dot = kumarV2_Q(t,q1q2q3q4)
%note: P1,P2,P3 go forwards in time
%Q1,Q2,Q3,Q4 go backwards in time, hence the odd signage in derivatives

kumarV2Const
nn=length(q1q2q3q4);

global tP_GG PP1_GG PP2_GG PP3_GG
tP=tP_GG;
PP1=PP1_GG;
PP2=PP2_GG;
PP3=PP3_GG;
tqLen=length(tP);
P1=zeros(nX,nX);
P2=zeros(nX,nX);
P3=zeros(nX,nX);
for i=1:nX
    for j=1:nX
        P1(i,j)=interp1(tP,reshape(PP1(i,j,:),[tqLen,1]),t);
        P2(i,j)=interp1(tP,reshape(PP2(i,j,:),[tqLen,1]),t);
        P3(i,j)=interp1(tP,reshape(PP3(i,j,:),[tqLen,1]),t);
    end
end

q1Evec=q1q2q3q4(1:1*nX^2);
q2Evec=q1q2q3q4(1*nX^2+1:2*nX^2);
q3Evec=q1q2q3q4(2*nX^2+1:3*nX^2);
q4Evec=q1q2q3q4(3*nX^2+1:4*nX^2);
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
Q1dot = -Q1*F-F'*Q1+Q1*G1*invR1*G1'*Q1-Q1*G2*inv(R22)*G2'*Q3+H2'*invV2*H2*P2*Q3'-H2'*invV2*H2*P2*Q4;

Q2dot = -Q2*F-F'*Q2-Q2*G2*invR2*G2'*Q2+Q2*G1*invR1*G1'*Q1+Q1'*G1*invR1*G1'*Q2-Q1'*G1*invR1*R21*invR1*G1'*Q1;

Q3dot = -Q1dot-(Q1+Q3)*F-F'*(Q1+Q3)+(Q1+Q3)'*G1*invR1*G1'*(Q1+Q3)+Q3*P2*H2'*invV2*H2+H2'*invV2*H2*P2*Q3';

Q4dot = (-F'+Q1'*G1*invR1*G1'-Q2*G2*invR2*G2'+H2'*invV2*H2*P2)*Q4+Q4*(-F+G1*invR1*G1'*Q1-G2*invR2*G2'*Q2+P2*H2'*invV2*H2)+...
    +Q3'*G1*invR1*G1'*Q3+Q2*G2*invR2*R12*invR2*G2'*Q2+Q3'*G2*invR2*G2'*Q2+Q2*G2*invR2*G2'*Q3;

Q1dot=-Q1dot;
Q2dot=-Q2dot;
Q3dot=-Q3dot;
Q4dot=-Q4dot;

q1dot=reshape(Q1dot,[nX^2,1]);
q2dot=reshape(Q2dot,[nX^2,1]);
q3dot=reshape(Q3dot,[nX^2,1]);
q4dot=reshape(Q4dot,[nX^2,1]);

q1q2q3q4_dot=[q1dot;q2dot;q3dot;q4dot];


end

