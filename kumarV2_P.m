function p123dot = kumarV2_P(t,p1p2p3)

kumarV2Const;

global tQ_GG QQ1_GG QQ2_GG QQ3_GG QQ4_GG flagRunV3
% if flagRunV3==1
%     derp=9001
%     kumarV3Const;
% end
tQ=tQ_GG;
QQ1=QQ1_GG;
QQ2=QQ2_GG;
QQ3=QQ3_GG;
QQ4=QQ4_GG;
tqLen=length(tQ);
Q1=zeros(nX,nX);
Q2=zeros(nX,nX);
Q3=zeros(nX,nX);
Q4=zeros(nX,nX);
for i=1:nX
    for j=1:nX
        Q1(i,j)=interp1(tQ,reshape(QQ1(i,j,:),[tqLen,1]),t);
        Q2(i,j)=interp1(tQ,reshape(QQ2(i,j,:),[tqLen,1]),t);
        Q3(i,j)=interp1(tQ,reshape(QQ3(i,j,:),[tqLen,1]),t);
        Q4(i,j)=interp1(tQ,reshape(QQ4(i,j,:),[tqLen,1]),t);
    end
end

invV1=inv(V1);
invV2=inv(V2);
invR1=inv(R11);
invR2=inv(R22);

P1=reshape(p1p2p3(1:nX^2),[nX,nX]);
P2=reshape(p1p2p3(nX^2+1:2*nX^2),[nX,nX]);
P3=reshape(p1p2p3(2*nX^2+1:3*nX^2),[nX,nX]);

P1dot = F*P1+P1*F'-P1*H1'*invV1*H1*P1-P1*H2'*invV2*H2*P1+W;

P2dot = F*P2+P2*F'-P2*H2'*invV2*H2*P2-G1*invR1*G1'*(Q1+Q3)*(P2-P3)-(P2-P3')*(Q1+Q3)'*G1*invR1*G1'+W;

P3dot = F*P3+P3*F'-P1*H1'*invV1*H1*P3-P1*H2'*invV2*H2*P3+P1*(Q1+Q3)'*G1*invR1*G1'-P3*(Q1+Q3)'*G1*invR1*G1'+...
    -P3*H2'*invV2*H2*P2+P1*H2'*invV2*H2*P2+W; %#ok<*MINV>

p123dot=[reshape(P1dot,[nX^2,1]); reshape(P2dot,[nX^2,1]); reshape(P3dot,[nX^2,1])];

end

