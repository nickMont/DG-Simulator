function Xdot = kumarV2Game(t,X)

kumarV2Const;

x1hat=X(1:nX);
x2hat=X(nX+1:2*nX);
xTrue=X(2*nX+1:3*nX);

global tmax_GG tprev_GG uPstore_GG uEstore_GG J1_GG J2_GG flagRunV3 flagNoSpies%#ok<NUSED>

if flagRunV3==1
    global tQe QQ2Pur QQ2Eva tQp %#ok
    tqLen=length(tQe);
    tqpLen=length(tQp);
    for i=1:nX
        for j=1:nX
            Q2(i,j)=interp1(tQe,reshape(QQ2Eva(i,j,:),[tqLen,1]),t);
            Q2p(i,j)=interp1(tQp,reshape(QQ2Pur(i,j,:),[tqpLen,1]),t);
        end
    end
    
    
    % up=vectorSaturationF(-inv(R11)*G1'*(Q1*x1hat+Q3*(x1hat-x2hat)),0,umaxPur)
    up=vectorSaturationF(-inv(R11)*G1'*Q2p*x1hat,0,umaxPur)
    ue=vectorSaturationF(-inv(R22)*G2'*Q2*x2hat,0,umaxEva)
    J1_GG=J1_GG+up'*R11*up*(t-tprev_GG)^2;
    J2_GG=J2_GG+ue'*R22*ue*(t-tprev_GG)^2;
    
    %uPstore_GG=up;
    %uEstore_GG=ue;
    
    xTrueDot=F*xTrue+G1*up-G2*ue+chol(W)*randn(nX,1);
    x1hatDot=F*x1hat+G1*up-G2*ue+chol(V1)*randn(nX,1);
    x2hatDot=F*x2hat+G1*up-G2*ue+chol(V2)*randn(nX,1);
    
    Xdot=[xTrueDot;x1hatDot;x2hatDot];
    tprev_GG=t;
else
    global tQ_GG QQ1_GG QQ2_GG QQ3_GG %#ok
    tQ=tQ_GG;
    QQ1=QQ1_GG;
    QQ2=QQ2_GG;
    QQ3=QQ3_GG;
    tmax=tmax_GG;
    
    
    Q1=zeros(nX,nX);
    Q2=zeros(nX,nX);
    Q3=zeros(nX,nX);
    tqLen=length(tQ);
    for i=1:nX
        for j=1:nX
            Q1(i,j)=interp1(tQ,reshape(QQ1(i,j,:),[tqLen,1]),t);
            Q2(i,j)=interp1(tQ,reshape(QQ2(i,j,:),[tqLen,1]),t);
            Q3(i,j)=interp1(tQ,reshape(QQ3(i,j,:),[tqLen,1]),t);
        end
    end
    
    if flagNoSpies==1
        up=vectorSaturationF(-inv(R11)*G1'*(Q1*x1hat),0,umaxPur)
        ue=vectorSaturationF(-inv(R22)*G2'*Q2*x2hat,0,umaxEva)
    else
        up=vectorSaturationF(-inv(R11)*G1'*(Q1*x1hat+Q3*(x1hat-x2hat)),0,umaxPur)
        ue=vectorSaturationF(-inv(R22)*G2'*Q2*x2hat,0,umaxEva)
    end
    J1_GG=J1_GG+up'*R11*up*(t-tprev_GG)^2;
    J2_GG=J2_GG+ue'*R22*ue*(t-tprev_GG)^2;
    
    %uPstore_GG=up;
    %uEstore_GG=ue;
    
    xTrueDot=F*xTrue+G1*up-G2*ue+chol(W)*randn(nX,1);
    x1hatDot=F*x1hat+G1*up-G2*ue+chol(V1)*randn(nX,1);
    x2hatDot=F*x2hat+G1*vectorSaturationF(-inv(R11)*G1'*Q1*x2hat,-umaxPur,umaxPur) - G2*ue+chol(V2)*randn(nX,1);
    
    Xdot=[xTrueDot;x1hatDot;x2hatDot];
    tprev_GG=t;

end

end

