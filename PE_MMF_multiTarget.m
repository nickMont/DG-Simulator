clear;clc;
%IMM for determining the target that an evader is trying to escape to
%NOTE: Model 2 is a min-distance one-step evasion strategy

nSim=10;

nX=4;
nU=2;
nZ=2;
dt=.5;
cd=.25; %drag coefficient, =0 to ignore
A_tr=[eye(nX/2) (1*dt-cd*dt^2/2)*eye(nX/2); zeros(nX/2,nX/2) (1-cd*dt)*eye(nX/2)];
B_tr=[dt^2/2*eye(nX/2); dt*eye(nX/2)];
Gammak=[dt^2/2*eye(nX/2); dt*eye(nX/2)];
H=[eye(nX/2) zeros(nX/2)];
P0=.2*eye(nX);
Q0=.05*eye(nX/2);
R0=.1*eye(nX/2);
PE=P0; PP=P0;

maxDist=2;

nmod_pur=1; %number of models considered
q1_pur=1.5;
q2_pur=[0 .5 1];
r1_pur=1.5;
ntrue_pur=2; %index of true model
MijE = .25*ones(nmod_pur,nmod_pur); %probability of mode switching

numTargets=2;
targetsList=zeros(nX/2,numTargets);
targetsList=[100 110
             50 75];
numControls=2;
ntrueEvaIndex=1*ones(nSim,1);
targetTrue_evaIndex=1*ones(nSim,1);
q1_eva=.05;
q2_eva=0;
r1_eva=[.1 3.0];
nmod_eva=numTargets*numControls; %number of control types considered
MijP = .25*ones(nmod_eva,nmod_eva); %probability of mode switching

flagBreakOnFlip=0;  %stop simulation on (estimated) collision if flag==1
flagUseDetm=1;  %zero noise if flag==1
flagUseXPrevInsteadOfXModel=1;  %use overall xhat estimate instead of model-
    %specific xhat estimate in MMKF if flag==1
    %Best performance at flag==0.
flagUseXPrevGroupedByModel=1;  %use a priori xhat estimates from the same 
    %target.  Essentially flagUseXPrevInsteadOfXModel==1 but with only
    %averaging over mu's corresponding to the same target
    %ONLY TRIGGERED IF FLAGUSEXPREV==1
%HANDLES BEST IF GROUPS==1
    
mu_min=10^-4;
normpdf_diag_DEBUG=.01; %with high confidence, diagonal elements of normpdf matrix go to zero

xEva=[5;50;0;0];
xPur=[0;0;0;.1];
umaxPur=2.5;
umaxEva=2;

captureThresh=1.5;
flagSign1Flip=0;
flagSign2Flip=0;

%set randn's seed RNG
n_seed = round(sum(clock*100));
randn('state', n_seed);

xhat0=xEva-xPur+(chol(P0))'*randn(nX,1);
if flagUseDetm==1
    xhat0=xEva-xPur; %deterministic behavior
end
ehat_prev_eva=xhat0;

Pstore_eva=zeros(nX,nX,nmod_eva,nSim);
for jj=1:nmod_eva
    Pstore_eva(:,:,jj,1)=P0;
end
ehathist_eva=zeros(nX,nmod_eva,nSim);

mukhist_eva=zeros(nmod_eva,nSim);
LambdaVec_eva=zeros(nmod_eva,nSim);
mukhist_eva(:,1)=1/nmod_eva*ones(nmod_eva,1);
muPrev_eva=mukhist_eva(:,1);

%----- Simulation parameters
tkhist = [0:nSim]'*dt;

%----- Random number seed
n_seed = round(sum(clock*100));
randn('state', n_seed);


global umax_GG xCurr_GG distMax_GG dt_GG u_GG nX_GG cd_GG %#ok<*NUSED>
dt_GG=dt;
nX_GG=nX;
cd_GG=cd;

e_true=zeros(nX,nSim);

LambdaTemp=zeros(nmod_eva,1);

optionsForFMC=optimset('Display','off');
numSimRuns=1;

%generate all possible combinations of targets and controls
possibleCombos=combvec(1:1:numControls,1:1:numTargets);

for i=1:nSim
    numSimRuns=numSimRuns+1;
    i
    
    if i==1
        xEvaEst=xEva;
        xhat=xhat0;
        xPurSelfEst=xPur;
        for j=1:nmod_eva
            ehathist_eva(:,j,1)=xhat0;
        end
    else
        xhat=ehat_prev_eva;
    end
    
    xhatE=xhat;
    xhatP=xhat;
    
    ntrue_eva=ntrueEvaIndex(i);
    targetTrue_eva=targetTrue_evaIndex(i);
    
    %Pursuer control
    umax_GG=umaxPur;        %#ok<NASGU>
    q2j=q2_pur(ntrue_pur); %index of true control QQ matrix
    QQ=[q1_pur*eye(nX/2) zeros(nX/2,nX/2); zeros(nX/2,nX/2) q2j*eye(nX/2)];
    RR=r1_pur*eye(nU);
    JJ=@(u) (A_tr*xhatP-B_tr*u)'*QQ*(A_tr*xhatP-B_tr*u) + u'*RR*u;
    [uPurTrue,JcP]=fmincon(JJ,zeros(nU,1),[],[], [],[], [],[], 'oneStepMaxT_constraint',optionsForFMC);
    
    %Evader control
    thisTarget=targetsList(:,targetTrue_eva);
    umax_GG=umaxEva;
    r1j=r1_eva(ntrue_eva); %index of true control QQ matrix
    QQ = [q1_eva*eye(nX/2) zeros(nX/2,nX/2); zeros(nX/2,nX/2) q2_eva*eye(nX/2)];
    RR = r1j*eye(nU);
    targetDist=[thisTarget;zeros(nX/2,1)]-xEvaEst;
    JJ=@(u) (A_tr*targetDist-B_tr*u)'*QQ*(A_tr*targetDist-B_tr*u) + u'*RR*u;
    [uEvaTrue,JcE]=fmincon(JJ,zeros(nU,1),[],[], [],[], [],[], 'oneStepMaxT_constraint',optionsForFMC);

    
    if flagUseDetm==1
        noiseX1=zeros(nX/2,1);
        noiseX2=noiseX1;
    else
        noiseX1=chol(Q0)*randn(nX/2,1);
        noiseX2=chol(Q0)*randn(nX/2,1);
    end
    
    xPur=A_tr*xPur+B_tr*uPurTrue+Gammak*noiseX1;
    xEva=A_tr*xEva+B_tr*uEvaTrue+Gammak*noiseX2;
    
    %Evader measurement of own state 
    if flagUseDetm==1
        noiseZ=zeros(nZ,1);
    else
        noiseZ=chol(R0)*randn(nZ,1);
    end
    zE=H*xEva+noiseZ;
    KK=PE*H'*(H*PE*H'+R0)^-1;
    xEvaEst=(A_tr*xEvaEst+B_tr*uEvaTrue)+KK*(zE-H*(A_tr*xEvaEst+B_tr*uEvaTrue));
    PE=(eye(nX)-KK*H)*PE*(eye(nX)-KK*H)'+KK*R0*KK';
    
    %Pursuer measurement of game state
    if flagUseDetm==1
        noiseZ=zeros(nZ,1);
    else
        noiseZ=chol(R0)*randn(nZ,1);
    end
    z=H*(xEva-xPur)+noiseZ;
    
    %Pursuer's location
    if flagUseDetm==1
        noiseZ=zeros(nZ,1);
    else
        noiseZ=chol(R0)*randn(nZ,1);
    end
    xPurPrev_SelfEst=xPurSelfEst;
    zP=H*xPur+noiseZ; KK=PP*H'*(H*PP*H'+R0)^-1;
    xPurSelfEst=(A_tr*xPurSelfEst-B_tr*uPurTrue)+KK*(zE-H*(A_tr*xPurSelfEst-B_tr*uPurTrue));
    PP=(eye(nX)-KK*H)*PP*(eye(nX)-KK*H)'+KK*R0*KK'; %Pursuer's location
    
    LambdaTemp=zeros(nmod_eva,1);
    for j=1:nmod_eva
        
        r1j=r1_eva(possibleCombos(1,j));
        targetInd=possibleCombos(2,j);
        thisTargetPoss=targetsList(:,targetInd);
        
        ehatkj_beforeMotion=ehathist_eva(:,j,i);
        if flagUseXPrevInsteadOfXModel==1
            if flagUseXPrevGroupedByModel==1
                ehatkj_beforeMotion=zeros(nX,1);
                modelInd=find(possibleCombos(2,:)==targetInd);
                muReNorm=sum(muPrev_eva(modelInd));
                for jj=1:length(modelInd)
                    ehatkj_beforeMotion=ehatkj_beforeMotion + muPrev_eva(jj)*ehathist_eva(:,jj,i)/muReNorm;
                end
            else
                ehatkj_beforeMotion=ehat_prev_eva;
            end
        end
        
        Pk=Pstore_eva(:,:,j,i);

        targetDist=[thisTargetPoss;zeros(nX/2,1)]-(ehatkj_beforeMotion+xPurPrev_SelfEst);
        
        umax_GG=umaxEva;
        QQ = [q1_eva*eye(nX/2) zeros(nX/2,nX/2); zeros(nX/2,nX/2) zeros(nX/2,nX/2)];
        RR = r1j*eye(nU);
        JJ=@(u) (A_tr*targetDist-B_tr*u)'*QQ*(A_tr*targetDist-B_tr*u) + u'*RR*u;
        [uEva,Jc]=fmincon(JJ,zeros(nU,1),[],[], [],[], [],[], 'oneStepMaxT_constraint',optionsForFMC)
        
        
        %propagate state differences
        Qk=Q0;
        R=R0;
        ebarj = A_tr*ehatkj_beforeMotion + B_tr*uEva - B_tr*uPurTrue;
        Pbarj = A_tr*Pk*A_tr'+Gammak*Qk*Gammak';
        
        %KF
        nuj=z-H*ebarj;
        Sk=H*Pbarj*H'+R;
        Sk_inv=inv(Sk);
        Wk=Pbarj*H'*Sk_inv;
        
        xhatkj=ebarj+Wk*nuj;
        Pj=Pbarj-Wk*Sk*Wk';
        
        PPstore(:,:,j,i+1)=Pj; %#ok<SAGROW>
        ehathist_eva(:,j,i+1)=xhatkj;
        
        muNu=zeros(nZ,1);
        normpdf_Eval = mvnpdf(nuj,muNu,Sk); %evaluates multivariate normal
        Lambdaj=normpdf_Eval;
        LambdaVec_eva(i,j)=Lambdaj;
        LambdaTemp(j)=Lambdaj;
    end
    
    muPrev_eva=mukhist_eva(:,i);
    muStackTemp=zeros(nmod_eva,1);
    for j=1:nmod_eva
        muStackTemp(j)=LambdaTemp(j)*muPrev_eva(j)/dot(LambdaTemp,muPrev_eva);
    end
    for j=1:nmod_eva
        if muStackTemp(j)<=mu_min
            muStackTemp(j)=mu_min;
        end
    end
    muPrev_eva=muStackTemp/sum(muStackTemp)
    mukhist_eva(:,i+1)=muPrev_eva;
    
    ehat_weighted_eva=zeros(nX,1);
    for j=1:nmod_eva
        ttt=muPrev_eva(j)*ehathist_eva(:,j,i+1);
        ehat_weighted_eva=ehat_weighted_eva+muPrev_eva(j)*ehathist_eva(:,j,i+1);
    end
    
    %If captured at a time step or thresholds have flipped
    if norm(xEva(1:nX/2)-xPur(1:nX/2))<=captureThresh
        break
    end
    if ehat_weighted_eva(1)*ehat_prev_eva(1)<=0
        flagSign1Flip=1;
    end
    if ehat_weighted_eva(2)*ehat_prev_eva(2)<=0
        flagSign1Flip=1;
    end
    if flagSign1Flip==1 && flagSign2Flip==1 && flagBreakOnFlip==1
        break
    end
    
    ehat_prev_eva=ehat_weighted_eva;
    


end

figure(1);clf;
colors='brgybrgybrgy';
for j=1:nmod_eva
    plot(tkhist(1:numSimRuns),mukhist_eva(j,1:numSimRuns),colors(j))
    hold on
end
axis([0 tkhist(end) 0 1.1])
legend('r1','r2','r3','r4')



