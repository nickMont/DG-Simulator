clear;clc;
%multiple model strategy determination, from the pursuer, of the evader
%NOTE: Model 2 is a min-distance one-step evasion strategy

%NEED TO UPDATE FROM NORMPDF TO MVNPDF

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

maxDist=2;

nmod_pur=2; %number of models considered
q1_pur=1.5;
q2_pur=[0 .5 1];
r1_pur=1.5;
ntrue_pur=2; %index of true model
MijE = .25*ones(nmod_pur,nmod_pur); %probability of mode switching

nmod_eva=2; %number of models considered
ntrueEvaIndex=ones(nSim,1);
q1_eva=2.5;
q2_eva=0;
r1_eva=[.1 3.0];
ntrue_eva=2; %index of true model
MijP = .25*ones(nmod_eva,nmod_eva); %probability of mode switching

flagBreakOnFlip=0;  %stop simulation on (estimated) collision if flag==1
flagUseDetm=1;  %zero initial noise if flag==1
flagUseXPrevInsteadOfXModel=0;  %use overall xhat estimate instead of model-
    %specific xhat estimate in MMKF if flag==1
    %Best performance at flag==0.
flagSwapAtTenEva=0;
flagSwapAtTwentyEva=0;  %swap time index

nSim=20;
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
xhat_prev_eva=xhat0;

Pstore_eva=zeros(nX,nX,nmod_eva,nSim);
Pstore_eva(:,:,1,1)=P0;
Pstore_eva(:,:,2,1)=P0;
xhathist_eva=zeros(nX,nmod_eva,nSim);

mukhist_eva=zeros(nmod_eva,nSim);
LambdaVec_eva=zeros(nmod_eva,nSim);
mukhist_eva(:,1)=1/nmod_eva*ones(nmod_eva,1);

%----- Simulation parameters
tkhist = [0:nSim]'*dt;

%----- Random number seed
n_seed = round(sum(clock*100));
randn('state', n_seed);


global umax_GG xCurr_GG distMax_GG dt_GG u_GG nX_GG cd_GG
dt_GG=dt;
nX_GG=nX;
cd_GG=cd;

e_true=zeros(nX,nSim);

LambdaTemp=zeros(nmod_eva,1);

optionsForFMC=optimset('Display','off');
numSimRuns=1;
for i=1:nSim
    numSimRuns=numSimRuns+1;
    i
    
    if i==1
        xhat=xhat0;
        for j=1:nmod_eva
            xhathist_eva(:,j,1)=xhat0;
        end
    else
        xhat=xhat_prev_eva;
    end
    
    xhatE=xhat;
    xhatP=xhat;
    
    ntrue_eva=ntrueEvaIndex(i);
    
    %Pursuer control
    umax_GG=umaxPur;        %#ok<NASGU>
    q2j=q2_pur(ntrue_pur); %index of true control QQ matrix
    QQ=[q1_pur*eye(nX/2) zeros(nX/2,nX/2); zeros(nX/2,nX/2) q2j*eye(nX/2)];
    RR=r1_pur*eye(nU);
    JJ=@(u) (A_tr*xhatP-B_tr*u)'*QQ*(A_tr*xhatP-B_tr*u) + u'*RR*u;
    [uPurTrue,JcP]=fmincon(JJ,zeros(nU,1),[],[], [],[], [],[], 'oneStepMaxT_constraint',optionsForFMC);
    
    %Evader control
    umax_GG=umaxEva;
    if ntrue_eva==1
        r1j=r1_eva(ntrue_eva); %index of true control QQ matrix
        QQ = -[q1_eva*eye(nX/2) zeros(nX/2,nX/2); zeros(nX/2,nX/2) q2_eva*eye(nX/2)];
        RR = r1j*eye(nU);
        JJ=@(u) (A_tr*xhatE-B_tr*u)'*QQ*(A_tr*xhatE-B_tr*u) + u'*RR*u;
        [uEvaTrue,JcE]=fmincon(JJ,zeros(nU,1),[],[], [],[], [],[], 'oneStepMaxT_constraint',optionsForFMC);
    elseif ntrue_eva==2
        xCurr_GG=xhatE;
        distMax_GG=maxDist;
        u_GG=uPurTrue;
        r1j=r1_eva(ntrue_eva); %index of true control QQ matrix
        QQ = -[q1_eva*eye(nX/2) zeros(nX/2,nX/2); zeros(nX/2,nX/2) q2_eva*eye(nX/2)];
        RR = r1j*eye(nU);
        JJ=@(u) (A_tr*xhatE-B_tr*u)'*QQ*(A_tr*xhatE-B_tr*u) + u'*RR*u;
        [uEvaTrue,JcE]=fmincon(JJ,zeros(nU,1),[],[], [],[], [],[], 'oneStep_maxTminD',optionsForFMC);
    end
    
    xPur=A_tr*xPur+B_tr*uPurTrue;
    xEva=A_tr*xEva+B_tr*uEvaTrue;
    
    z=H*(xEva-xPur);
    
    LambdaTemp=zeros(nmod_eva,1);
    for j=1:nmod_eva
        
        xhatkj_beforeMotion=xhathist_eva(:,j,i);
        if flagUseXPrevInsteadOfXModel==1
            xhatkj_beforeMotion=xhat_prev_eva;
        end
        Pk=Pstore_eva(:,:,j,i);
        
        umax_GG=umaxEva;
        r1j=r1_eva(j); %index of true control QQ matrix
        QQ = -[q1_eva*eye(nX/2) zeros(nX/2,nX/2); zeros(nX/2,nX/2) zeros(nX/2,nX/2)];
        RR = r1j*eye(nU);
        JJ=@(u) (A_tr*xhat-B_tr*u)'*QQ*(A_tr*xhat-B_tr*u) + u'*RR*u;
        [uEva,Jc]=fmincon(JJ,zeros(nU,1),[],[], [],[], [],[], 'oneStepMaxT_constraint',optionsForFMC);
        
        %propagate state differences
        Qk=Q0;
        R=R0;
        ebarj = A_tr*xhatkj_beforeMotion + B_tr*uEva - B_tr*uPurTrue;
        Pbarj = A_tr*Pk*A_tr'+Gammak*Qk*Gammak';
        
        %KF
        nuj=z-H*ebarj;
        Sk=H*Pbarj*H'+R;
        Sk_inv=inv(Sk);
        Wk=Pbarj*H'*Sk_inv;
        
        xhatkj=ebarj+Wk*nuj;
        Pj=Pbarj-Wk*Sk*Wk';
        
        PPstore(:,:,j,i+1)=Pj;
        xhathist_eva(:,j,i+1)=xhatkj;
        
        normpdf_Eval = normpdf(nuj,0,Sk); %returns a diagonal matrix of P(x1|z),P(x2|z),...
        for kk=1:nZ
            if normpdf_Eval(kk,kk)<=normpdf_diag_DEBUG
                normpdf_Eval(kk,kk)=normpdf_diag_DEBUG;
            end
        end
        Lambdaj=trace(normpdf_Eval);
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
    muPrev_eva=muStackTemp/sum(muStackTemp);
    mukhist_eva(:,i+1)=muPrev_eva;
    
    xhat_weighted_eva=zeros(nX,1);
    for j=1:nmod_eva
        ttt=muPrev_eva(j)*xhathist_eva(:,j,i+1);
        xhat_weighted_eva=xhat_weighted_eva+muPrev_eva(j)*xhathist_eva(:,j,i+1);
    end
    
    %If captured at a time step or thresholds have flipped
    if norm(xEva(1:nX/2)-xPur(1:nX/2))<=captureThresh
        break
    end
    if xhat_weighted_eva(1)*xhat_prev_eva(1)<=0
        flagSign1Flip=1;
    end
    if xhat_weighted_eva(2)*xhat_prev_eva(2)<=0
        flagSign1Flip=1;
    end
    if flagSign1Flip==1 && flagSign2Flip==1 && flagBreakOnFlip==1
        break
    end
    
    xhat_prev_eva=xhat_weighted_eva;
    


end

% figure(1);clf;
% colors='brgybrgybrgy';
% for j=1:nmod_eva
%     plot(tkhist(1:numSimRuns),mukhist_eva(j,1:numSimRuns),colors(j))
%     hold on
% end
% axis([0 tkhist(end) 0 1.1])
% legend('r1','r2')



