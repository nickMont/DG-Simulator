clear;clc;
%Plays the stubborn game in which the evader is playing suboptimally by
%repeating previous plays with probability Ps

%NOTES
%Need to handle evader's filter for pursuer state


Ps=.4;
Po=.5;

lookAheadTime=1;

nSim=20;

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

nmod_eva=12;
thetaEva=0:30:330; %Evader controls are maximum thrust in 15deg increments
C=300; %normalization constant

flagBreakOnFlip=0;  %stop simulation on (estimated) collision if flag==1
flagUseDetm=0;  %zero noise if flag==1
flagUseXPrevInsteadOfXModel=1;  %use overall xhat estimate instead of model-
    %specific xhat estimate in MMKF if flag==1
    %Best performance at flag==0.
    
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
%possibleCombos=combvec(1:1:numControls,1:1:numTargets);

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
    
    
    %generate greedy open-loop controls
    for indentOpenLoop=1:1
        ttf=lookAheadTime;
        thetaPur = 0:5:355;
        n0=length(thetaPur)^ttf;
        uconP=zeros(nU,1,ttf,n0);
        if nX==2
            if ttf==1
                uconP(1,1,1,:)=u1p;
            elseif ttf >= 2
                possibleCombPrev=combvec(u1p,u1p);
                if ttf >= 3
                    for i1=3:ttf
                        possibleCombPrev=combvec(possibleCombPrev,u1p);
                    end
                end
                uconP(1,1,:,:)=possibleCombPrev;
            end
        elseif nX>=4
            if nX==4
                uconP=zeros(nU,1,ttf,n0);
                comb1vec=umaxPur*[cosd(thetaPur);sind(thetaPur)];
            else
                printf('Only supports x,y dimension')
            end
            if ttf==1
                uconP(:,:,1,:)=comb1vec;
            elseif ttf>=2
                possibleCombPrev=comb1vec;
                for i1=2:ttf
                    possibleCombPrev=combvec(comb1vec,possibleCombPrev);
                end
                for i1=1:length(possibleCombPrev)
                    for i2=0:floor(length(possibleCombPrev(:,i1))/nU)-1
                        nblock=i2*nU+1:(i2+1)*nU;
                        uconP(:,1,i2+1,i1)=possibleCombPrev(nblock,i1);
                    end
                end
            end
        end
        
        
        n0=length(thetaEva)^ttf;
        uconE=zeros(nU,1,ttf,n0);
        if nX==2
            if ttf==1
                uconE(1,1,1,:)=u1e;
            elseif ttf >= 2
                possibleCombPrev=combvec(u1e,u1e);
                if ttf >= 3
                    for i1=3:ttf
                        possibleCombPrev=combvec(possibleCombPrev,u1e);
                    end
                end
                uconE(1,1,:,:)=possibleCombPrev;
            end
        elseif nX>=4
            if nX==4
                uconE=zeros(nU,1,ttf,n0);
                comb1vec=umaxEva*[cosd(thetaEva);sind(thetaEva)];
            else
                comb1vec=combvec(combvec(u1p,u2p),u3p);
            end
            if ttf==1
                uconE(:,:,1,:)=comb1vec;
            elseif ttf>=2
                possibleCombPrev=comb1vec;
                for i1=2:ttf
                    possibleCombPrev=combvec(comb1vec,possibleCombPrev);
                end
                for i1=1:length(possibleCombPrev)
                    for i2=0:floor(length(possibleCombPrev(:,i1))/nU)-1
                        nblock=i2*nU+1:(i2+1)*nU;
                        uconE(:,1,i2+1,i1)=possibleCombPrev(nblock,i1);
                    end
                end
            end
        end
    end
    
    nmod_pur=length(uconP);
    nmod_eva=length(uconE);
    hP=ttf*ones(nmod_pur,1); %horizon length for pursuer, must define AFTER nmod calculation
    hE=ttf*ones(nmod_eva,1);
    KmatP=[];
    KmatE=[];
    typesE=ones(ttf,nmod_eva);
    typesP=ones(ttf,nmod_pur);
    
    uPurTrue=zeros(nU,1); uEvaTrue=uPurTrue; uEvaExpected=uPurTrue; uPurExpected=uPurTrue;
    
    %if noise cost is flat across all time steps then ignore noise
    minTE=min(typesE); minTP=min(typesP); maxTE=max(typesE); maxTP=max(typesP);
    if max(hP)==min(hP) && max(hE)==min(hE) && max(maxTE)==min(minTE)...
            && max(maxTP)==min(minTP) && min(minTP)==min(minTE)
        noiseIsSunkCost=0; %if all strategies for all times are same type AND same horizon
    else
        noiseIsSunkCost=1;
    end
    
    strategiesE=struct('matrices',KmatE,'constant',uconE,'horizon',hE,'types',typesE);
    strategiesP=struct('matrices',KmatP,'constant',uconP,'horizon',hP,'types',typesP);
    %KmatE, KmatP are nu x nx x T x nmodE/nmodP feedback matrices,
    %where T is the maximum time horizon considered
    [JpPur,JePur]=generateCostMatrices(strategiesP,strategiesE,xhatP,noiseIsSunkCost);
    [JpEva,JeEva]=generateCostMatrices(strategiesP,strategiesE,xhatE,noiseIsSunkCost);
    
    
    %convert cost matrices to payoff matrices
    VpPur=C-JpPur;
    VePur=C-JePur;
    VpEva=C-JpEva;
    VeEva=C-JeEva;
    [eqLocP,nashReturnFlagP,~]=findRDEq(VpPur,VePur);
    [eqLocE,nashReturnFlagE,~]=findRDEq(VpEva,VeEva);
    
    %Generate optimal pursuer control from NE
    for indentProcessNE_P=1:1
        if nashReturnFlagP>=1 %if there IS an RDEq
            uClassPp=eqLocP(1,1);
            uClassEp=eqLocP(2,1);
            uPurTrueOpt = vectorSaturationF(uconP(:,:,1,uClassPp),0,umaxPur);
            uEvaExpectedOpt = vectorSaturationF(uconE(:,:,1,uClassEp),0,umaxEva);
        elseif nashReturnFlagP==0 %suboptimal
            fprintf('Running LH2')
            nashMixedP=LH2(VpPur,VePur);
            KpM=zeros(size(KmatP(:,:,1,1)));
            KeM=zeros(size(KmatE(:,:,1,1)));
            nashP=nashMixedP{1};
            nashE=nashMixedP{2};
            u0p=zeros(nU,1);
            u0e=zeros(nU,1);
            for kk=1:length(nashP)
                if nashP(kk) > 0
                    if typesP(kk)==0
                        u0p=u0p+nashP(kk)*KmatP(:,:,1,kk)*xhatP;
                    else
                        u0p=u0p+nashP(kk)*uconP(:,:,1,kk);
                    end
                end
            end
            for kk=1:length(nashE)
                if nashE(kk) > 0
                    if typesE(kk)==0
                        u0e=u0e+nashE(kk)*KmatE(:,:,1,kk)*xhatP;
                    else
                        u0e=u0e+nashE(kk)*uconE(:,:,1,kk);
                    end
                end
            end
            uPurTrueOpt=vectorSaturationF(u0p,0,umaxPur);
            uEvaExpectedOpt=vectorSaturationF(u0e,0,umaxEva);
        end
    end
    
    %Genderate optimal evader control from NE
    for indendProcessNE_E=1:1
        if nashReturnFlagE>=1
            uClassPe=eqLocE(1,1);
            uClassEe=eqLocE(2,1);
            uPurExpected = vectorSaturationF(uconP(:,:,1,uClassPe),0,umaxPur);
            uEvaTrueOpt = vectorSaturationF(uconE(:,:,1,uClassEe),0,umaxEva);
        elseif nashReturnFlagE==0 %suboptimal
            fprintf('Running LH2')
            nashMixedE=LH2(VpEva,VeEva);
            KpM=zeros(size(KmatP(:,:,1,1)));
            KeM=zeros(size(KmatE(:,:,1,1)));
            nashP=nashMixedE{1};
            nashE=nashMixedE{2};
            u0p=zeros(nU,1);
            u0e=zeros(nU,1);
            for kk=1:length(nashP)
                if nashP(kk) > 0
                    if typesP(kk)==0
                        u0p=u0p+nashP(kk)*KmatP(:,:,1,kk)*xhatE;
                    else
                        u0p=u0p+nashP(kk)*uconP(:,:,1,kk);
                    end
                end
            end
            for kk=1:length(nashE)
                if nashE(kk) > 0
                    if typesE(kk)==0
                        u0e=u0e+nashE(kk)*KmatE(:,:,1,kk)*xhatE;
                    else
                        u0e=u0e+nashE(kk)*uconE(:,:,1,kk);
                    end
                end
            end
            uPurExpectedOpt=vectorSaturationF(u0p,0,umaxPur);
            uEvaTrueOpt=vectorSaturationF(u0e,0,umaxEva);
        end
    end
    
    %HANDLE PROBABILITIES HERE
    if i==1
        Pprev=Po;
        Pm=(1-Pprev)/(nmod_eva-1);
        pprob=Pm*ones(nmod_eva,1);
        pprob(uClassEe)=Pprev;
        uChooseIndex=randp(pprob,1);
        uEvaTrue=uconE(:,1,1,uChooseIndex);
        uPrevIndex=uChooseIndex;
    else
        if uPrevIndex==uClassEe %if the previous play was also this-time
                %optimal, Pm is split among more possible controls
            Pprev=Ps+Po;
            Pm=(1-Pprev)/(nmod_eva-1);
            pprob=Pm*ones(nmod_eva,1);
            pprob(uPrevIndex)=Pprev;
            uChooseIndex=randp(pprob,1);
            uEvaTrue=uconE(:,1,1,uChooseIndex);
            uPrevIndex=uChooseIndex;
        else
            Pm=(1-Ps-Po)/(nmod_eva-2);
            pprob=Pm*ones(nmod_eva,1);
            pprob(uPrevIndex)=Ps;
            pprob(uClassEe)=Po;
            uChooseIndex=randp(pprob,1);
            uEvaTrue=uconE(:,1,1,uChooseIndex);
            uPrevIndex=uChooseIndex;
        end
        
    end
    MijEnn=(muPrev_eva'./pprob)';
    MijE=MijEnn./sum(MijEnn,2); %normalizing
    
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
        
%         ehatkj_beforeMotion=ehathist_eva(:,j,i);
%         if flagUseXPrevInsteadOfXModel==1
%             ehatkj_beforeMotion=ehat_prev_eva;
%         end
%         
%         Pk=Pstore_eva(:,:,j,i);
        
        uEva=umaxEva*[cosd(thetaEva(j));sind(thetaEva(j))];
        
        muI=muPrev_eva(j);
        
        %Set up IMM mu values, entertaining the possibility of model
        %switching
        cc=zeros(nmod_eva,1);
        for l=1:nmod_eva
            ccDum=0;
            for i2=1:nmod_eva
                ccDum=ccDum+MijE(j,i2)*muPrev_eva(i2);
            end
            cc(l)=ccDum;
        end
        muij=zeros(nmod_eva,1);
        for j2=1:nmod_eva
            muij(j2)=MijE(j,j2)*1/cc(j2)*muPrev_eva(j2);
        end
        
        ehatkj_beforeMotion=zeros(nX,1);
        Pk=zeros(nX,nX);
        for l2=1:nmod_eva
            ehatkj_beforeMotion=ehatkj_beforeMotion+ehathist_eva(:,j,i)*muij(l2);
        end
        for i3=1:nmod_eva
            xDfTemp=ehatkj_beforeMotion-ehathist_eva(:,i3,i);
            Pk=Pk+muij(i3)*(Pstore_eva(:,:,i3,i)+xDfTemp*xDfTemp');
        end
        
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

















