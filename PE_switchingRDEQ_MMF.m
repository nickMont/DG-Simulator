clear;clc;
%closes the loop for control switching by risk-dominance criterion

%NOTES
% Swap so as to use MC to generate OPTIMAL play and non-MC to generate ACTUAL play

tmax=10;
lookAheadTime=1; %number of steps to look ahead in optimization; set 0 to use t_remaining instead
dt=1;
MCmax=50;
nSim=ceil(tmax/dt); %simulation time within MC
%nSim=1;

flagDiminishHorizonNearNSim=1;  %if ==1, consider lookAheadTime
         %if TRemain>lookAheadTime, use TRemain if TRemain<lookAheadTime
flagUseDetm=0;  %zero noise if flag==1
flagUseFeedback=0; %generate control history via feedback
flagUseModeControlInsteadOfMean=1;

global destMatGG destNumGG QdestGG;
destNumGG=1;
destMatGG=[[1;5] [1;-1] [10;1]];
kTruth=1; %truth index
QdestGG=10*[eye(2) zeros(2,2); zeros(2,4)];
[~,numTargets]=size(destMatGG);

nsimTrunc=zeros(MCmax,1);
captureIndex=zeros(MCmax,1);
flagPlotAsGo=0;  %=1 to plot motion
flagPlotFinal=0; %=1 to plot final dist/speed history
plotint=.1;

nX=4;   %full size of state space (even number)
nU=nX/2;
nV=nX/2;
nZ=nX/2;

MC_rdeq_Max=100; %length of MC used to generate payoffs
if flagUseDetm==1
    MC_rdeq_Max=1;
end

mstay=.9; mswap=.1;
Mij=mstay*eye(numTargets)+(mswap/(numTargets-1))*(ones(numTargets,numTargets)-eye(numTargets,numTargets));

eStore=zeros(nX,nSim,MCmax);

mukMC=zeros(numTargets,nSim+1,MCmax);

for MCL=1:MCmax
    
    
    uphist=[];
    uehist=[];
    upDethist=[];
    ueDethist=[];
    
    C=300; %cost offset to guarantee positivity of V=C-J
    cdd=0;
    eyeHalfNX=eye(nX/2);
    zerosHalfNX=zeros(nX/2,nX/2);
    A_tr=[eyeHalfNX (1*dt-cdd*dt^2/2)*eyeHalfNX; zerosHalfNX (1-cdd*dt)*eyeHalfNX];
    B_tr=[dt^2/2*eyeHalfNX; dt*eyeHalfNX];
    Gammak=[dt^2/2*eyeHalfNX; dt*eyeHalfNX];
    H0=[eyeHalfNX zerosHalfNX];
    Hstack=[eyeHalfNX zeros(nX/2,nX/2) zeros(nX/2,nX) ; zeros(nX/2,nX) eyeHalfNX zeros(nX/2,nX/2)];
    P0=diag([.1*ones(1,nX/2) .1*ones(1,nX/2)]);
    P0e=P0; P0p=P0;
    Q0p=.01*eyeHalfNX; %process noise
    Q0e=Q0p;
    Q0s=Q0p*sqrt(2);
    R0p=.01*eye(nZ); %measurement noise
    R0e=R0p;
    Astack=[A_tr zeros(nX,nX); zeros(nX,nX) A_tr]; Bstack=[B_tr zeros(nX,nU); zeros(nX,nU) B_tr];
    Gammastack=[Gammak zeros(nX,nV); zeros(nX,nV) Gammak];
    Qstack=[Q0p zeros(nV,nV); zeros(nV,nV) Q0e]; Rstack=[R0e zeros(nZ,nZ);zeros(nZ,nZ) R0e];
    cholQ0P_T=chol(Q0p)';
    cholQ0E_T=chol(Q0e)';
    cholR0P_T=chol(R0p)';
    cholR0E_T=chol(Rstack)';
    P_pur=P0; P_eva=[P0p zeros(nX,nX); zeros(nX,nX) P0e]; %pursuer uses a stacked observation matrix
    PpPur=P0;
    PePur=P0;
    
    %cost matrices
    QfPur=8*[eyeHalfNX zerosHalfNX; zerosHalfNX zerosHalfNX];
    QstepPur=zeros(nX,nX);
    RPur=110*eye(nU);  %NOTE: discrete R = continuous R*dt^2 %it IS dt^2, NOT dt^2/2
    RPurScaled=RPur*dt^2;
    
    QfEva=5*[eyeHalfNX zerosHalfNX; zerosHalfNX zerosHalfNX];
    QstepEva=zeros(nX,nX);
    REva=125*eye(nU);
    REvaScaled=REva*dt^2;
    global QfEvaGG QfPurGG QstepEvaGG QstepPurGG REvaGG ...
        RPurGG umaxPurGG umaxEvaGG nXGG QnoiseGG dtGG cddGG nUGG%#ok
    QfEvaGG=QfEva; QfPurGG=QfPur; QstepEvaGG=QstepEva; QnoiseGG=Q0p; cddGG=cdd;
    QstepPurGG=QstepPur; REvaGG=REvaScaled; RPurGG=RPurScaled; nXGG=nX; dtGG=dt; nUGG=nU;
    muPrev_Eva=1/numTargets*ones(numTargets,1);

    flagBreakOnFlip=0;  %stop simulation on (estimated) collision if flag==1
    flagUseXPrevInsteadOfXModel=0;  %use overall xhat estimate instead of model-
    %specific xhat estimate in MMKF if flag==1
    %Best performance at flag==0.    
    
    mu_min=10^-2;
    normpdf_diag_DEBUG=.01; %with high confidence, diagonal elements of normpdf matrix go to zero
    
    if nX==2
        xEva=[0;0];
        xPur=[-1;.1];
    else
        xEva=[0;0;0;0];
        xPur=[-.5;0;0;0];
    end
    umaxPur=.5;
    umaxEva=.5;
    umaxPurGG=umaxPur; umaxEvaGG=umaxEva;
    JpurRunning=0;  %FOR COMPARISON WITH KUMARV2 ONLY
    JevaRunning=0;
    
    captureThresh=.5;
    flagSign1Flip=0;
    flagSign2Flip=0;
    
    %set randn's seed RNG
    n_seed = round(sum(clock*100));
    randn('state', n_seed);
    
    ehat0=xEva-xPur+(chol(P0))'*randn(nX,1);
    if flagUseDetm==1
        ehat0=xEva-xPur; %deterministic behavior
    end
    ehat_prev_pur=ehat0;
    
    PPstore=zeros(nX,nX,numTargets,nSim+1);
    for jk=1:numTargets
        PPstore(:,:,jk,1)=P0;
    end
    
    muPrev_pur=1/numTargets*ones(numTargets,1);
    mukhist=muPrev_pur;    
    mu_min=0.001;
    
    %----- Simulation parameters
    tkhist = [0:nSim]'*dt;
    
    %----- Simulation parameters
    tkhist = [0:nSim]'*dt;
    
    global umax_GG %#ok<TLEV,NUSED>
    
    LambdaTemp=zeros(numTargets,1);
    
    numSimRuns=0;
    i=1;
    
    if flagPlotAsGo==1
        figure(42);clf;
        gameStatePlot=scatter([xEva(1) xPur(1)],[xEva(2) xPur(2)],'red');
    end
    
    for i=1:nSim
        numSimRuns=numSimRuns+1;
        i %#ok<NOPTS>
        
        if i==1
            xEpur=xEva; xEeva=xEva; xPpur=xPur; xPeva=xPur;
            for jk=1:numTargets
                xEhathist(:,jk,1)=xEpur;
            end
        end
        
        
        %generate greedy open-loop controls
        if flagUseFeedback ~= 1
            if lookAheadTime==0
                ttf=ceil(tmax/dt + 1 - i);
            elseif flagDiminishHorizonNearNSim==1
                tRemain=nSim-i+1;
                if tRemain>=lookAheadTime
                    ttf=lookAheadTime;
                else
                    ttf=tRemain;
                end
            else
                ttf=lookAheadTime;
            end
            u1pX = -.05:.025:.25; %.177
            u1pY = -.25:.05:.25;
            u1eX = -.05:.025:.2; %.097
            u1eY = -.2:.05:.2;
            
            uconP=uConMat(ttf,umaxPur,u1pX,u1pY);
            uconE=uConMat(ttf,umaxEva,u1eX,u1eY);
        end
        
        uEvaExpectedMat=zeros(nU,numTargets); uPurTrue=zeros(nU,1);
        uPurExpected=zeros(nU,1); uEvaTrue=zeros(nU,1);
        for jk=1:numTargets
            
            destNumGG=jk;
            
            nmod_pur=length(uconP);
            nmod_eva=length(uconE);
            hP=ttf*ones(nmod_pur,1); %horizon length for pursuer, must define AFTER nmod calculation
            hE=ttf*ones(nmod_eva,1);
            KmatP=[];
            KmatE=[];
            typesE=ones(ttf,nmod_eva);
            typesP=ones(ttf,nmod_pur);
            
            uPurTrue=zeros(nU,1); uEvaExpected=uPurTrue;
            uclassVecPpur=zeros(nmod_pur,1); uclassVecEpur=zeros(nmod_eva,1);
            uclassVecPeva=zeros(nmod_pur,1); uclassVecEeva=zeros(nmod_eva,1);
            
            %counter variables for Truth/eXpected for mode control instead of
            %mean control
            uPTct=zeros(nmod_pur,1); uPXct=zeros(nmod_pur,1);
            uETct=zeros(nmod_eva,1); uEXct=zeros(nmod_eva,1);
            
            cctemp=1/MC_rdeq_Max; %normalization constant
            
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
            [JpPur,JePur]=generateCostMatrices_MMF(strategiesP,strategiesE,xPpur,xEpur);
            VpPur = -JpPur;
            VePur = -JePur;
            [eqLocP,nashReturnFlagP,~]=findRDEq(VpPur,VePur);
            
            if nashReturnFlagP>=1 %if there IS an RDEq
                uClassPp=eqLocP(1,1);
                uClassEp=eqLocP(2,1);
                [uPurTrueTemp,uEvaExpectedTemp] = processNashType1(eqLocP,umaxPur,umaxEva,uconP,uconE);
            elseif nashReturnFlagP==0 %suboptimal
                fprintf('Running LH2')
                %split for efficiency
                [uPurTrueTemp,uEvaExpectedTemp] = processNashType0(VpPur,VePur,umaxPur,umaxEva,uconP,uconE);
            end
            
            if jk==kTruth
                [JpEva,JeEva]=generateCostMatrices_MMF(strategiesP,strategiesE,xPeva,xEeva);
                VpEva = -JpEva;
                VeEva = -JeEva;
                [eqLocE,nashReturnFlagE,~]=findRDEq(VpEva,VeEva);
                if nashReturnFlagE>=1
                    [uPurExpected,uEvaTrue] = processNashType1(eqLocE,umaxPur,umaxEva,uconP,uconE);
                elseif nashReturnFlagE==0 %suboptimal
                    fprintf('Running LH2')
                    [uPurExpected,uEvaTrue] = processNashType0(VpEva,VeEva,umaxPur,umaxEva,uconP,uconE);
                end
            end
            
            uPurTrue=uPurTrue+uPurTrueTemp*muPrev_Eva(jk);
            uEvaExpectedMat(:,jk)=uEvaExpectedTemp;
        end
        uEt=uEvaTrue
        uEm=uEvaExpectedMat
        
        
        uphist=[uphist uPurTrue];
        uehist=[uehist uEvaTrue];
        
        %process noise
        if flagUseDetm==1 %deterministic behavior
            nP=zeros(nX,1);
            nE=zeros(nX,1);
        else
            nP=Gammak*cholQ0P_T*randn(nV,1);
            nE=Gammak*cholQ0E_T*randn(nV,1);
        end
        
        
        %Propagate states
        xPur=A_tr*xPur+B_tr*uPurTrue+nP;
        xEva=A_tr*xEva+B_tr*uEvaTrue+nE;
        eTrue=xEva-xPur;
        
        eStore(:,i,MCL)=eTrue;
        
        JpurRunning=JpurRunning+uPurTrue'*RPurScaled*uPurTrue;
        JevaRunning=JevaRunning+uEvaTrue'*REvaScaled*uEvaTrue;
        
        %Measurement
        zEva=Hstack*[xPur(1:nX);xEva(1:nX)]+cholR0E_T*randn(2*nZ,1);
        zPur=Hstack*[xPur(1:nX);xEva(1:nX)]+cholR0E_T*randn(2*nZ,1);
        if flagUseDetm==1
            zEva=Hstack*[xPur(1:nX);xEva(1:nX)];
            zPur=Hstack*[xPur(1:nX);xEva(1:nX)]; %deterministic behavior
        end
        
        %KF for evader
        ustackEva=[uPurExpected; uEvaTrue]; 
        xpurxevaEva=[xPeva;xEeva];
        [xpurxevaEva,P_eva]=linearKFStep(xpurxevaEva,zEva,Astack,Bstack,Gammastack,P_eva,Qstack,ustackEva,Hstack,Rstack);
        xPeva=xpurxevaEva(1:nX); xEeva=xpurxevaEva(nX+1:end);
        
        %KF for pursuer
        [xPpur,PpPur]=linearKFStep(xPpur,zPur(1:nZ),A_tr,B_tr,Gammak,PpPur,Q0p,uPurTrue,H0,R0p);
        
        %IMM
        lambdaTemp=zeros(numTargets,1);
        xhatkj=zeros(nX,numTargets);
        for jk=1:numTargets
            
            muI=muPrev_pur(jk);
            %Set up IMM mu values, entertaining the possibility of model
            %switching
            cc=zeros(numTargets,1);
            for i2=1:numTargets
                ccDum=0;
                for i3=1:numTargets
                    ccDum=ccDum+Mij(i3,jk)*muPrev_Eva(i3);
                end
                cc(i2)=ccDum;
            end
            muij=zeros(numTargets,1);
            for j2=1:numTargets
                muij(j2)=Mij(j2,jk)*1/cc(j2)*muPrev_Eva(j2);
            end
            xhatkj_beforeMotion=zeros(nX,1);
            Pk=zeros(nX,nX);
            for j4=1:numTargets
                xhatkj_beforeMotion=xhatkj_beforeMotion+xEhathist(:,jk,i)*muij(j4);
            end
            for j5=1:numTargets
                xDfTemp=xhatkj_beforeMotion-xEhathist(:,j5,i);
                Pk=Pk+muij(j5)*(PPstore(:,:,j5,i)+xDfTemp*xDfTemp');
            end
            zEpur=zPur(nZ+1:end);
            uE=uEvaExpectedMat(:,jk);
            xbar=A_tr*xhatkj_beforeMotion+B_tr*uE;
            Pbarj = A_tr*PePur*A_tr'+Gammak*Q0e*Gammak';
            
            nuj=zEpur-H0*xbar;
            Sk=H0*Pbarj*H0'+R0p;
            Sk_inv=inv(Sk);
            Wk=Pbarj*H0'*Sk_inv;
            xhatkj(:,jk)=xbar+Wk*nuj;
            Pj=Pbarj-Wk*Sk*Wk';
            xEhathist(:,jk,i+1)=xhatkj(:,jk);
            
            PPstore(:,:,jk,i+1)=Pj;
            
            muNu=zeros(nZ,1);
            normpdf_Eval = mvnpdf(nuj,muNu,Sk); %evaluates multivariate normal
            LambdaTemp(jk)=normpdf_Eval;
        end
        muStackTemp=zeros(numTargets,1);
        for jk=1:numTargets
            muStackTemp(jk)=LambdaTemp(jk)*muPrev_Eva(jk)/dot(LambdaTemp,muPrev_Eva);
        end
        for jk=1:numTargets
            if muStackTemp(jk)<=mu_min
                muStackTemp(jk)=mu_min;
            end
        end
        muPrev_Eva=muStackTemp/sum(muStackTemp);
        xEpur=zeros(nX,1);
        for jk=1:numTargets
            xEpur=xEpur+xhatkj(:,jk)*muPrev_Eva(jk);
        end
        mukhist=[mukhist muPrev_Eva];
        
        
        if flagPlotAsGo==1
            figure(42);delete(gameStatePlot)
            pause(plotint)
            gameStatePlot=scatter([xEva(1) xPur(1)],[xEva(2) xPur(2)],'red')
        end
        
        
    end
    nsimTrunc(MCL)=numSimRuns;
    
    mukMC(:,:,MCL)=mukhist;
    
    
end

muMean=mean(mukMC,3);

colors='rbgky';
figure(1);clf;
for jk=1:numTargets
    figure(1);
    hold on
    plot(dt*(0:nSim),muMean(jk,:),colors(jk))
end
xlabel('Time')
ylabel('Model probability \mu')
legend('Target 1','Target 2','Target 3')
title('Model probability for multiple target filter')

