clear;clc;
%closes the loop for control switching by risk-dominance criterion

%NOTES
% approximate filtering

tmax=10;
lookAheadTime=1; %number of steps to look ahead in optimization; set 0 to use t_remaining instead
dt=1;
MCmax=1;
nSim=ceil(tmax/dt); %simulation time within MC
nSim=1;

flagDiminishHorizonNearNSim=1;  %if ==1, consider lookAheadTime
         %if TRemain>lookAheadTime, use TRemain if TRemain<lookAheadTime
flagUseDetm=1;  %zero noise if flag==1
flagUseFeedback=0; %generate control history via feedback

muAllP=zeros(MCmax,nSim,2);
muAllE=zeros(MCmax,nSim,2);
nsimTrunc=zeros(MCmax,1);
captureIndex=zeros(MCmax,1);
flagPlotAsGo=0;  %=1 to plot motion
flagPlotFinal=0; %=1 to plot final dist/speed history
plotint=.1;

nX=2;   %full size of state space (even number)
nU=nX/2;
nV=nX/2;
nZ=nX/2;

eStore=zeros(nX,nSim,MCmax);

for MCL=1:MCmax
    
    
    uphist=[];
    uehist=[];
    
    C=300; %cost offset to guarantee positivity of V=C-J
    cd=0; %drag coefficient, =0 to ignore
    eyeHalfNX=eye(nX/2);
    zerosHalfNX=zeros(nX/2,nX/2);
    A_tr=[eyeHalfNX (1*dt-cd*dt^2/2)*eyeHalfNX; zerosHalfNX (1-cd*dt)*eyeHalfNX];
    B_tr=[dt^2/2*eyeHalfNX; dt*eyeHalfNX];
    Gammak=[dt^2/2*eyeHalfNX; dt*eyeHalfNX];
    H=[eyeHalfNX zeros(nX/2)];
    P0=diag([.1*ones(1,nX/2) .1*ones(1,nX/2)]);
    Q0=.01*eyeHalfNX; %process noise
    R0P=.05*eye(nZ); %measurement noise
    R0E=R0P;
    cholQ0_T=chol(Q0)';
    cholR0P_T=chol(R0P)';
    cholR0E_T=chol(R0E)';
    
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
        RPurGG umaxPurGG umaxEvaGG nXGG QnoiseGG dtGG cdGG nUGG%#ok
    QfEvaGG=QfEva; QfPurGG=QfPur; QstepEvaGG=QstepEva; QnoiseGG=Q0; cdGG=cd;
    QstepPurGG=QstepPur; REvaGG=REvaScaled; RPurGG=RPurScaled; nXGG=nX; dtGG=dt; nUGG=nU;
    
    
    nmod_pur=2; %number of models considered
    MijE = 1/nmod_pur*ones(nmod_pur,nmod_pur); %probability of mode switching
    
    nmod_eva=2; %number of models considered
    MijP = 1/nmod_pur*ones(nmod_eva,nmod_eva); %probability of mode switching

    flagBreakOnFlip=0;  %stop simulation on (estimated) collision if flag==1
    flagUseXPrevInsteadOfXModel=0;  %use overall xhat estimate instead of model-
    %specific xhat estimate in MMKF if flag==1
    %Best performance at flag==0.    
    
    mu_min=10^-2;
    normpdf_diag_DEBUG=.01; %with high confidence, diagonal elements of normpdf matrix go to zero
    
    if nX==2
        xEva=[5.000;0];
        xPur=[0;.1];
    else
        xEva=[5.000;3;0;0];
        xPur=[0;0;0;.1];
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
    
    Pstore_pur=zeros(nX,nX,nmod_pur,nSim);
    Pstore_pur(:,:,1,1)=P0;
    Pstore_pur(:,:,2,1)=P0;
    xhathist_pur=zeros(nX,nmod_pur,nSim);
    
    mukhist_pur=[];
    mu0=1/nmod_pur*ones(nmod_pur,1);
    LambdaVec_pur=zeros(nmod_pur,nSim);
    mukhist_pur=[mukhist_pur mu0];
    muPrev_pur=mu0;
    xhathist_weighted_pur=ehat0;
    
    
    %set randn's seed RNG
    n_seed = round(sum(clock*100));
    randn('state', n_seed);
    
    ehat0=xEva-xPur+(chol(P0))'*randn(nX,1);
    if flagUseDetm==1
        ehat0=xEva-xPur; %deterministic behavior
    end
    xhat_prev_eva=ehat0;
    
    Pstore_eva=zeros(nX,nX,nmod_eva,nSim);
    Pstore_eva(:,:,1,1)=P0;
    Pstore_eva(:,:,2,1)=P0;
    xhathist_eva=zeros(nX,nmod_eva,nSim);
    
    mukhist_eva=[];
    mu0=1/nmod_eva*ones(nmod_eva,1);
    LambdaVec_eva=zeros(nmod_eva,nSim);
    mukhist_eva=[mukhist_eva mu0];
    muPrev_eva=mu0;
    xhathist_weighted_eva=ehat0;
    
    %----- Simulation parameters
    tkhist = [0:nSim]'*dt;
    
    %----- Simulation parameters
    tkhist = [0:nSim]'*dt;
    
    %----- Random number seed
    n_seed = round(sum(clock*100));
    randn('state', n_seed);
    
    global umax_GG %#ok<TLEV,NUSED>
    
    e_true=zeros(nX,nSim);
    
    LambdaTemp=zeros(nmod_pur,1);
    
    optionsForFMC=optimset('Display','off');
    numSimRuns=1;
    i=1;
    
    if flagPlotAsGo==1
        figure(42);clf;
        gameStatePlot=scatter([xEva(1) xPur(1)],[xEva(2) xPur(2)],'red');
    end
    
    for i=1:nSim
        numSimRuns=numSimRuns+1;
        i %#ok<NOPTS>
        
        if i==1
            xhatP=ehat0;
            for j=1:nmod_pur
                xhathist_pur(:,j,1)=ehat0;
            end
            xhatE=ehat0;
            for j=1:nmod_eva
                xhathist_eva(:,j,1)=ehat0;
            end
        else
            xhatP=ehat_prev_pur;
            xhatE=xhat_prev_eva;
        end
        
        %generate greedy closed-loop controls
        for indentFeedback=1:1
            if flagUseFeedback==1
                %generate greedy feedback matrices for matrix case
                ttf=ceil(tmax/dt - i);
                k1p = .05:.1:.25;
                k2p = [0 .1];
                possibleCombAtOneTimeStep=combvec(k1p,k2p);
                possibleCombosPrev=possibleCombAtOneTimeStep;
                for i0=1:ttf
                    possibleCombosPrev=combvec(possibleCombosPrev, possibleCombAtOneTimeStep);
                end
                nnp=length(possibleCombosPrev);
                KpmatList=zeros(nU,nX,ttf,nnp);
                for i1=1:nnp
                    for i2=1:ttf
                        %grab the i2'th block of possibleCombosPrev
                        if nX==2
                            nblock=nU*(i2-1);
                            KpmatList(:,:,i2,i1)=[possibleCombosPrev(nblock+1,i1)...
                                possibleCombosPrev(nblock+2,i1)];
                        elseif nX==4
                            nblock=nU*(i2-1);
                            k1Ind=possibleCombosPrev(nblock+1,i1);
                            k2Ind=possibleCombosPrev(nblock+2,i1);
                            KpmatList(:,:,i2,i1)=[k1Ind 0 k2Ind 0; 0 k1Ind 0 k2Ind];
                        end
                    end
                end
                k1e=k1p;
                k2e=k2p;
                possibleCombAtOneTimeStep=combvec(k1e,k2e);
                possibleCombosPrev=possibleCombAtOneTimeStep;
                for i0=1:ttf
                    possibleCombosPrev=combvec(possibleCombosPrev,possibleCombAtOneTimeStep);
                end
                nne=length(possibleCombosPrev);
                KematList=zeros(nU,nX,ttf,nne);
                for i1=1:nne
                    for i2=1:ttf
                        %grab the i2'th block of possibleCombosPrev
                        if nX==2
                            nblock=nU*(i2-1);
                            KematList(:,:,i2,i1)=[possibleCombosPrev(nblock+1,i1)...
                                possibleCombosPrev(nblock+2,i1)];
                        elseif nX==4
                            nblock=nU*(i2-1);
                            k1Ind=possibleCombosPrev(nblock+1,i1);
                            k2Ind=possibleCombosPrev(nblock+2,i1);
                            KematList(:,:,i2,i1)=[k1Ind 0 k2Ind 0; 0 k1Ind 0 k2Ind];
                        end
                    end
                end
                sizeKmat=size(KpmatList);
                nmod_pur=nnp;
                nmod_eva=nne;
                hP=ttf*ones(nmod_pur,1); %horizon length for pursuer, must define AFTER nmod calculation
                hE=ttf*ones(nmod_eva,1);
                KmatP=KpmatList;
                KmatE=KematList;
                typesE=zeros(ttf,nmod_eva);
                typesP=zeros(ttf,nmod_pur);
                uconE=[]; %negation is handled in cost generation code
                uconP=[];
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
            
            k1 = 1.0e-04 *[-0.4066   -0.0004]*ehat_prev_pur;
            k2 = 1.0e-04 *[-0.3578   -0.0003]*ehat_prev_pur;
            
            u1p=[(k1-.0005):.0001:(k1+.005) (k2-.0005):.0001:(k2+.005) .01:.001:.2]; %.177
            u1e=u1p; %.097
            
            uconP=uConMat(ttf,umaxPur,u1p);
            uconE=uConMat(ttf,umaxEva,u1e);
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
        if nashReturnFlagP>=1 %if there IS an RDEq
            uClassPp=eqLocP(1,1);
            uClassEp=eqLocP(2,1);
            if typesP(uClassPp)==1
                uPurTrue = vectorSaturationF(uconP(:,:,1,uClassPp),0,umaxPur);
            else
                uPurTrue = vectorSaturationF(KmatP(:,:,1,uClassPp)*xhatP,0,umaxPur);
            end
            if typesE(uClassEp)==1
                uEvaExpected = vectorSaturationF(uconE(:,:,1,uClassEp),0,umaxEva);
            else
                uEvaExpected = vectorSaturationF(KmatE(:,:,1,uClassEp)*xhatE,0,umaxEva);
            end
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
            uPurTrue=vectorSaturationF(u0p,0,umaxPur);
            uEvaExpected=vectorSaturationF(u0e,0,umaxEva);
        end
        
        if nashReturnFlagE>=1
            uClassPe=eqLocE(1,1);
            uClassEe=eqLocE(2,1);
            if typesP(uClassPe)==1
                uPurExpected = vectorSaturationF(uconP(:,:,1,uClassPe),0,umaxPur);
            else
                uPurExpected = vectorSaturationF(KmatP(:,:,1,uClassPe)*xhatE,0,umaxPur);
            end
            if typesE(uClassEe)==1
                uEvaTrue = vectorSaturationF(uconE(:,:,1,uClassEe),0,umaxEva);
            else
                uEvaTrue = vectorSaturationF(KmatE(:,:,1,uClassEe)*xhatE,0,umaxEva);
            end
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
            uPurExpected=vectorSaturationF(u0p,0,umaxPur);
            uEvaTrue=vectorSaturationF(u0e,0,umaxEva);

        end
        
        uphist=[uphist uPurTrue];
        uehist=[uehist uEvaTrue];
        
        %process noise
        if flagUseDetm==1 %deterministic behavior
            nP=zeros(nX,1);
            nE=zeros(nX,1);
        else
            nP=Gammak*cholQ0_T*randn(nV,1);
            nE=Gammak*cholQ0_T*randn(nV,1);
        end
        
        
        %Propagate states
        xPur=A_tr*xPur+B_tr*uPurTrue+nP;
        xEva=A_tr*xEva+B_tr*uEvaTrue+nE;
        eTrue=xEva-xPur;
        
        eStore(:,i,MCL)=eTrue;
        
        JpurRunning=JpurRunning+uPurTrue'*RPurScaled*uPurTrue;
        JevaRunning=JevaRunning+uEvaTrue'*REvaScaled*uEvaTrue;
        
        xhat_prev_eva=eTrue+0.5*chol(P0)*randn(nX,1);
        ehat_prev_pur=eTrue+0.5*chol(P0)*randn(nX,1);
        
        %Various outputs
%        kHist=KmatE(:,:,:,uClassEe)
%        kE=KmatE(:,:,1,uClassEe)
%        kP=KmatP(:,:,1,uClassPp)
        uP=uPurTrue
        uE=uEvaTrue
%        xPurLoc=xPur
%        xEvaLoc=xEva
        e=eTrue;
        vp=VpPur;
        ve=VePur;
        
        
        if flagPlotAsGo==1
            figure(42);delete(gameStatePlot)
            pause(plotint)
            gameStatePlot=scatter([xEva(1) xPur(1)],[xEva(2) xPur(2)],'red')
        end
        
        
        %Measurement
        ze=H*(xEva-xPur)+cholR0E_T*randn(nZ,1);
        zp=H*(xEva-xPur)+cholR0P_T*randn(nZ,1);
        if flagUseDetm==1
            ze=H*(xEva-xPur);
            zp=H*(xEva-xPur); %deterministic behavior
        end
        

        
        %capture times
        captureIndex(MCL)=i; %will be overwritten continuously until
        %the player is captured (the loop breaks at capture, after
        %incrementing captureIndex).
        
        %If captured at a time step or thresholds have flipped
        if norm(eTrue(1:nX/2))<=captureThresh
            fprintf('Captured \n')
            break
        end
        %If captured at a time step or thresholds have flipped
        if norm(eTrue(1:nX/2))<=captureThresh
            fprintf('Captured \n')
            break
        end
    end
    nsimTrunc(MCL)=numSimRuns;
    
    %for comparison with KumarV2 ONLY
    JpurFinal_minimized=JpurRunning+eTrue'*QfPur*eTrue %minimized
    JevaFinal_maximized=-JevaRunning+eTrue'*QfEva*eTrue %maximized
    
    
end
% 
% fastestEnd=min(captureIndex);
% 
% finalEndPlot=finalEnd-1;
% 
% if MCmax>1 %plotting rules differ for multiple vs. single simulation lengths
%     figure(1);clf;
%     colors='brgybrgybrgy';
%     for j=1:nmod_pur
%         plot(tkhist(1:finalEndPlot),meanMuCalcP(j,1:finalEndPlot),colors(j))
%         hold on
%     end
%     title('Pursuer strategy')
%     axis([0 tkhist(finalEndPlot) 0 1.1])
%     legend('q1','q2')
%     
%     figure(2);clf;
%     colors='brgybrgybrgy';
%     for j=1:nmod_eva
%         plot(tkhist(1:finalEndPlot),meanMuCalcE(j,1:finalEndPlot),colors(j))
%         hold on
%     end
%     title('Evader strategy')
%     axis([0 tkhist(finalEndPlot) 0 1.1])
%     legend('r1','r2')
% else
%     finalEndPlot=captureIndex(1)
%     figure(1);clf;
%     colors='brgybrgybrgy';
%     for j=1:nmod_pur
%         plot(tkhist(1:finalEndPlot),meanMuCalcP(j,1:finalEndPlot),colors(j))
%         hold on
%     end
%     title('Pursuer strategy')
%     axis([0 tkhist(finalEndPlot) 0 1.1])
%     legend('q1','q2')
%     
%     figure(2);clf;
%     colors='brgybrgybrgy';
%     for j=1:nmod_eva
%         plot(tkhist(1:finalEndPlot),meanMuCalcE(j,1:finalEndPlot),colors(j))
%         hold on
%     end
%     title('Evader strategy')
%     axis([0 tkhist(finalEndPlot) 0 1.1])
%     legend('r1','r2')
% end
if flagPlotFinal==1
    mState=mean(eStore,3);
    figure(2)
    subplot(2,1,1)
    plot(1:1:10,mState(1,:));
    title('Monte Carlo simulation')
    xlabel('Time step')
    ylabel('Distance')
    subplot(2,1,2)
    plot(1:1:10,mState(2,:));
    xlabel('Time step')
    ylabel('Relative speed')
end







