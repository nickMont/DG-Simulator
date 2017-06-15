clear;clc;
%closes the loop for control switching by risk-dominance criterion

%NOTES
% Swap so as to use MC to generate OPTIMAL play and non-MC to generate ACTUAL play

tmax=10;
lookAheadTime=1; %number of steps to look ahead in optimization; set 0 to use t_remaining instead
dt=1;
MCmax=1;
nSim=ceil(tmax/dt); %simulation time within MC
nSim=1;

Jfall=9001;

global boxBotTGG boxBotBGG boxTopBGG destGG QdestGG;
boxBotTGG=.1; boxTopBGG=-.1; boxBotBGG=-.5; destGG=[1;0]; QdestGG=10*eye(2);
%--- boxBotT
%
%--- boxTopB
%---
%--- boxBotB
typeConst=0; %0 for top path, 1 for bottom path
Cbox=boxBotTGG-boxTopBGG;
if typeConst==0
    xTop=boxBotTGG;
else
    Cbox=2*Cbox;
    xTop=boxBotBGG;
end

flagDiminishHorizonNearNSim=1;  %if ==1, consider lookAheadTime
         %if TRemain>lookAheadTime, use TRemain if TRemain<lookAheadTime
flagUseDetm=0;  %zero noise if flag==1
flagUseFeedback=0; %generate control history via feedback
flagUseModeControlInsteadOfMean=1;

muAllP=zeros(MCmax,nSim,2);
muAllE=zeros(MCmax,nSim,2);
nsimTrunc=zeros(MCmax,1);
captureIndex=zeros(MCmax,1);
flagPlotAsGo=0;  %=1 to plot motion
flagPlotFinal=0; %=1 to plot final dist/speed history
plotint=.1;

nX=2;   %full size of state space (even number)
nU=1;
nV=nU;
nZ=2*nX;

eStore=zeros(nX,nSim,MCmax);

intMcNum=10; %length of MC inside of costmatrix generation

for MCL=1:MCmax
    
    
    uphist=[];
    uehist=[];
    
    Coff=300; %cost offset to guarantee positivity of V=C-J
    cd=0; %drag coefficient, =0 to ignore
    eyeNX=eye(nX);
    zerosNX=zeros(nX,nX);
    A_tr=eyeNX;
    B_tr=dt*eyeNX;
    thetaPpur=0; thetaEeva=0;
    Hstack=eye(nZ); %each player measures each player's positions
    P0=.1*eyeNX;
    P0stack=[P0 zerosNX; zerosNX P0];
    Q0p=.1; %process noise, pursuer
    Q0e=.005;
    Q0stack=[Q0p 0; 0 Q0e];
    R0P=.05*eye(nZ); %measurement noise
    R0E=R0P;
    cholQ0p_T=chol(Q0p)';
    cholQ0e_T=chol(Q0e)';
    cholR0P_T=chol(R0P)';
    cholR0E_T=chol(R0E)';
    
    %cost matrices
    QfPur=8*[eyeNX];
    QstepPur=zeros(nX,nX);
    RPur=110*eye(nU);  %NOTE: discrete R = continuous R*dt^2
    RPurScaled=RPur*dt^2;
    
    QfEva=5*[eyeNX];
    QstepEva=zeros(nX,nX);
    REva=125*eye(nU);
    REvaScaled=REva*dt^2;
    global QfEvaGG QfPurGG QstepEvaGG QstepPurGG REvaGG QnoisePGG QnoiseEGG...
        RPurGG umaxPurGG umaxEvaGG nXGG QnoiseGG dtGG cdGG nUGG%#ok
    QfEvaGG=QfEva; QfPurGG=QfPur; QstepEvaGG=QstepEva; QnoisePGG=Q0p; cdGG=cd; QnoiseEGG=Q0e;
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
        xEva=[0;0];
        xPur=[-1;0];
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
    ehat_prev_eva=ehat0;
    
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
            xPpur=xPur;
            xEpur=xEva;
            xPeva=xPur;
            xEeva=xEva;
        end
        x0p=[xPpur;xEpur]; x0e=[xPeva;xEeva]; %a priori estimates for xP,xE for EKF
        
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
            
            u1p=.05:.05:.25; %.177
            u1e=.05:.05:.25; %.097
            
            uconP=uConMat(ttf,umaxPur,u1p);
            uconE=uConMat(ttf,umaxEva,u1e);
        end
        
        
        nmod_pur=length(uconP);
        nmod_eva=length(uconE);
        hP=ttf*ones(nmod_pur,1); %horizon length for pursuer, must define AFTER nmod calculation
        hE=ttf*ones(nmod_eva,1);
        KmatP=[];
        KmatE=[];
        typesE=typeConst*ones(ttf,nmod_eva);
        typesP=typeConst*ones(ttf,nmod_pur);
        
        uclassVecPpur=zeros(nmod_pur,1); uclassVecEpur=zeros(nmod_eva,1);
        uclassVecPeva=zeros(nmod_pur,1); uclassVecEeva=zeros(nmod_eva,1);
        
        %counter variables for Truth/eXpected for mode control instead of
        %mean control
        uPTct=zeros(nmod_pur,1); uPXct=zeros(nmod_pur,1);
        uETct=zeros(nmod_eva,1); uEXct=zeros(nmod_eva,1);
        
        %if noise cost is flat across all time steps then ignore noise
        %find true optimal play on last iteration
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
        [JpPur,JePur]=generateCostMatrices_DS(strategiesP,strategiesE,Jfall,xPpur,xEpur,intMcNum);
        [JpEva,JeEva]=generateCostMatrices_DS(strategiesP,strategiesE,Jfall,xPeva,xEeva,intMcNum);
        
        %convert cost matrices to payoff matrices
        VpPur=Coff-JpPur;
        VePur=Coff-JePur;
        VpEva=Coff-JpEva;
        VeEva=Coff-JeEva;
        [eqLocP,nashReturnFlagP,~]=findRDEq(VpPur,VePur);
        [eqLocE,nashReturnFlagE,~]=findRDEq(VpEva,VeEva);
        
        %Process equilibria into controls
        if nashReturnFlagP>=1 %if there IS an RDEq
            uClassPp=eqLocP(1,1);
            uClassEp=eqLocP(2,1);
            [uPurTrueTemp,uEvaExpectedTemp] = processNashType1(eqLocP,umaxPur,umaxEva,uconP,uconE);
        elseif nashReturnFlagP==0 %suboptimal
            fprintf('Running LH2')
            %split for efficiency
            [uPurTrueTemp,uEvaExpectedTemp] = processNashType0(VpPur,VePur,umaxPur,umaxEva,uconP,uconE);
        end
        if nashReturnFlagE>=1
            [uPurExpectedTemp,uEvaTrueTemp] = processNashType1(eqLocE,umaxPur,umaxEva,uconP,uconE);
        elseif nashReturnFlagE==0 %suboptimal
            fprintf('Running LH2')
            [uPurExpectedTemp,uEvaTrueTemp] = processNashType0(VpEva,VeEva,umaxPur,umaxEva,uconP,uconE);
        end
        
        yP=xTop-Cbox/2-xPpur(2);
        yE=xTop-Cbox/2-xEpur(2);
        thetaPpur=real(controlAngle(yP,uPurTrueTemp,pi/6));
        thetadum=real(controlAngle(yE,uEvaExpectedTemp,pi/6));
        uPurTrue=uPurTrueTemp*[cos(thetaPpur);sin(thetaPpur)];
        uEvaExpected=uEvaExpectedTemp*[cos(thetadum);sin(thetadum)];
        
        yP=xTop-Cbox/2-xPeva(2);
        yE=xTop-Cbox/2-xEeva(2);
        thetadum=real(controlAngle(yP,uPurExpectedTemp,pi/6));
        thetaEeva=real(controlAngle(yE,uEvaTrueTemp,pi/6));
        uPurExpected=uPurExpectedTemp*[cos(thetadum);sin(thetadum)];
        uEvaTrue=uEvaTrueTemp*[cos(thetaEeva);sin(thetaEeva)];
        
        uphist=[uphist uPurTrue];
        uehist=[uehist uEvaTrue];
        
        %process noise
        if flagUseDetm==1 %deterministic behavior
            nP=0;
            nE=0;
        else
            nP=cholQ0p_T*randn;
            nE=cholQ0e_T*randn;
        end
        uPurTrueNoisy=uPurTrueTemp*[cos(thetaPpur+nP);sin(thetaPpur+nP)];
        uEvaTrueNoisy=uEvaTrueTemp*[cos(thetaEeva+nE);sin(thetaEeva+nE)];
        
        %Propagate states
        xPur=A_tr*xPur+B_tr*uPurTrueNoisy;
        xEva=A_tr*xEva+B_tr*uEvaTrueNoisy;
        
        %Measurement
        zEva=Hstack*[xPur;xEva]+cholR0E_T*randn(nZ,1);
        zPur=Hstack*[xPur;xEva]+cholR0E_T*randn(nZ,1);
        Astack=[eyeNX zerosNX; zerosNX eyeNX]; Bstack=[eyeNX zerosNX; zerosNX eyeNX];
        uCombPur=[uPurTrue;uEvaExpected];
        uCombEva=[uPurExpected;uEvaTrue];
        GammaLinStackP=[0 0; .01+uPurTrue(2) 0; 0 0; 0 .01+uEvaExpected(2)];
        GammaLinStackE=[0 0; .01+uPurExpected(2) 0; 0 0; 0 .01+uEvaTrue(2)];
        [xstackP,P_pur]=linearKFStep(x0p,zPur,Astack,Bstack,GammaLinStackP,P0stack,Q0stack,uCombPur,Hstack,R0P);
        [xstackE,P_eva]=linearKFStep(x0e,zEva,Astack,Bstack,GammaLinStackE,P0stack,Q0stack,uCombEva,Hstack,R0E);
        
        xPpur=xstackP(1:nX); xEpur=xstackP(nX+1:end); xPeva=xstackE(1:nX); xEeva=xstackE(nX+1:end);
        
        %Various outputs
%        kHist=KmatE(:,:,:,uClassEe)
%        kE=KmatE(:,:,1,uClassEe)
%        kP=KmatP(:,:,1,uClassPp)
        uP=uPurTrueNoisy
        uE=uEvaTrueNoisy
%        xPurLoc=xPur
%        xEvaLoc=xEva
        
        
        if flagPlotAsGo==1
            figure(42);delete(gameStatePlot)
            pause(plotint)
            gameStatePlot=scatter([xEva(1) xPur(1)],[xEva(2) xPur(2)],'red')
        end
        
        
    end
    nsimTrunc(MCL)=numSimRuns;
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
    plot(1:1:nSim,mState(1,:));
    title('Monte Carlo simulation')
    xlabel('Time step')
    ylabel('Distance')
    subplot(2,1,2)
    plot(1:1:nSim,mState(2,:));
    xlabel('Time step')
    ylabel('Relative speed')
end


figure(3);clf
subplot(2,1,1)
plot(dt*(1:1:numSimRuns),eStore(1,1:numSimRuns),'r')
hold on
plot(dt*(1:1:numSimRuns),eStore(2,1:numSimRuns),'b')
hLeg = legend('$\|e\|$','$\|\dot{e}\|$');
set(hLeg,'Interpreter','latex');
set(hLeg,'Interpreter','latex');
title('Relative pos and velocity')
subplot(2,1,2)
stairs(dt*(1:1:numSimRuns),uphist,'r')
hold on
stairs(dt*(1:1:numSimRuns),upDethist,'k')
hold on
stairs(dt*(1:1:numSimRuns),uehist,'b')
hold on
stairs(dt*(1:1:numSimRuns),ueDethist,'g')
hLeg = legend('$\|u_p\|$','$\|u_p\|_{opt}$','$\|u_e\|$','$\|u_e\|_{opt}$');
set(hLeg,'Interpreter','latex');
set(hLeg,'Interpreter','latex');
title('Applied control and optimal control')


% CODE GRAVEYARD
% generating/concatenating feedback controls
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



