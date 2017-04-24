clear;clc;clf;
%multiple model strategy determination, from the evader's perspective, of the pursuer

MCL=100;
nmod_tot=3;
muAll=zeros(30,nmod_tot,MCL);

for tMC=1:MCL
    
    nX=4;
    nU=2;
    nZ=2;
    dt=0.5;
    cd=0; %drag coefficient, =0 to ignore
    A_tr=[eye(nX/2) (1*dt-cd*dt^2/2)*eye(nX/2); zeros(nX/2,nX/2) (1-cd*dt)*eye(nX/2)];
    B_tr=[dt^2/2*eye(nX/2); dt*eye(nX/2)];
    Gammak=[dt^2/2*eye(nX/2); dt*eye(nX/2)];
    H=[eye(nX/2) zeros(nX/2)];
    P0=diag([.05*ones(1,nX/2) .1*ones(1,nX/2)]);
    Q0=.001*eye(nX/2);
    R0=.005*eye(nZ);
    
    nmod_pur=nmod_tot; %number of models considered
    q1_pur=5;
    q2_pur=[0  15 10000];
    r1_pur=1.5;
    ntrue_pur=2; %index of true model
    Mij = .25*ones(nmod_pur,nmod_pur); %probability of mode switching
    
    flagBreakOnFlip=0;  %stop simulation on (estimated) collision if flag==1
    flagUseDetm=0;  %zero/low noise if flag==1
    flagUseXPrevInsteadOfXModel=0;  %use overall xhat estimate instead of model-
    %specific xhat estimate in MMKF if flag==1
    %Best performance at flag==0.
    flagSwapAtTen=0;
    flagSwapAtTwenty=1;  %swap time index
    
    nSim=30;
    mu_min=10^-2;
    normpdf_diag_DEBUG=.01; %with high confidence, diagonal elements of normpdf matrix go to zero
    
    xEva=[20;10;0;0];
    xPur=[0;0;0;0];
    umaxPur=2.25;
    umaxEva=2;
    
    captureThresh=.5;
    flagSign1Flip=0;
    flagSign2Flip=0;
    
    %set randn's seed RNG
    n_seed = round(sum(clock*100));
    randn('state', n_seed);
    
    xhat0=xEva-xPur+(chol(P0))'*randn(nX,1);
    if flagUseDetm==1
        xhat0=xEva-xPur; %deterministic behavior
    end
    xhat_prev_pur=xhat0;
    
    Pstore_pur=zeros(nX,nX,nmod_pur,nSim);
    Pstore_pur(:,:,1,1)=P0;
    Pstore_pur(:,:,2,1)=P0;
    xhathist_pur=zeros(nX,nmod_pur,nSim);
    
    mukhist_pur=[];
    mu0=1/nmod_pur*ones(nmod_pur,1);
    LambdaVec_pur=zeros(nmod_pur,nSim);
    mukhist_pur=[mukhist_pur mu0];
    muPrev_pur=mu0;
    xhathist_weighted_pur=xhat0;
    
    %----- Simulation parameters
    tkhist = [0:nSim]'*dt;
    
    %----- Random number seed
    n_seed = round(sum(clock*100));
    randn('state', n_seed);
    
    global umax_GG
    umax_GG=umaxPur;
    
    e_true=zeros(nX,nSim);
    
    LambdaTemp=zeros(nmod_pur,1);
    
    optionsForFMC=optimset('Display','off');
    numSimRuns=1;
    for i=1:nSim
        numSimRuns=numSimRuns+1;
        i
        
        if i==10 && flagSwapAtTen==1
            if ntrue_pur==1
                ntrue_pur=2;
            else
                ntrue_pur=1;
            end
        end
        if i==20 && flagSwapAtTwenty==1
            if ntrue_pur==1
                ntrue_pur=2;
            else
                ntrue_pur=1;
            end
        end
        
        if i==1
            xhat=xhat0;
            for j=1:nmod_pur
                xhathist_pur(:,j,1)=xhat0;
            end
        else
            xhat=xhat_prev_pur;
        end
        
        %Evader control
        nDiff=xEva(1:nX/2)-xPur(1:nX/2);
        uEvaTrue=umaxEva*nDiff/norm(nDiff); %run away as quickly as possible
        %uEva=[0;0];
        
        %Pursuer control
        umax_GG=umaxPur;
        q2j=q2_pur(ntrue_pur); %index of true control QQ matrix
        QQ=[q1_pur*eye(nX/2) zeros(nX/2,nX/2); zeros(nX/2,nX/2) q2j*eye(nX/2)];
        RR=r1_pur*eye(nU);
        JJ=@(u) (A_tr*xhat-B_tr*u)'*QQ*(A_tr*xhat-B_tr*u) + u'*RR*u;
        [uPurTrue,Jc]=fmincon(JJ,zeros(nU,1),[],[], [],[], [],[], 'oneStepMaxT_constraint',optionsForFMC);
        
        xPur=A_tr*xPur+B_tr*uPurTrue;
        xEva=A_tr*xEva+B_tr*uEvaTrue;
        
        z=H*(xEva-xPur)+(chol(R0))'*randn(nZ,1);
        if flagUseDetm==1
            z=H*(xEva-xPur); %deterministic behavior
        end
        
        LambdaTemp=zeros(nmod_pur,1);
        for j=1:nmod_pur
            
            %         xhatkj_beforeMotion=xhathist_pur(:,j,i);
            %         if flagUseXPrevInsteadOfXModel==1
            %             xhatkj_beforeMotion=xhat_prev_pur;
            %         end
            %         Pk=Pstore_pur(:,:,j,i);
            
            muI=muPrev_pur(j);
            
            %Set up IMM mu values, entertaining the possibility of model
            %switching
            cc=zeros(nmod_pur,1);
            for l=1:nmod_pur
                ccDum=0;
                for i2=1:nmod_pur
                    ccDum=ccDum+Mij(i2,j)*muPrev_pur(i2);
                end
                cc(l)=ccDum;
            end
            muij=zeros(nmod_pur,1);
            for j2=1:nmod_pur
                muij(j2)=Mij(j2,j)*1/cc(j2)*muPrev_pur(j2);
            end
            
            xhatkj_beforeMotion=zeros(nX,1);
            Pk=zeros(nX,nX);
            for l2=1:nmod_pur
                xhatkj_beforeMotion=xhatkj_beforeMotion+xhathist_pur(:,j,i)*muij(l2);
            end
            for i3=1:nmod_pur
                xDfTemp=xhatkj_beforeMotion-xhathist_pur(:,i3,i);
                Pk=Pk+muij(i3)*(Pstore_pur(:,:,i3,i)+xDfTemp*xDfTemp');
            end
            
            umax_GG=umaxPur;
            q2j=q2_pur(j); %q2(j) refers to the jth model for control
            QQ=[q1_pur*eye(nX/2) zeros(nX/2,nX/2); zeros(nX/2,nX/2) q2j*eye(nX/2)];
            RR=r1_pur*eye(nU);
            JJ=@(u) (A_tr*xhatkj_beforeMotion-B_tr*u)'*QQ*(A_tr*xhatkj_beforeMotion-B_tr*u) + u'*RR*u;
            [uPur,Jc]=fmincon(JJ,zeros(nU,1), [],[], [],[], [],[], 'oneStepMaxT_constraint',optionsForFMC);
            
            %propagate state differences
            Qk=Q0;
            R=R0;
            ebarj = A_tr*xhatkj_beforeMotion + B_tr*uEvaTrue - B_tr*uPur;
            Pbarj = A_tr*Pk*A_tr'+Gammak*Qk*Gammak';
            
            %KF
            nuj=z-H*ebarj;
            Sk=H*Pbarj*H'+R;
            Sk_inv=inv(Sk);
            Wk=Pbarj*H'*Sk_inv;
            
            xhatkj=ebarj+Wk*nuj;
            Pj=Pbarj-Wk*Sk*Wk';
            
            PPstore(:,:,j,i+1)=Pj;
            xhathist_pur(:,j,i+1)=xhatkj;
            
            normpdf_Eval = normpdf(nuj,0,Sk); %returns a diagonal matrix of P(x1|z),P(x2|z),...
            for kk=1:nZ
                if normpdf_Eval(kk,kk)<=normpdf_diag_DEBUG
                    normpdf_Eval(kk,kk)=normpdf_diag_DEBUG;
                end
            end
            Lambdaj=trace(normpdf_Eval);
            LambdaVec_pur(i,j)=Lambdaj;
            LambdaTemp(j)=Lambdaj;
        end
        
        muPrev_pur=mukhist_pur(:,i);
        muStackTemp=zeros(nmod_pur,1);
        for j=1:nmod_pur
            muStackTemp(j)=LambdaTemp(j)*muPrev_pur(j)/dot(LambdaTemp,muPrev_pur);
        end
        for j=1:nmod_pur
            if muStackTemp(j)<=mu_min
                muStackTemp(j)=mu_min;
            end
        end
        muPrev_pur=muStackTemp/sum(muStackTemp);
        mukhist_pur=[mukhist_pur muPrev_pur]; %#okAGROW
        
        xhat_weighted_pur=zeros(nX,1);
        for j=1:nmod_pur
            ttt=muPrev_pur(j)*xhathist_pur(:,j,i+1);
            xhat_weighted_pur=xhat_weighted_pur+muPrev_pur(j)*xhathist_pur(:,j,i+1);
        end
        xhathist_weighted_pur=[xhathist_weighted_pur xhat_weighted_pur];
        
        %If captured at a time step or thresholds have flipped
        if norm(xhat_weighted_pur(1:nX/2))<=captureThresh
            fprintf('Captured \n')
            break
        end
        if xhat_weighted_pur(1)*xhat_prev_pur(1)<=0
            flagSign1Flip=1;
        end
        if xhat_weighted_pur(2)*xhat_prev_pur(2)<=0
            flagSign1Flip=1;
        end
        if flagSign1Flip==1 && flagSign2Flip==1 && flagBreakOnFlip==1
            break
        end
        
        xhat_prev_pur=xhat_weighted_pur;
        
        
        
    end
    
%     figure(1);clf;
%     colors='brgybrgybrgy';
%     for j=1:nmod_pur
%         plot(tkhist(1:numSimRuns),mukhist_pur(j,1:numSimRuns),colors(j))
%         hold on
%     end
%     axis([0 tkhist(end) 0 1.1])
%     legend('Collision','Rendezvous')
    
    muAll(1:numSimRuns,:,tMC)=mukhist_pur(:,1:numSimRuns)';


end

muMean=(mean(muAll,3))';

figure(1);clf;
colors='brgybrgybrgy';
for j=1:nmod_pur
    plot(1:1:numSimRuns,muMean(j,1:numSimRuns),colors(j))
    hold on
end
axis([1 numSimRuns 0 1.1])
if nmod_pur==2
    legend('Collision','Rendezvous')
else
    legend('u1','u2','u3','u4')
end

