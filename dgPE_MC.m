clear;clc;close all;
%Monte carlo DG main

% %NOTES
%x~pursuer
%y~evader

nmin=1;
nmax=3; %try to find solutions with time steps of size nmin through nmax
monteCarloSize=25;
nX=2*2;
nW1=nX; %noise values to apply at each epoch
nW2=1; %generate noise for entire remaining trajectory or single epoch?
nU=nX/2;
r0=.01; %sensor noise floor
r1=.05; %scale factor

flagRunMC=0; %1 to run algorithm, 0 to avoid
flagRunLDY=1;
flagPlot=1;
flagGoToTen=1;
flagFM_options=0;
flagUseMinimaxControl=0; %use worst-case control if 1; use mean control instead if 0

xPV0=[10;0;0 ; 0;0;0];
%xPV0=[10;0];
xPV0=[10;10; 0;0];
evaPV_xEst0=[0;0;.1;.1];
evaPV_xEst=evaPV_xEst0;
evaPV_xEstest0=[2/sqrt(2);2/sqrt(2);0;0];
xPV=xPV0;
dtset=.5; %update period, seconds
%set blanks for fdyn_eva or fdyn_pur for unknown pursuer/evader dynamics
set_problemSpecs(nX,nU,'fdyn_eva','fdyn_pur',xPV0,dtset);

oneDiagNx=diag(ones(nX/2,1));
A_tr=[oneDiagNx dtset*oneDiagNx
    zeros(nX/2,nX/2) oneDiagNx];
B_tr=[.5*dtset^2*oneDiagNx
    dtset*oneDiagNx];
QQ=zeros(nX,nX); %assume no dynamic uncertainty for now
H=[1 0 0 0
    0 1 0 0]; %position sensor
RR=@(d) r0*diag([ones(nX/2,1)])+r1*[abs(d(1)) 0;0 abs(d(2))];
PP0=1*[.5*eye(2) zeros(2,2); zeros(2,2) .2*eye(2)]; %initial uncertainty in evaPV_xEst
PP=PP0;

j=1;

if length(xPV0)~=nX
    fprintf('xPV0 does not have the right number of states\n');
end
if nmin>nmax
    fprintf('nmin>nmax')
end

xhistReal=[];
uhistReal=[];
xvuhistReal=[];
Jhist=zeros(nmax-nmin+1,1);
numtimes=0;
xvuHistAllMC=zeros(nX+nU,nmax,nmax-nmin+1);  %each row corresponds to a run of length nstep
numSimRuns=0;
cont=true;
xhistReal=[];
uhistReal=[];
xvuhistReal=[];
Jcprev=9001;



% %Pursuer motion

% if flagRunMC==1
%     evaPV_xEstest=evaPV_xEstest0;
%     %MC for optimal control without information gain step
%     while cont && numSimRuns<=8
%         numSimRuns=numSimRuns+1;
%         Jstore=[]; %[monteCarloSize x nmax-nmin]
%         ustore=[]; %block matrix of nmax-nmin blocks of size nU x monteCarloSize
%         %n=1 [uMC1 uMC2 uMC3
%         %n=2  uMC1 uMC2 uMC3];
%         
%         %Investigate possible step lengths of nmin to nmax
%         for nstep=nmin:1:nmax
%             numtimes=numtimes+1;
%             
%             n = nstep; %time steps to rendezvous
%             
%             xvuIter=zeros(nX+nU,monteCarloSize);
%             Jrow=[];
%             urow=[];
%             
%             
%             for j=1:monteCarloSize
%                 
%                 %noisemat=[RR(abs(xPV(1:2)-[0;0])) zeros(2,2); zeros(2,4)]*randn(nW1,nW2);
%                 noisemat=PP*randn(nW1,n);
%                 set_W(noisemat);
%                 evaPV_xEste=evaPV_xEstest; %temporary evader PV estimate
%                 for L=1:n
%                     evaPV_xEste=A_tr*evaPV_xEste+noisemat(:,L);
%                 end
%                 finalPosEst=evaPV_xEste(1:2);
%                 finalPVEst=evaPV_xEste;
%                 
%                 %generation of xu0, the initial estimate for full x- and u- state history
%                 xu0=[];
%                 xu0t=[];
%                 xPV0_est=xPV;
%                 v0=xPV0_est(1+nX/2:end);
%                 xprev=xPV0_est(1:nX/2);
%                 for L=1:n
%                     xnew=xprev+dtset*v0;
%                     vnew=v0;
%                     xu0t=[xu0t;xnew;vnew;zeros(nU,1)];
%                     xprev=xnew;
%                     vprev=vnew;
%                 end
%                 xu0=xu0t;
%                 
%                 %xPV0 cannot be initialized at the start as xPV0 changes with step size
%                 set_problemSpecs(nX,nU,'fdyn_eva','fdyn_pur',xPV0_est,dtset);
%                 
%                 Arendezvous=[zeros(nX/2,(nX+nU)*(n-1)) diag(ones(nX/2,1)) zeros(nX/2,nX/2) zeros(nX/2,nU) ];
%                 b=finalPosEst;
%                 
%                 %call solver
%                 if flagFM_options==1
%                     options = optimoptions(@fmincon,'MaxFunctionEvaluations',5000,'Display','off');
%                     [xvuhist,J]=fmincon('J_PEMC',xu0,[],[],Arendezvous,b,[],[],'dynamicsConstraint_PEMC',options);
%                     uC=xvuhist((n-1)*(nX+nU)+nX+1:end);
%                 else
%                     [xvuhist,J]=fmincon('J_PEMC',xu0,[],[],Arendezvous,b,[],[],'dynamicsConstraint_PEMC');
%                     uC=xvuhist((n-1)*(nX+nU)+nX+1:end);
%                 end
%                 
%                 urow=[urow uC];
%                 Jrow=[Jrow J];
%                 
%             end
%             
%             ustore=[ustore; urow];
%             Jstore=[Jstore;Jrow];
%             
%         end
%         
%         if flagUseMinimaxControl==1
%             [maxJ,maxdex]=max(Jstore');
%             [~,mindex]=min(maxJ);
%             uchosen=ustore( (mindex-1)*nU+1:(mindex-1)*nU+nU , maxdex(mindex));
%             Jc=Jstore(mindex,maxdex(mindex));
%         else
%             meanJ=mean(Jstore')';
%             umean=mean(ustore')';
%             [~,mindex]=min(meanJ);
%             uchosen=umean( (mindex-1)*nU+1 : (mindex-1)*nU+nU);
%             Jc=meanJ(mindex);
%         end
%         
%         
%         %advance dynamics per selected u
%         xPV=A_tr*xPV+B_tr*uchosen;
%         
%         xhistReal=[xhistReal xPV(1:nX/2)];
%         uhistReal=[uhistReal uchosen];
%         xvuhistReal=[xvuhistReal; xPV; uchosen];
%         
%         evaPV_xEst=A_tr*evaPV_xEst;
% 
%         %estimate ePV
%         ePVbar=A_tr*evaPV_xEstest;
%         PPapest=A_tr*PP*A_tr'+QQ;
%         y=evaPV_xEst(1:2)-H*ePVbar; %residuals; can add noise here for simulation
%         RR_est=RR(xPV(1:nX/2)-evaPV_xEstest(1:nX/2));
%         Sk=H*PPapest*H'+RR_est;
%         W=PPapest*H'*inv(Sk);
%         evaPV_xEstest=ePVbar+W*y;
%         PP=(eye(nX,nX)-W*H)*PPapest;
%         
%         if abs(Jcprev-Jc)<=.01 %if cost unjustified
%             cont=false;
%         end
%         if norm(xPV(1:nX/2)-evaPV_xEst(1:nX/2))<=.5 %if captured
%             cont=false;
%         end
%         
%         Jcprev=Jc;
%     end
% end
% 
% 
% if flagRunLDY==1
%     %motion planning for optimal control with information gain step
%     cont=true;
%     evaPV_xEst=evaPV_xEst0;
%     xPV=xPV0;
%     PP=PP0;
%     global PPmat_GG xest_GG QQ_GG ePV_GG r1_GG r0_GG
%     QQ_GG=QQ;
%     r1_GG=r1;
%     r0_GG=r0;
%     
%     evaPV_xEstest=evaPV_xEstest0;
%     numSimRuns=0;
%     Jcprev=9001;
%     xhistLDY=[];
%     uhistLDY=[];
%     xvuhistLDY=[];
%     while cont && numSimRuns<=10
%         numSimRuns=numSimRuns+1;
%         
%         PPmat_GG=PP;
%         xest_GG=xPV;
%         ePV_GG=evaPV_xEstest;
%         
%         lookAheadLength=1:5;
%         u0=zeros(nU,1);
%         
%         Jstore=[];
%         ustore=[];
%         
%         for i=lookAheadLength
%             u0=fsolve('rendezvous',zeros(nU*i,1)); %time series of [uT1; uT2; ...; uTF]
%             %call solver
%             if flagFM_options==1
%                 options = optimoptions(@fmincon,'MaxFunctionEvaluations',5000,'Display','off');
%                 [uc,Jc,exitFlag,~]=fmincon('J_PE_I_pur',u0,[],[],[],[],[],[],'rendezvousConstraint',options);
%             else
%                 [uc,Jc,exitFlag,~]=fmincon('J_PE_I_pur',u0,[],[],[],[],[],[],'rendezvousConstraint');
%             end
%             %         %Ideal u without rendezvous constraint
%             %         [uUN,JUN]=fminunc('J_PE_I',u0);
%             %         uUN1=uUN(1:nU)
%             uc1=uc(1:nU); %next-step control
%             if exitFlag>=0
%                 Jstore=[Jstore Jc]; %#ok<AGROW>
%                 ustore=[ustore uc1]; %#ok<AGROW>
%             else
%                 Jstore=[Jstore 9001^2];
%                 ustore=[ustore zeros(nU,1)];
%             end
%         end
%         
%         [minJ,mindex]=min(Jstore);
%         uchosen=ustore(:,mindex);
%         
%         %advance true pursuer
%         xPV=A_tr*xPV+B_tr*uchosen;
%         xhistLDY=[xhistLDY xPV(1:2)]; %#ok<AGROW>
%         uhistLDY=[uhistLDY uchosen]; %#ok<AGROW>
%         xvuhistLDY=[xvuhistLDY; xPV; uchosen]; %#ok<AGROW>
%         
%         %advance true ePV
%         evaPV_xEst=A_tr*evaPV_xEst;
%         
%         %estimate ePV
%         ePVbar=A_tr*evaPV_xEstest;
%         PPapest=A_tr*PP*A_tr'+QQ;
%         y=evaPV_xEst(1:2)-H*ePVbar; %residuals; can add noise here for simulation
%         RR_est=RR(xPV(1:nX/2)-evaPV_xEst(1:nX/2));
%         Sk=H*PPapest*H'+RR_est;
%         W=PPapest*H'*inv(Sk);
%         evaPV_xEstest=ePVbar+W*y;
%         PP=(eye(nX,nX)-W*H)*PPapest;
%         
%         if abs(Jcprev-Jc)<=.01 %if cost unjustified
%             cont=false;
%         end
%         if norm(xPV(1:nX/2)-evaPV_xEst(1:nX/2))<=.5 %if captured
%             cont=false;
%         end
%         
%     end
% end

% J_mc=J_PEMC(xvuhistReal)
% J_ldy=J_PEMC(xvuhistReal_ldY)

% xFullLDY=[xPV0(1:nX/2) xhistLDY];
% plot(xFullLDY(1,:), xFullLDY(2,:))
% axis([0 10 0 10])
% 


nmax=10;
yEva=[10;10;0;0];  %PV
xPurEst=[9;10;0;0];
for i=1:nmax
    
    %motion planning for optimal control with information gain step
    global PPmat_GG xEva_GG QQ_GG xPurEst_GG r1_GG r0_GG n0_GG nmax_GG uPurGuess_GG xest_GG ePV_GG;
    QQ_GG=QQ;
    r1_GG=r1;
    r0_GG=r0;
    
    nmax_GG=nmax;
    
    PPmat_GG=PP;
    xEva_GG=yEva;
    xPurEst_GG=xPurEst;
    n0_GG=i;
    
    xest_GG=yEva;
    ePV_GG=xPurEst; % swapped but needed to maintain symmetry in J_lDY
    
    u0=zeros((1+nmax-i)*nU,1);
    uPurGuess_GG=zeros((1+nmax-i)*nU,1);
    [uc,Jc,exitFlag,~]=fmincon('J_PE_I_eva',u0,[],[],[],[],[],[],'evasionConstraint');

    uc





end

















