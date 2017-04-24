clear;clc;close all;
%Monte carlo DG main

nmin=6;
nmax=15; %try to find solutions with time steps of size nmin through nmax
monteCarloSize=5;
nX=2*2;
nW1=1; %noise values to apply at each epoch
nW2=nW1; %generate noise for entire remaining trajectory or single epoch?
nU=nX/2;
xSTD=.1;
vSTD=0;

flagPlot=1;
flagGoToTen=1;
flagFM_options=1;

xPV0=[10;0;0 ; 0;0;0];
%xPV0=[10;0];
xPV0=[10;10; 0;0];
dtset=1; %update period, seconds
%set blanks for fdyn_eva or fdyn_pur for unknown pursuer/evader dynamics
set_problemSpecs(nX,nU,'fdyn_eva','fdyn_pur',xPV0,dtset);

oneDiagNx=diag(ones(nX/2,1));
A_tr=[oneDiagNx dtset*oneDiagNx
    zeros(nX/2,nX/2) oneDiagNx];
B_tr=[.5*dtset^2*oneDiagNx
    dtset*oneDiagNx];

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
for nstep=nmin:1:nmax
    numtimes=numtimes+1;
    
    
    n = nstep %time steps to rendezvous
    
    xPV=xPV0;
    
    xhistReal=[];
    uhistReal=[];
    xvuhistReal=[];
    for i=1:n
        
        xvuIter=zeros(nX+nU,monteCarloSize);
        for j=1:monteCarloSize
            
            %generation of xu0, the initial estimate for full x- and u- state history
            xu0=[];
            xu0t=[];
            xPV0_est=xPV+[xSTD*randn(nX/2,1);vSTD*randn(nX/2,1)];
            v0=xPV0_est(1+nX/2:end);
            xprev=xPV0_est(1:nX/2);
            if n-i>=1 && flagGoToTen==1 %goes to [10,0,0] at n-1 step
                for L=i:n-2
                    xnew=xprev+dtset*v0;
                    vnew=v0;
                    xu0t=[xu0t;xnew;vnew;zeros(nU,1)];
                    xprev=xnew;
                    vprev=vnew;
                end
                x_inter=xPV0(1:nX/2); x_inter(1)=10;
                um2=2/dtset^2*(x_inter-(xprev+dtset*xPV0(nX/2+1:nX)));
                vm2=xPV0(nX/2+1:nX)+dtset*um2;
                evaFinalPos=feval('fdyn_eva',n);
                ufinal=2/dtset^2*(evaFinalPos-(xprev+dtset*vm2));
                noiseCurr=randn(nW1,nW2);
                set_W(noiseCurr);  %can create a full noise trajectory; currently not in use
                xu0=[xu0t;x_inter;vm2;um2;evaFinalPos;v0+ufinal*dtset;ufinal];
            else %wait (u=0) and go to x=xEva at last step
                for L=i:n-1
                    xnew=xprev+dtset*v0;
                    vnew=v0;
                    xu0t=[xu0t;xnew;vnew;zeros(nU,1)];
                    xprev=xnew;
                    vprev=vnew;
                end
                evaFinalPos=feval('fdyn_eva',n);
                ufinal=2/dtset^2*(evaFinalPos-(xprev+dtset*xPV0(nX/2+1:nX)));
                noiseCurr=randn(nW1,nW2);
                set_W(noiseCurr);  %can create a full noise trajectory; currently not in use
                xu0=[xu0t;evaFinalPos;v0+ufinal*dtset;ufinal];
                
            end
            
            %xPV0 cannot be initialized at the start as xPV0 changes with step size
            set_problemSpecs(nX,nU,'fdyn_eva','fdyn_pur',xPV0_est,dtset);
            
            %solve for control
            Arendezvous=[zeros(nX/2,(nX+nU)*(n-i)) diag(ones(nX/2,1)) zeros(nX/2,nX/2) zeros(nU,nU) ];
            
            
            %call solver
            if flagFM_options==1
                options = optimoptions(@fmincon,'MaxFunctionEvaluations',5000,'Display','off');
                [xvuhist,J]=fmincon('J_pursuer',xu0,[],[],Arendezvous,evaFinalPos,[],[],'dynamicsConstraint',options);
            else
                [xvuhist,J]=fmincon('J_pursuer',xu0,[],[],Arendezvous,evaFinalPos,[],[],'dynamicsConstraint');
            end
            noiseHist(:,:,i,j)=noiseCurr;
            
            xvuIter(:,j)=xvuhist(1:nX+nU);
            
        end
        
        if monteCarloSize>=2
            uchosen=(mean((xvuIter(nX+1:end,:))'))';
        else
            uchosen=xvuIter(nX+1:end,:);
        end
        
        %advance dynamics per selected u
        xPV=A_tr*xPV+B_tr*uchosen;
        
        xhistReal=[xhistReal xPV(1:nX/2)];
        uhistReal=[uhistReal uchosen];
        xvuhistReal=[xvuhistReal; xPV; uchosen];
    end
    Jtot=J_pursuer(xvuhistReal);
    xvuHistAllMC(:,1:n,numtimes)=reshape(xvuhistReal,[nX+nU,n]);
    
    Jhist(numtimes)=Jtot;
    
end


if flagPlot==1
    j=0;
    for i=nmin:nmax
        j=j+1;
        hold on
        figure(1)
        plot([xPV0(1) xvuHistAllMC(1,1:end,j)],[xPV0(2) xvuHistAllMC(2,1:end,j)])
        hold on
        plot(5+cos(0:.1:2*pi),5+sin(0:.1:2*pi));
        axis([-2 12 -2 12])
    end
end












