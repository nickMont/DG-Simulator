function [Jpur,Jeva] = generateCostMatrices_DS(strategiesP,strategiesE,Jfall,xP0,xE0,numMC)
%inputs:
%     strategiesP (pursuer) is a struct containing
%         'matrices': nu x nx x T x nmodP  feedback matrices
%             where T is the maximum number of time steps considered by any
%             control
%             nu x nx zero matrix at (:,:,i) if types(i)==1
%         'constant': T x nmodP  constant max thrust vectors such that
%             u=uconst*unit(x)
%         'horizon': nmodP x 1  vector of integers containing finite-time horizon
%             length of cost function to evaluate, from 1 to T
%         'types':  T x nmodP vector of 0s,1s.  0 denotes top path, 1
%             denotes bottom path
%     strategiesE is the same for all evader strategies
%     xhat0 is the ehat of the player computing J matrices
%     noiseTypeFlag refers to how the cost due to noise is to be
%         calculated.
%             noiseTypeFlag=0 means ignore noise terms in J (sunk cost)
%             noiseTypeFlag=1 means to calculate them
%             This may also be done with an internal Monte Carlo but I prefer to be
%             efficient with my code.
% outputs:
%     Jpur, Jeva are cost matrices (NOT payoff matrices) for PE players
%         row player is pursuer, column player is evader
%         size nmodP x nmodE.
%         These matrices are E(J) (as the game is stochastic).  Thus, the
%         noise cost is higher for higher time horizons.

%notes:
%When calculating J, if the two time horizons for different strategies do
%not match, then J=9001^2;
%nU=nX/2 is assumed

thetamax=pi/6;

nmodP=length(strategiesP.horizon);
nmodE=length(strategiesE.horizon);

global QfEvaGG QfPurGG QstepEvaGG QstepPurGG REvaGG RPurGG umaxPurGG destGG QdestGG ...
   umaxEvaGG nXGG QnoisePGG QnoiseEGG dtGG nUGG boxBotTGG boxBotBGG boxTopBGG
QfEva=QfEvaGG; QfPur=QfPurGG; QstepEva=QstepEvaGG; dt=dtGG; 
QstepPur=QstepPurGG; REva=REvaGG; RPur=RPurGG; Qpdyn=QnoisePGG;
Qedyn=QnoiseEGG; umaxEva=umaxEvaGG; umaxPur=umaxPurGG; nX=nXGG;
nU=nUGG; boxBotT=boxBotTGG; boxBotB=boxBotBGG; boxTopB=boxTopBGG;
dest=destGG; Qdest=QdestGG;

tWidth=boxBotT-boxTopB;
bWidth=2*tWidth;

Jpur=zeros(nmodP,nmodE);
Jeva=zeros(nmodP,nmodE);

eyeNX=eye(nX);
zerosNX=zeros(nX,nX);

A_tr=[eyeNX];
B_tr=[dt*eyeNX];

TT=min(strategiesP.horizon(1),strategiesE.horizon(1)); %time check


for iiMC=1:numMC
    
    noisevecP=sqrt(Qpdyn)*randn(TT,1); noisevecE=sqrt(Qedyn)*randn(TT,1);
    
    for i=1:nmodP
        for j=1:nmodE
            xP=xP0; xE=xE0; Jpurloc=0; Jevaloc=0;
            
            TT=min(strategiesP.horizon(i),strategiesE.horizon(j)); %time check
            if TT ~= max(strategiesP.horizon(i),strategiesE.horizon(j)) %only compare strategies of equal length
                Jpur(i,j)=9001;
                Jeva(i,j)=9001;
            else
                TT=max(TT,1); %ensure that TT>=1;
                
                type=strategiesP.types(:,i);
                
                for k=1:TT
                    %use terminal vs. stepping cost
                    if k==TT
                        QuseP=QfPur;
                        QuseE=QfEva;
                    else
                        QuseP=QstepPur;
                        QuseE=QstepEva;
                    end
                    
                    if type==0
                        xTop=boxBotT;
                        C=tWidth;
                    else
                        xTop=boxBotB;
                        C=bWidth;
                    end
                    
                    
                    uconP=strategiesP.constant(:,:,:,i);
                    uconE=strategiesE.constant(:,:,:,j);
                    
                    yP=xTop-C/2-xP(2);
                    yE=xTop-C/2-xE(2);
                    
                    thetaP=controlAngle(yP,uconP,thetamax);
                    thetaE=controlAngle(yE,uconE,thetamax);
                    
                    vP=noisevecP(k); vE=noisevecE(k);
                    up=uconP*[cos(thetaP+vP);sin(thetaP+vP)];
                    ue=uconE*[cos(thetaE+vE);sin(thetaE+vE)];
                    
                    up=vectorSaturationF(up,0,umaxPur);
                    ue=vectorSaturationF(ue,0,umaxEva);
                    
                    %advance dynamics
                    xP=A_tr*xP+B_tr*up;
                    xE=A_tr*xE+B_tr*ue;
                    e=xP-xE;
                    
                    JfallLocP=0; JfallLocE=0;
                    if type==0
                        if xP(2)>boxBotT
                            JfallLocP=Jfall;
                        elseif xP(2)<boxTopB
                            JfallLocP=Jfall;
                        end
                        if xE(2)>boxBotT
                            JfallLocE=Jfall;
                        elseif xE(2)<boxTopB
                            JfallLocE=Jfall;
                        end
                    else
                        if xP(2)>boxBotB
                            JfallLocP=Jfall;
                        end
                        if xE(2)>boxBotE
                            JfallLocP=Jfall;
                        end
                    end
                    
                    %calculate costs
                    Jpurloc=Jpurloc + e'*QuseP*e + RPur*uconP^2 + JfallLocP;
                    Jevaloc=Jevaloc - e'*QuseE*e + REva*uconE^2 + JfallLocE;
                    if k==TT
                        Jevaloc = Jevaloc + (xE-dest)'*Qdest*(xE-dest);
                    end
                    
                end
                
                Jpur(i,j)=Jpur(i,j)+Jpurloc;
                Jeva(i,j)=Jeva(i,j)+Jevaloc;
            end
            
            
        end
    end
end






end

