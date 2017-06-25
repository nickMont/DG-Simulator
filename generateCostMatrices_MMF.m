function varargout = generateCostMatrices_MMF(strategiesP,strategiesE,xP0,xE0,answerFlag)
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
%     answerFlag:  ==0 -> return mean payoff
%                  ==1 -> return most-used play
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


nmodP=length(strategiesP.horizon);
nmodE=length(strategiesE.horizon);

if nargin==4
    answerFlag=0;
end

global QfEvaGG QfPurGG QstepEvaGG QstepPurGG REvaGG RPurGG umaxPurGG QdestGG ...
   umaxEvaGG nXGG QnoisePGG QnoiseEGG dtGG nUGG destNumGG destMatGG cddGG
QfEva=QfEvaGG; QfPur=QfPurGG; QstepEva=QstepEvaGG; dt=dtGG; 
QstepPur=QstepPurGG; REva=REvaGG; RPur=RPurGG; Qpdyn=QnoisePGG;
Qedyn=QnoiseEGG; umaxEva=umaxEvaGG; umaxPur=umaxPurGG; nX=nXGG;
nU=nUGG; cdd=cddGG;
Qdest=QdestGG;

dest=[destMatGG(:,destNumGG);0;0];


Jpur=zeros(nmodP,nmodE);
Jeva=zeros(nmodP,nmodE);

eyeHalfNX=eye(nX/2); zerosHalfNX=zeros(nX/2,nX/2);

A_tr=[eyeHalfNX (1*dt-cdd*dt^2/2)*eyeHalfNX; zerosHalfNX (1-cdd*dt)*eyeHalfNX];
B_tr=[dt^2/2*eyeHalfNX; dt*eyeHalfNX];

TT=min(strategiesP.horizon(1),strategiesE.horizon(1)); %time check

playedP=zeros(nmodP,1); playedE=zeros(nmodE,1);
for i=1:nmodP
    for j=1:nmodE
        xP=xP0; xE=xE0; Jpurloc=0; Jevaloc=0;
        
        TT=min(strategiesP.horizon(i),strategiesE.horizon(j)); %time check
        if TT ~= max(strategiesP.horizon(i),strategiesE.horizon(j)) %only compare strategies of equal length
            Jpur(i,j)=9001;
            Jeva(i,j)=9001;
        else
            TT=max(TT,1); %ensure that TT>=1;
            
            for k=1:TT
                %use terminal vs. stepping cost
                if k==TT
                    QuseP=QfPur;
                    QuseE=QfEva;
                else
                    QuseP=QstepPur;
                    QuseE=QstepEva;
                end
                
                uconP=strategiesP.constant(:,:,:,i);
                uconE=strategiesE.constant(:,:,:,j);
                
                up=vectorSaturationF(uconP,0,umaxPur);
                ue=vectorSaturationF(uconE,0,umaxEva);
                
                %advance dynamics
                xP=A_tr*xP+B_tr*up;
                xE=A_tr*xE+B_tr*ue;
                e=xP-xE;

                %calculate costs
                Jpurloc=Jpurloc + e'*QuseP*e + uconP'*RPur*uconP;
                Jevaloc=Jevaloc - e'*QuseE*e + uconE'*REva*uconE;
                if k==TT
                    Jevaloc = Jevaloc + (xE-dest)'*Qdest*(xE-dest); %cost not payoff
                end
                
            end
            
            Jpur(i,j)=Jpur(i,j)+Jpurloc;
            Jeva(i,j)=Jeva(i,j)+Jevaloc;
        end
        
        
    end
end



if answerFlag==0
    varargout{1}=Jpur; varargout{2}=Jeva;
elseif answerFlag==1
    varargout{1}=playedP; varargout{2}=playedE;
end




end

