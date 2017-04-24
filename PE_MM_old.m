clear;clc;
%multiple model strategy determination

nX=4;
nU=2;
nZ=2;
dt=.5;
A_tr=[eye(nX/2) dt*eye(nX/2); zeros(nX/2,nX/2) eye(nX/2)];
B_tr=[dt^2/2*eye(nX/2); dt*eye(nX/2)];

xEva=[5;10;0;0];
xPur=[0;0;0;0];
umaxPur=3;
umaxEva=2;

H=[eye(nX/2) zeros(nX/2)];
P0=.2*eye(nX);
Q0=.2*eye(nX/2);
R0=.2*eye(nX/2);

nmod=2; %number of models considered
q1=1;
q2=[0 5];
r1=.5;
ntrue=1; %index of true model

nSim=50;

%copied MM code
mu_min_thresh=10^-15;

%----- Simulation parameters
Nsim = 10;
tkhist = [0:Nsim-1]'*dt;

%----- Random number seed
n_seed = round(sum(clock*100));
randn('state', n_seed);

global umax_GG
umax_GG=umaxPur;

%----- Storage matrices
xtruekhist = zeros(Nsim,nX)';
ztruekhist = zeros(Nsim,nZ)';
xhatkhist = zeros(Nsim,nX)';
Pkhist = zeros(Nsim,nX);
% xhatkMhist is an array of all the mode-conditioned estimates
xhatkMhist = zeros(Nsim,nX,nmod);
PkMhist = zeros(Nsim,nX,nX,nmod);
mukhist = zeros(Nsim,nmod);
nukMhist = zeros(Nsim,nZ,nmod);
LambdakMhist = zeros(Nsim,nmod);

Gammak=[eye(nX/2);eye(nX/2)];


Hkp1 = H;
Hk = H;
Rkp1=R0;
x1 = xPur-xEva;

%----- Generate truth-model states and measurements
xk = x1;
switchIndex = [1 40];
switchHist=[1 2];
ii=1;


%----- Run the Multiple-model filter
% Set up initial states and covariances
P1 = P0;
Qk=Q0;
xtruekhist(:,1)=x1;
xhat1 = x1 + (chol(P1))'*randn(nX,1);

for j=1:nmod
    PkMhist(1,:,:,j) = P1;
    % Initialize the mode probabilities as equally probable
    mukhist(1,j) = 1/nmod;
    xhatkMhist(1,:,j) = xhat1';
    xhatkp1=xhat1;
end
for k=1:Nsim-1
    kp1 = k + 1;

    xhat=xhatkp1;
    
    ukEva=0*ones(nU,1)/sqrt(2);
    
    q2j=q2(j);
    QQ=[q1*eye(nX/2) zeros(nX/2,nX/2); zeros(nX/2,nX/2) q2j*eye(nX/2)];
    RR=r1*eye(nU);
    JJ=@(u) (A_tr*xhat+B_tr*u)'*QQ*(A_tr*xhat+B_tr*u) + u'*RR*u;
    [ukPur,Jc]=fmincon(JJ,zeros(nU,1),[],[], [],[], [],[], 'oneStepMaxT_constraint');
    
    xtruekhist(:,k+1)=A_tr*xtruekhist(:,k) + B_tr*ukPur - B_tr*ukEva;
    zkp1=H*xtruekhist(:,k+1);
    
    for j=1:nmod
        
        uk=zeros(nU,1);
        xhatkj = xhatkMhist(k,:,j)';
        Pkj = squeeze(PkMhist(k,:,:,j));
        
        %determine control
        q2j=q2(j); %q2(j) refers to the jth model for control
        QQ=[q1*eye(nX/2) zeros(nX/2,nX/2); zeros(nX/2,nX/2) q2j*eye(nX/2)];
        RR=r1*eye(nU);
        JJ=@(u) (A_tr*xhatkj+B_tr*u)'*QQ*(A_tr*xhatkj+B_tr*u) + u'*RR*u;
        [ukPur,Jc]=fmincon(JJ,zeros(nU,1),[],[], [],[], [],[], 'oneStepMaxT_constraint');
 
        
        % Propagation step
        xbarkp1j = A_tr*xhatkj+B_tr*ukPur+B_tr*ukEva;
        Pbarkp1j = A_tr*Pkj*A_tr'+Gammak*Qk*Gammak';
        
        % Measurement update
        nukp1j = zkp1 - Hkp1*xbarkp1j;
        Skp1j = Hk*Pbarkp1j*Hk'+Rkp1;
        invSkp1j = inv(Skp1j);
        Wkp1j = Pbarkp1j*Hk'*inv(Skp1j);
        xhatkp1j = xbarkp1j + Wkp1j*nukp1j;
        Pkp1j = Pbarkp1j - Wkp1j*Skp1j*Wkp1j';
        
        % Calculate the likelihood of the current innovation
        normpdf_Eval = normpdf(nukp1j,0,Skp1j); %returns a diagonal matrix of P(x1|z),P(x2|z),...
        Lambdakp1j=trace(normpdf_Eval);
        
        % Store state estimate and covariance, innovation, and likelihood
        xhatkMhist(kp1,:,j) = xhatkp1j';
        PkMhist(kp1,:,:,j) = Pkp1j;
        nukMhist(kp1,:,j) = nukp1j';
        LambdakMhist(kp1,j) = Lambdakp1j;
    end
    
    % Update the mode probabilities
    mukvec = mukhist(k,:)';
    Lambdakp1vec = LambdakMhist(kp1,:)';
    % mukp1vec is the nmod-by-1 vector that contains the mode probabilities
    % for the nmod modes at time index kp1; it must satisfy sum(mukp1vec) = 1.
    mukp1vec = Lambdakp1vec.*mukvec/sum(Lambdakp1vec.*mukvec);
    for j=1:nmod
        if mukp1vec(j) <= mu_min_thresh
            mukp1vec(j) = mu_min_thresh;
        end
    end
    %renormalize after flooring
    mukhist(kp1,:) = mukp1vec'/sum(mukp1vec);
    
    % Calculate the combined state estimate and covariance
    xMdum = reshape(squeeze(xhatkMhist(kp1,:,:)),nX,nmod);
    mudum = mukhist(kp1,:)';
    % Calculate the MMSE estimate averaged over all the modes
    xhatkp1 = (xMdum*mudum);
    % Calculate the estimation error covariance of xhatkp1
    Pkp1 = zeros(nX,nX);
    for j=1:nmod
        Pkp1j = squeeze(PkMhist(kp1,:,:,j));
        xhatkp1j = xhatkMhist(kp1,:,j)';
        mukp1j = mukhist(kp1,j);
        Pkp1 = Pkp1 + mudum(j)*(Pkp1j + (xhatkp1j-xhatkp1)*(xhatkp1j-xhatkp1)');
    end
    xhatkhist(:,kp1) = xhatkp1;
    Pkhist(kp1,:) = diag(Pkp1)';
    
   
end

xhatkhist=xhatkhist';
xtruekhist=xtruekhist';

%----- Display results
figure(1);clf;
plot(tkhist,mukhist);shg;
xlabel('Time (s)');
ylabel('Probability');
legend('\mu_1', '\mu_2','\mu_3');
title('Mode probability time histories');
ylim([0 1.1]);