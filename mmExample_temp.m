% mmExample_temp.m
% Multiple-model estimation example: The solar panel deployment problem

clear;clc;
%----- Simulation parameters
dt = 0.1;
Nsim = 1000;
tkhist = [0:Nsim-1]'*dt;
% nmod is the number of models considered
nmod = 3;
nx = 2;
nz = 1;
nu = 1;

%----- Random number seed
n_seed = round(sum(clock*100));
randn('state', n_seed);

%----- Storage matrices
xtruekhist = zeros(Nsim,nx);
ztruekhist = zeros(Nsim,nz);
xhatkhist = zeros(Nsim,nx);
Pkhist = zeros(Nsim,nx);
% xhatkMhist is an array of all the mode-conditioned estimates
xhatkMhist = zeros(Nsim,nx,nmod);
PkMhist = zeros(Nsim,nx,nx,nmod);
mukhist = zeros(Nsim,nmod);
nukMhist = zeros(Nsim,nz,nmod);
LambdakMhist = zeros(Nsim,nmod);
FkM = zeros(nx,nx,nmod);
GkM = zeros(nx,nu,nmod);

%----- Set up multiple system models
cSet = [0.1;0.5;1];  hxSet = [1;2;3];
avec = cSet./hxSet;
bvec = 1./hxSet;
for j=1:nmod
    a = avec(j);
    b = bvec(j);
    % Discrete-time state transition matrix
    FkM(:,:,j) = ????
    % Discrete-time control matrix
    GkM(:,:,j) = ????
end

%----- Truth model
% switchHist gives the model switching time history and switchIndex gives
% the indices at and after which each model obtains.  For the static case,
% switchHist is either 1, 2, or 3, and switchIndex = 1;
switchHist = [1];
switchIndex = [1];
Qk = 0.001*diag([0.1 1]); Rq = chol(Qk);
vtruekhist = (Rq'*randn(nx,Nsim))';
Rkp1 = 0.1; Rr = chol(Rkp1);
wtruekhist = (Rr'*randn(nz,Nsim))';
utruekhist = 2*randn(Nsim,nu);
Hkp1 = [1 0];
x1 = [0;0.1];

%----- Generate truth-model states and measurements
xk = x1;
ii = 1;
for k=1:Nsim-1
    if k==switchIndex(ii)
        Fk = FkM(:,:,switchHist(ii));
        Gk = GkM(:,:,switchHist(ii));
        if ii<length(switchHist)
            ii = ii+1;
        end
    end
    kp1 = k + 1;
    uk = utruekhist(k,:)';
    vk = vtruekhist(k,:)';
    wkp1 = wtruekhist(kp1,:)';
    xkp1 = Fk*xk + Gk*uk + vk;
    zkp1 = Hkp1*xkp1 + wkp1;
    xtruekhist(k,:) = xk';
    ztruekhist(kp1,:) = zkp1';
    xk = xkp1;
end

%----- Run the Multiple-model filter
% Set up initial states and covariances
P1 = 10*eye(nx);
xhat1 = x1 + (chol(P1))'*randn(nx,1);
for j=1:nmod
    PkMhist(1,:,:,j) = P1;
    % Initialize the mode probabilities as equally probable
    mukhist(1,j) = ????
    xhatkMhist(1,:,j) = xhat1';
end
for k=1:Nsim-1
    kp1 = k + 1;
    zkp1 = ztruekhist(kp1,:)';
    uk = utruekhist(k,:)';
    for j=1:nmod
        % Propagation step
        Fkj = FkM(:,:,j);
        Gkj = GkM(:,:,j);
        xhatkj = xhatkMhist(k,:,j)';
        Pkj = squeeze(PkMhist(k,:,:,j));
        xbarkp1j = ????
        Pbarkp1j = ????
        
        % Measurement update
        nukp1j = zkp1 - Hkp1*xbarkp1j;
        Skp1j = ????
        invSkp1j = inv(Skp1j);
        Wkp1j = ????
        xhatkp1j = ????
        Pkp1j = ????
        
        % Calculate the likelihood of the current innovation
        Lambdakp1j = ????
        
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
    mukp1vec = ????
    mukhist(kp1,:) = mukp1vec';
    
    % Calculate the combined state estimate and covariance
    xMdum = reshape(squeeze(xhatkMhist(kp1,:,:)),nx,nmod);
    mudum = mukhist(kp1,:)';
    % Calculate the MMSE estimate averaged over all the modes
    xhatkp1 = (xMdum*mudum);
    % Calculate the estimation error covariance of xhatkp1
    Pkp1 = zeros(nx,nx);
    for j=1:nmod
        Pkp1j = squeeze(PkMhist(kp1,:,:,j));
        xhatkp1j = xhatkMhist(kp1,:,j)';
        mukp1j = mukhist(kp1,j);
        Pkp1 = ????
    end
    xhatkhist(kp1,:) = xhatkp1';
    Pkhist(kp1,:) = diag(Pkp1)';
end

%----- Display results
figure(1);clf;
plot(tkhist,mukhist);shg;
xlabel('Time (s)');
ylabel('Probability');
legend('\mu_1', '\mu_2','\mu_3');
title('Mode probability time histories');
ylim([0 1.1]);

figure(2);clf;
subplot(211)
iidum = 1:Nsim-1;
plot(tkhist(iidum), xtruekhist(iidum,1) - xhatkhist(iidum,1));
hold on;
plot(tkhist(iidum), sqrt(Pkhist(iidum,1)), 'k');
plot(tkhist(iidum), -sqrt(Pkhist(iidum,1)), 'k');
ylabel('\Delta \theta_z (rad)');
title('Estimation errors and covariances');

subplot(212)
iidum = 1:Nsim-1;
plot(tkhist(iidum), xtruekhist(iidum,2) - xhatkhist(iidum,2));
hold on;
plot(tkhist(iidum), sqrt(Pkhist(iidum,2)), 'k');
plot(tkhist(iidum), -sqrt(Pkhist(iidum,2)), 'k');
ylabel('\Delta \omega_z (rad/s)');
xlabel('Time (s)');



  
  


