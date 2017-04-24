function [ J ] = J_pursuer(xu_hist)

% xu_hist is a n*(nX+nU)x1 vector arranged as [x1x2x3u1u2u3_t1
%   x1x2x3u1u2u3_t2 ...]
% time is determined by n
% can call this function multiple times to determine lowest n-cost

% Must handle rendezvous constraint in fmincon


[nX,nU,feva,~,~,dt]=get_problemSpecs;
sizeN=length(xu_hist);
n=sizeN/(nX+nU);

n_integer=isinteger(n); %check for integer time steps
if n_integer
    error('Non-integer number of time steps')
end

Rfull = 0 * diag(ones(n*nX/2,1));
Qfull = .15 * diag(ones(n*nU,1));
an=.1;

xhist=[];  %nX x n
uhist=[];  %nU x n
for i=1:n
    xhist=[xhist xu_hist((i-1)*(nX+nU)+1 : (i-1)*(nX+nU)+nX)];
    uhist=[uhist xu_hist((i-1)*(nX+nU)+nX+1 : (i-1)*(nX+nU)+nX+nU)];
end

posErrHist=[];
ureshape=[];
for i=1:n
    posErrHist=[posErrHist; xhist(1:nX/2,i)-feval(feva,i)];
    ureshape=[ureshape;uhist(:,i)];
end

Je=posErrHist'*Rfull*posErrHist;
Ju=ureshape'*Qfull*ureshape;
J=Je+Ju+an*n*dt;


end

