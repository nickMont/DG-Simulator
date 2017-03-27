clear;clc;
%Implements the suboptimal solution of Kumar 1980 in which a player assumes
%that he is being spied upon

global nX_GG nU_GG MM_GG tM_GG tPN_GG PP_GG NN_GG
nX_GG=4; %size of full state space, even
nU_GG=nX_GG/2;

tmax=25;

kumarConst;
P0=Q1f; %defined in kumarConst
N0=zeros(nX,nX);
M0=P0; %initial estimate covariance

i1=reshape(P0,[nX^2,1]);
i2=reshape(N0,[nX^2,1]);
i3=reshape(M0,[nX^2,1]);

%initial guess of M
tM=0:1:tmax;
MM=[];
for i=1:length(tM) MM(:,:,i)=M0; end

for i=1:5
    tM_GG=tM;
    MM_GG=MM;
    
    [tPN,pnvec]=ode23s('rhodesV2Ode',[0 tmax],[i1 i2]);
    
    PP=[]; NN=[];
    for j=1:length(tPN)
        PP(:,:,j)=reshape(pnvec(i,1:nX^2),[nX,nX]); %#ok
        NN(:,:,j)=reshape(pnvec(i,nX^2+1:2*nX^2),[nX,nX]); %#ok
    end
    tPN_GG=tPN;
    PP_GG=PP;
    NN_GG=NN;
    
    [tMr,mmvec]=ode23s('rhodesV2_M',[0 tmax],i3);
    MM=[];
    for j=1:length(tMr)
        MM(:,:,j)=reshape(mmvec(i,:),[nX,nX]);
    end
    tM=tmax-tMr; %the other ode goes backwards in time


end

pEvec=[];
nEvec=[];
mEvec=[];
fbM=[];
fbpos=[];
fbvel=[];
for i=1:length(tPN)
    pEvec(i,:)=pnvec(i,1:nX^2);
    nEvec(i,:)=pnvec(i,nX^2+1:2*nX^2);
    fbM(:,:,i)=-inv(R22)*G2'*reshape(pEvec(i,:),[nX,nX]);
    fbpos(i)=fbM(1,1,i);
    fbvel(i)=fbM(end,end,i);
end


clf;
figure(1)
plot(tmax-tPN,fbpos)
figure(2)
plot(tmax-tPN,fbpos)



