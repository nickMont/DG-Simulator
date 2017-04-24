clear;clc;
%Implements the suboptimal solution of Kumar 1980 in which a player assumes
%that he is being spied upon

global nX_GG nU_GG
nX_GG=4; %size of full state space, even
nU_GG=nX_GG/2;

kumarConst;
P0=Q1f; %defined in kumarConst
N0=zeros(nX,nX);
M0=P0; %initial estimate covariance

i1=reshape(P0,[nX^2,1]);
i2=reshape(N0,[nX^2,1]);
i3=reshape(M0,[nX^2,1]);

tmax=25;
[t,y]=ode23s('rhodesOde',[0 tmax],[i1 i2 i3]);

pEvec=[];
nEvec=[];
mEvec=[];
fbM=[];
fbpos=[];
fbvel=[];
for i=1:length(t)
    pEvec(i,:)=y(i,1:nX^2);
    nEvec(i,:)=y(i,nX^2+1:2*nX^2);
    mEvec(i,:)=y(i,2*nX^2+1:3*nX^2);
    fbM(:,:,i)=-inv(R22)*G2'*reshape(pEvec(i,:),[nX,nX]);
    fbpos(i)=fbM(1,1,i);
    fbvel(i)=fbM(end,end,i);
end


clf;
figure(1)
plot(tmax-t,fbpos)
figure(2)
plot(tmax-t,fbpos)



