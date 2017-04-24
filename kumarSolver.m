clear;clc;
%Implements the suboptimal solution of Kumar 1980 in which a player assumes
%that he is being spied upon

global nX_GG
nX_GG=4; %size of full state space, even

kumarConst;
P10=P0;
P20=P0;
P30=P0;
Q1T=Q1f;
Q2T=Q2f;
Q3T=zeros(nX,nX);
Q4T=zeros(nX,nX);

i1=reshape(P10,[nX^2,1]);
i2=reshape(P20,[nX^2,1]);
i3=reshape(P30,[nX^2,1]);
i4=reshape(Q1T,[nX^2,1]);
i5=reshape(Q2T,[nX^2,1]);
i6=reshape(Q3T,[nX^2,1]);
i7=reshape(Q4T,[nX^2,1]);

tmax=25;
[t,y]=ode23s('kumarOde',[0 tmax],[i1 i2 i3 i4 i5 i6 i7]);

p1Evec=[];
p2Evec=[];
p3Evec=[];
q1Evec=[];
q2Evec=[];
q3Evec=[];
q4Evec=[];
fbM=[];
fbpos=[];
fbvel=[];
for i=1:length(t)
    p1Evec(i,:)=y(i,1:nX^2);
    p2Evec(i,:)=y(i,nX^2+1:2*nX^2);
    p3Evec(i,:)=y(i,2*nX^2+1:3*nX^2);
    q1Evec(i,:)=y(i,3*nX^2+1:4*nX^2);
    q2Evec(i,:)=y(i,4*nX^2+1:5*nX^2);
    q3Evec(i,:)=y(i,5*nX^2+1:6*nX^2);
    q4Evec(i,:)=y(i,6*nX^2+1:7*nX^2);
    fbM(:,:,i)=-inv(R22)*G2'*reshape(q2Evec(i,:),[nX,nX]);
    fbpos(i)=fbM(1,1,i);
    fbvel(i)=fbM(end,end,i);
end


clf;
figure(1)
plot(tmax-t,fbpos)
figure(2)
plot(tmax-t,fbpos)



