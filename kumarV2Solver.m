clear;clc;
%Implements the suboptimal solution of Kumar 1980 in which a player assumes
%that he is being spied upon

clear tP PP1 PP2 PP3 QQ1 QQ2 QQ3 QQ4 tQ qqvec tPr ppvec

global nX_GG tP_GG PP1_GG PP2_GG PP3_GG tQ_GG QQ1_GG QQ2_GG QQ3_GG ...
    QQ4_GG tmax_GG tprev_GG uPstore_GG uEstore_GG flagRunV3 tQe flagNoSpies %#ok<NUSED>

flagAssumeDualSpies=1; %assume that BOTH players are spying if ==1

flagRunV3=0;
nX_GG=2; %size of full state space, even
nX=nX_GG;
tprev_GG=0;

kumarV2Const;

tmax_GG=2;
tmax=tmax_GG;

refineQsteps=2; %number of times to refine Q after an initial guess

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

%initial guess of M
tP=0:1:tmax;
PP1=[];
for i=1:length(tP) PP1(:,:,i)=P10; PP2(:,:,i)=P20; PP3(:,:,i)=P30; end

for i=1:refineQsteps %two iterations are sufficient, from my own testing
    tP_GG=tP;
    PP1_GG=PP1;
    PP2_GG=PP2;
    PP3_GG=PP3;
    
    [tQ,qqvec]=ode23s('kumarV2_Q',[0 tmax],[i4 i5 i6 i7]);
    
    QQ1=[]; QQ2=[]; QQ3=[]; QQ4=[];
    for j=1:length(tQ)
        QQ1(:,:,j)=reshape(qqvec(i,1:nX^2),[nX,nX]); %#ok
        QQ2(:,:,j)=reshape(qqvec(i,nX^2+1:2*nX^2),[nX,nX]); %#ok
        QQ3(:,:,j)=reshape(qqvec(i,2*nX^2+1:3*nX^2),[nX,nX]); %#ok
        QQ4(:,:,j)=reshape(qqvec(i,3*nX^2+1:4*nX^2),[nX,nX]); %#ok
    end
    tQ_GG=tmax-tQ;
    QQ1_GG=QQ1;
    QQ2_GG=QQ2;
    QQ3_GG=QQ3;
    QQ4_GG=QQ4;
    
    [tPr,ppvec]=ode23s('kumarV2_P',[0 tmax],[i1 i2 i3]);
    PP1=[]; PP2=[]; PP3=[];
    for j=1:length(tPr)
        PP1(:,:,j)=reshape(ppvec(i,1:nX^2),[nX,nX]);
        PP2(:,:,j)=reshape(ppvec(i,nX^2+1:2*nX^2),[nX,nX]);
        PP3(:,:,j)=reshape(ppvec(i,2*nX^2+1:3*nX^2),[nX,nX]);
    end
    tP=tPr; %the other ode goes backwards in time
    
    
end

tQe=tQ_GG;

global QQ2Eva QQ2Pur tQp
QQ2Eva=QQ2;

if flagAssumeDualSpies==1
    
    kumarV3Const;
    flagRunV3=1;
    
    tmax_GG=2;
    tmax=tmax_GG;
    
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
    
    %initial guess of M
    clear tP PP1 PP2 PP3 QQ1 QQ2 QQ3 QQ4 tQ qqvec tPr ppvec
    tP=0:1:tmax;
    PP1=[];
    for i=1:length(tP) PP1(:,:,i)=P10; PP2(:,:,i)=P20; PP3(:,:,i)=P30; end
    
    for i=1:refineQsteps %two iterations are sufficient, from my own testing
        tP_GG=tP;
        PP1_GG=PP1;
        PP2_GG=PP2;
        PP3_GG=PP3;
        
        [tQ,qqvec]=ode23s('kumarV2_Q',[0 tmax],[i4 i5 i6 i7]);
        
        QQ1=[]; QQ2=[]; QQ3=[]; QQ4=[];
        for j=1:length(tQ)
            QQ1(:,:,j)=reshape(qqvec(i,1:nX^2),[nX,nX]); %#ok
            QQ2(:,:,j)=reshape(qqvec(i,nX^2+1:2*nX^2),[nX,nX]); %#ok
            QQ3(:,:,j)=reshape(qqvec(i,2*nX^2+1:3*nX^2),[nX,nX]); %#ok
            QQ4(:,:,j)=reshape(qqvec(i,3*nX^2+1:4*nX^2),[nX,nX]); %#ok
        end
        tQ_GG=tmax-tQ;
        QQ1_GG=QQ1;
        QQ2_GG=QQ2;
        QQ3_GG=QQ3;
        QQ4_GG=QQ4;
        
        [tPr,ppvec]=ode23s('kumarV2_P',[0 tmax],[i1 i2 i3]);
        PP1=[]; PP2=[]; PP3=[];
        for j=1:length(tPr)
            PP1(:,:,j)=reshape(ppvec(i,1:nX^2),[nX,nX]);
            PP2(:,:,j)=reshape(ppvec(i,nX^2+1:2*nX^2),[nX,nX]);
            PP3(:,:,j)=reshape(ppvec(i,2*nX^2+1:3*nX^2),[nX,nX]);
        end
        tP=tPr; %the other ode goes backwards in time
        
        
    end
    
    QQ2Pur=QQ2; tQp=tQ_GG;
    
end

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
for i=1:length(tQ)
    q1Evec(i,:)=qqvec(i,0*nX^2+1:1*nX^2);
    q2Evec(i,:)=qqvec(i,1*nX^2+1:2*nX^2);
    q3Evec(i,:)=qqvec(i,2*nX^2+1:3*nX^2);
    q4Evec(i,:)=qqvec(i,3*nX^2+1:4*nX^2);
    fbM(:,:,i)=-inv(R22)*G2'*reshape(q2Evec(i,:),[nX,nX]);
    fbpos(i)=fbM(1,1,i);
    fbvel(i)=fbM(end,end,i);
end

global J1_GG J2_GG
J1_GG=0;
J2_GG=0;

xPur1=[0;.1];
xEva1=[5;0];
xTrue0=(xPur1-xEva1);
x1hat0=xTrue0; x2hat0=xTrue0;
X0=[xTrue0;x1hat0;x2hat0];

[t,X]=ode23s('kumarV2Game',[0 tmax],X0);
xEndTrue=X(end,2*nX+1:3*nX)
JpurToMin=J1_GG+xEndTrue*Q1f*xEndTrue'
JevaToMax=J2_GG+xEndTrue*Q2f*xEndTrue'
Jdiff=JpurToMin-JevaToMax


% clf;
% figure(1)
% plot(tmax-tQ,fbpos)



