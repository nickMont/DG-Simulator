clear;clc;
%Drunken spider pursuit
%sets bridge width and runs simulation to see percentage of solutions
%crossing bridge

%This script calls PE_switchingRDEQ_DS.  As such, comment out clear;clc;,
%Jfall=..., and boxTopBGG=...;boxBotTGG=...; to run this script

bridgeWidthVec=[0.001 0.0025 0.005 0.0075 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1];
%bridgeWidthVec=[0.01 0.1];
JfallVec=[0 5 100 9001];
%JfallVec=10;
M=length(JfallVec);
N=length(bridgeWidthVec);
global boxBotTGG boxTopBGG

percentPurTop=zeros(M,N); percentEvaTop=zeros(M,N);

for ik=1:M
    Jfall=JfallVec(ik);
    for jk=1:length(bridgeWidthVec)
        boxBotTGG=0+bridgeWidthVec(jk)/2;
        boxTopBGG=0-bridgeWidthVec(jk)/2;
        
        PE_switchingRDEQ_DS
        
        nnp=floor(length(playedPpur)/2);
        nne=floor(length(playedEpur)/2);
        percentPurTop(ik,jk)=sum(playedPpur(1:nnp))/sum(playedPpur);
        percentEvaTop(ik,jk)=sum(playedEpur(1:nne))/sum(playedEpur);
        
    end
end


colors='rgkby';
figure(1);clf;
for ik=1:M
    hold on
    plot(bridgeWidthVec,percentPurTop(ik,:),colors(ik))
    % hold on
    % plot(bridgeWidthVec,percentEvaTop,'g')
    % hold off
end
legend('J_{fall}=0','J_{fall}=5','J_{fall}=100','J_{fall}=9000')
xlabel('Width of bridge (m)')
ylabel('Percent of solutions taking bridge')
title('Path chosen to cross the lake')
axis([0 max(bridgeWidthVec) 0 1])



