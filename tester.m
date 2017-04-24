% 
% a=[1 2 3]
% b=[5 6]
% possibleCombosAtTimeKStore=[];
% possibleCombosAtTimeK=combvec(a,b)
% possibleCombosAtTimeKPrev=possibleCombosAtTimeK
% for i=2:3
%     possibleCombosAtTimeKPrev=combvec(possibleCombosAtTimeKPrev,possibleCombosAtTimeK)
% end
% 
% 

% clf;
% load eICZero
% load eICOne
% mState=mean(eICZero,3);
% mState2=mean(eICOne,3);
% figure(1)
% subplot(2,1,1)
% plot(0:1:10,[5 mState(1,:)],'blue'); hold on;
% plot(0:1:10,[5 mState2(1,:)],'red');
% title('Monte Carlo simulation')
% legend('v_0=0','v_0=-0.1')
% xlabel('Time step')
% ylabel('Distance (km)')
% subplot(2,1,2)
% plot(0:1:10,[0 mState(2,:)],'blue'); hold on;
% plot(0:1:10,[-.1 mState2(2,:)],'red');
% xlabel('Time step')
% ylabel('Relative speed (km/T)')

x0=.1;
Qu=5;Ru=5;Qv=-5;Rv=1;
u1=1;u2=0; v1=1;v2=0;
Vpmat=10-[(x0+v1-u1)^2*Qu+u1^2*Ru (x0+v1-u1)^2*Qu+u1^2*Ru
    (x0+v1-u2)^2*Qu+u2^2*Ru (x0+v2-u2)^2*Qu+u2^2*Ru]
Vemat=10-[(x0+v1-u1)^2*Qv+v1^2*Rv (x0+v2-u1)^2*Qv+v2^2*Rv
    (x0+v1-u2)^2*Qv+v1^2*Rv (x0+v2-u2)^2*Qv+v2^2*Rv]
[coord,flag]=findNashEq(Vpmat,Vemat)











