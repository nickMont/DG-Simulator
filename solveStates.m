clear;clc;

initializeQuad

[ewxvuVec,violation]=fsolve('trajectoryConstraint',[0;0;0;0;0;0; 0;0;0;0;0;0; 5.5;5;5;5.5; .2;0;0;0;0;0; .5;0;0;.5;0;0]);
ewxvuDum=[ewxvuVec;0;0;0;0];
ewxvu16N=reshape(ewxvuDum,16,length(ewxvuDum)/16);
uHist=ewxvu16N(13:16,:);

violation

