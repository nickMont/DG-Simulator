function [uP,uE] = processNashType0(Vp,Ve,umaxPur,umaxEva,uconP,uconE)

%IF nashReturnFlag==0
%split into two functions for efficiency
nashMixedP=LH2(Vp,Ve);
nashP=nashMixedP{1};
nashE=nashMixedP{2};
u0p=zeros(nU,1);
u0e=zeros(nU,1);
for kk=1:length(nashP)
    if nashP(kk) > 0
        u0p=u0p+nashP(kk)*uconP(:,:,1,kk);
    end
end
for kk=1:length(nashE)
    if nashE(kk) > 0
        u0e=u0e+nashE(kk)*uconE(:,:,1,kk);
    end
end
uP=vectorSaturationF(u0p,0,umaxPur);
uE=vectorSaturationF(u0e,0,umaxEva);

end

