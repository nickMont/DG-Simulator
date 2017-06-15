function [uP,uE] = processNashType1(eqLoc,umaxPur,umaxEva,uconP,uconE)
%IF nashReturnFlag==1
%split into two functions for efficiency
uClassPp=eqLoc(1,1);
uClassEp=eqLoc(2,1);

uP = vectorSaturationF(uconP(:,:,1,uClassPp),0,umaxPur);

uE = vectorSaturationF(uconE(:,:,1,uClassEp),0,umaxEva);


end

