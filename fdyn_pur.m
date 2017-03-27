function xPur = fdyn_pur(n)

[nX,nU,~,~,xPV0,dt]=get_problemSpecs;

xPur=[0;0;1;1]+n*[.5;.5;0;0];


end

