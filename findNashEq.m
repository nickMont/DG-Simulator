function [eqLoc,flag,locC,locR] = findNashEq(A,B)
%takes in payoff matrices, returns matrix of Nash Equilibria, where
% each column represents NE coordinates

%inputs: A,B mxn payoff matrices for players 1,2
%              where 1 is row and 2 is col
%outputs:  eqLoc: 2xY matrix of Y pure nash equilibria for the game
%                 where eqLoc(1,:) is row location and (2,:) is column
%                 location on (A,B)
%          flag:  =1 if NE found
%                 =0 if no NE found (eqLoc=[])
%                 =-1 if error
%          locC,locR: column/row best plays


eqLoc=[];
[m,n]=size(A);
[mc,nc]=size(B);
if m~=mc 
    flag=-1;
elseif n~=nc
    flag=-1;
else %if no size mismatch
    
    locR=[];
    %row player
    for i=1:n
        col=A(:,i);
        mx=max(col);
        maxdex=(find(col==mx))';
        newEq=[maxdex; i*ones(1,length(maxdex))];
        locR=[locR newEq]; %#ok
    end
    
    %column player
    locC=[];
    for i=1:m
        ro=B(i,:);
        mx=max(ro);
        maxdex=find(ro==mx);
        newEq=[i*ones(1,length(maxdex));maxdex];
        locC=[locC newEq]; %#ok
    end
    
    [~,potentialNE]=size(locC);
    
    for i=1:potentialNE
        ind1=find(locR(1,:)==locC(1,i));
        ind2=find(locR(2,:)==locC(2,i));
        matchInd=intersect(ind1,ind2);
        eqLoc=[eqLoc locR(:,matchInd)]; %#ok
    end
    
    if isempty(eqLoc)
        flag=0;
    else
        flag=1;
    end
    
end



end

