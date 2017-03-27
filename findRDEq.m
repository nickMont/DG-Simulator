function [rdEq,flag,numEquivRDEQ] = findRDEq(A,B)
%takes in payoff matrices for row player A and column player B
%rdEq is the risk-dominant EQ point as a column.  If Y points are tied
%   for deviation penalty, this will return Y columns
%flag  =  2: RDEq found
%      =  1: one Nash Equilibrium
%      =  0: no Nash equilibrium
%      = -1: error
%numEquivRDEQ = number of points with identical minimax cost of deviation
%   (ie number of identical RDEQ points)

rdEq=[];
[m,n]=size(A);
maxL=max(m,n);
[neLoc,flag1]=findNashEq(A,B);
if flag1==0
    flag=0;
    numEquivRDEQ=0;
elseif flag1==-1
    flag=-1;
    numEquivRDEQ=0;
else  %if there are equilibria and there are no errors
    [~,numEq]=size(neLoc);
    if numEq==1 %Nash Eq is the only Eq
        rdEq=neLoc;
        flag=1;
        numEquivRDEQ=1;
    else
        flag=2;
        
        pairsToCheck=combnk(1:1:numEq,2);
        %returns two columns.  Reduce game to 2x2 payoffs defined by rows
        %and then compare RDEq there
        
        [L,~]=size(pairsToCheck);
        pairs=[];
        rdm=[];
        
        for i=1:L  %all possible point combinations
            p1=neLoc(:,pairsToCheck(i,1));
            p2=neLoc(:,pairsToCheck(i,2));
            pairs=[pairs  p1 p2]; %#ok
            tA=[A(p1(1),p1(2)) A(p1(1),p2(2))
                A(p2(1),p1(2)) A(p2(1),p2(2))];
            tB=[B(p1(1),p1(2)) B(p1(1),p2(2))
                B(p2(1),p1(2)) B(p2(1),p2(2))];
            %p1 is top left, p2 is bottom right, no matter how they were
            %originally arranged
            rd1=(tA(1,1)-tA(2,1))*(tB(1,1)-tB(1,2));
            rd2=(tA(2,2)-tA(1,2))*(tB(2,2)-tB(2,1));
            rdm=[rdm rd1 rd2]; %#ok
        end
        
        uniquePairs=(unique(pairs','rows'))';
        [~,numUnique]=size(uniquePairs);
        maxRD=zeros(1,numUnique);
        for i=1:numUnique
            indi = find(pairs(1,:)==uniquePairs(1,i) & pairs(2,:)==uniquePairs(2,i));
            rdThisPoint=rdm(indi); %#ok<FNDSB>
            maxRD(i)=max(rdThisPoint);
        end
        
        minValMaxRD=min(unique(maxRD));
        mindex=find(maxRD==minValMaxRD);
        numEquivRDEQ=length(mindex);
        rdEq=uniquePairs(:,mindex);
        
    end
end




end

