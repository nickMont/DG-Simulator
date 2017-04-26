function uconP = uConMat(ttf,umaxPur,u1p,u2p,u3p)
%nargin sets nU
%uconP = nU x 1 x nMod x ttf  matrix of nMod possible control combinations


umaxPur=abs(umaxPur); %ensure that umaxPur>=0

nU=nargin-2;

u1L=length(u1p);
n0=length(u1p)^ttf;
uconP=zeros(nU,1,ttf,n0);
if nU==1
    for i=1:u1L
        u1p(i)=saturationF(u1p(i),-umaxPur,umaxPur);
    end
    if ttf==1
        uconP(1,1,1,:)=u1p;
    elseif ttf >= 2
        possibleCombPrev=combvec(u1p,u1p);
        if ttf >= 3
            for i1=3:ttf
                possibleCombPrev=combvec(possibleCombPrev,u1p);
            end
        end
        uconP(1,1,:,:)=possibleCombPrev;
    end
elseif nU>=2
    %generate possible u-combinations
    if nU==2
        uconP=zeros(u1L,1,ttf,n0^2);
        comb1vec=combvec(u1p,u2p);
        
    else %nU=3
        comb1vec=combvec(combvec(u1p,u2p),u3p);
    end
    
    %handle saturations
    for i=1:length(comb1vec)
        comb1vec(:,i)=vectorSaturationF(comb1vec(:,i),-umaxPur,umaxPur);
    end
    
    %generate stacked uconP matrix
    if ttf==1
        uconP(:,:,1,:)=comb1vec;
    elseif ttf>=2
        possibleCombPrev=comb1vec;
        for i1=2:ttf
            possibleCombPrev=combvec(comb1vec,possibleCombPrev);
        end
        for i1=1:length(possibleCombPrev)
            for i2=0:floor(length(possibleCombPrev(:,i1))/nU)-1
                nblock=i2*nU+1:(i2+1)*nU;
                uconP(:,1,i2+1,i1)=possibleCombPrev(nblock,i1);
            end
        end
    end
end


end

