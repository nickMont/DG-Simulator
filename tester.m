
a=[1 2 3]
b=[5 6]
possibleCombosAtTimeKStore=[];
possibleCombosAtTimeK=combvec(a,b)
possibleCombosAtTimeKPrev=possibleCombosAtTimeK
for i=2:3
    possibleCombosAtTimeKPrev=combvec(possibleCombosAtTimeKPrev,possibleCombosAtTimeK)
end


