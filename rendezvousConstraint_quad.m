function violation = rendezvousConstraint_quad(uhist)
%takes in uhist and determines whether or not the dynamics are feasible for
%a rendezvous at xF

global xF_GG ewxv0_GG

n=length(uhist)/4;
ewxv=ewxv0_GG;
xF=xF_GG;
for i=1:n
    uloc=uhist((i-1)*4+1:(i-1)*4+4);
    ewxv=rendezvousConstraint_quad([ewxv;uloc]);
end
violation=xF-ewxv(7:9);


end

