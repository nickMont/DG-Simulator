function J = J_quad_softConstraint(uhist)

global ewxv0est_GG xF_GG
xF=xF_GG;
ewxv=ewxv0est_GG;

n=length(uhist)/4;

Qterminal=500*eye(3);
Rfull=.15*eye(4*n);


for i=1:n
    ewxv=ewxvDiscreteStep([ewxv; uhist( (i-1)*4+1:(i-1)*4+4)]);
end

Jt=(ewxv(7:9)-xF)'*Qterminal*(ewxv(7:9)-xF);
Ju=uhist'*Rfull*uhist;

J=Ju+Jt;


end

