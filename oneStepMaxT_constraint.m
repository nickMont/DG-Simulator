function [cineq,ceq] = oneStepMaxT_constraint(u)

global umax_GG
umax=umax_GG;

ceq=[];
cineq=norm(u)-umax;



end

