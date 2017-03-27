function [ceq] = trajectoryConstraint(ewxvuLong)
%breaks ewxvu_long into subintervals over which dynamics are linked
%ewxvu_long is the (16*N-4)x1 long vector of stacked ewxvu for each time
%interval such that the u of each part is between intervals
%[ewxvT1 uT1 ewxvT2 uT2...ewxvTn-1 uTn-1 ewxvT]


ewxv0=[0;0;0; 0;0;0; 0;0;0; 0;0;0];
xF=[.5;0;0];

%add blank zeros to end of ewxvulong so that it can be reshaped as a 16xN
ewxvu_long2=[ewxvuLong;0;0;0;0]; 

len=length(ewxvu_long2);
N=len/16;
ewxvu=reshape(ewxvu_long2,16,N);

%initial point
ceq_start=ewxvu(1:12,1)-ewxv0;

%trajectory constraint
ceq_dyn=[];
for i=2:N
    ceq_dyn=[ceq_dyn; ewxvu(1:12,i)-ewxvDiscreteStep(ewxvu(:,i-1))];
end

%terminal point constraint
ceq_end=ewxvu(7:9,N)-xF;

ceq=[ceq_start;ceq_dyn;ceq_end];


end


%fsolve('trajectoryConstraint',[0;0;0;0;0;0; 0;0;0;0;0;0; 5;5;5;5; 0;0;0;0;0;0; .5;0;0;0;0;0])

