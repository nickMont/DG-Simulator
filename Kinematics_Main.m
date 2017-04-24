%Quadcopter dynamics simulator
%UT Austin RNL
clear;clc;
%close all;
digits(32)

%%%%Temporary notes on current testing:
%Control feasiability block needs improvement on maxed thrust case
%Indexing for point clouds in ML needs work


%Simulation time
t0 = 0;
uprate = 10; %Hz, rate of control generation
t_int = 1/uprate; %time interval between steps
plot_int = .002;
tmax = 5;
plottime=tmax; %final time value to print in main simulation graphic

axes=1000*[-1 1 -1 1 -1 1];

N=2;
nu=3; %numinputs
ns=6; %num states, double integrator
x=zeros(ns,N,2);
x(:,1,1)=[0;0;0;0;0;0];
x(:,2,1)=[300;0;0;0;0;0];

%Stack A,B (dynamics matrices) on top of one another
A=zeros(ns,ns,N);
A(:,:,1)=[zeros(3,3) eye(3);zeros(3,3) zeros(3,3)];
A(:,:,2)=[zeros(3,3) eye(3);zeros(3,3) zeros(3,3)];
B=zeros(ns,nu,N);
B(:,:,1)=[zeros(3); eye(3)];
B(:,:,2)=[zeros(3); eye(3)];

%Stack Q,R,K (lqr matrices) on top of one another
Q=zeros(ns,ns,N);
Q(:,:,1)=eye(ns);
Q(:,:,2)=eye(ns);
R=zeros(nu,nu,N);
c1=.5;
c2=.35;
R(:,:,1)=c1*eye(nu);
R(:,:,2)=c2*eye(nu);
a=1; %terminal cost is here a*dot(dx,dx)
K=zeros(ns,ns,N); %terminal cost
K(:,:,1)=a*eye(ns);
K(:,:,2)=a*eye(ns);

if c1>c2
    s=-1;
else
    s=1;
end

figure(1)
clf;
view(3)
scatter3(x(1,:,1), x(2,:,1), x(3,:,1))
axis(axes)


n=0;
for i = t0:t_int:tmax %note: i will be time
    n=n+1;
    
    u=zeros(nu,N);
    %generate controls via LQR
    for j=1:N
        %k=lqr(A(:,:,j),B(:,:,j),Q(:,:,j),R(:,:,j),zeros(length(x(:,j,n),nu));
        %f=@(t,x) (A-B*k)*x
        
        %wlog, p~1 and e~2
        % http://www.dtic.mil/dtic/tr/fulltext/u2/609772.pdf
        if j==1
            u(:,j)=s*c1*(tmax-i)*(x(1:3,1,n)-x(1:3,2,n)+(x(4:6,1,n)-x(4:6,2,n))*(tmax-i))/(1/a+(c1-c2)*(tmax-i)^3/3);
        elseif j==2 
            u(:,j)=s*c2*(tmax-i)*(x(1:3,1,n)-x(1:3,2,n)+(x(4:6,1,n)-x(4:6,2,n))*(tmax-i))/(1/a+(c1-c2)*(tmax-i)^3/3);
        end
        
    end
    
    
    for j=1:N
        AA=A(:,:,j);
        BB=B(:,:,j);
        uu=u(:,j);
        f=@(t,x) AA*x+BB*uu;
        
        [tt,xt]=ode45(f,[i i+t_int], x(:,j,n));
        
        x(:,j,n+1)=(  xt(end,:) )' ;

    end
    
    x(:,:,n+1)
    
    if i<=plottime
        pause(plot_int)
        figure(1)
        clf;
        scatter3(x(1,:,n), x(2,:,n), x(3,:,n))
%         for j=1:N
%             view(3)
%             hold on
%             plot3(x(1,j,n), x(2,j,n), x(3,j,n))
%             
%         end
        axis(axes)
    end
    
    
end



