function [xbar_up,Pbar_up,chi_pts,chi_bar_pts] = unscented_transform(f,xbar,P,a,B,K)

n=length(xbar);
S=chol(P)';
lambda=a^2*(n+K)-n;

chi_pts=zeros(n,1+2*n);
chi_pts(:,1)=xbar;

%Generate chi-values and then propagate.
for i=1:n
   chi_pts(:,1+i) = xbar+sqrt(n+lambda)*S(:,i);
end
for i=n+1:2*n
    chi_pts(:,1+i)=xbar-sqrt(n+lambda)*S(:,i-n);
end

%propagate through f
chi_bar_pts=zeros(n,1+2*n);
for i=1:2*n+1
    chi_bar_pts(:,i)=feval(f,chi_pts(:,i));    
end

%Recalculate xbar,Pbar
W0m=lambda/(lambda+n);
Wim=1/(2*(lambda+n));
W0c=lambda/(lambda+n)+1-a^2+B;
Wic=1/(2*(lambda+n));

xbar_sum=zeros(n,1);
P_sum=zeros(n,n);
for i=1:2*n+1
    if i==1
        Wm=W0m;
    else
        Wm=Wim;
    end
    xbar_sum=xbar_sum+Wm*chi_bar_pts(:,i);
end
xbar_up = xbar_sum;

for i=1:2*n+1
    if i==1
        Wc=W0c;
    else
        Wc=Wic;
    end
    P_sum=P_sum+ Wc*(chi_bar_pts(:,i)-xbar_up)*(chi_bar_pts(:,i)-xbar_up)';
end
Pbar_up=P_sum;


end

