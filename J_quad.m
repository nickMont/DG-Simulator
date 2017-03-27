function J = J_quad(ewxvuHist)
ewxvu_long2=[ewxvuHist;0;0;0;0]; 

dest=[.5;0;0];

len=length(ewxvu_long2);
N=len/16;
ewxvu=reshape(ewxvu_long2,16,N);

rr=.005;
qq=.1;
RR=diag(ones(3,1))*rr;
QQ=diag(ones(4,1))*qq;

ewxvuE=ewxvu;
ewxvuE(7:9,:)=ewxvuE(7:9,:)-dest*ones(1,N);


xTRx=trace(ewxvuE(7:9,:)'*RR*ewxvuE(7:9,:));
uTQu=trace(ewxvuE(13:16,:)'*QQ*ewxvuE(13:16,:));

J=xTRx+uTQu;


end

