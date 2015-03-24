function wald = chowgmm(y,x,z,w,t,heter);
%Calculates Chow test statistic by GMM estimation for a fixed time-break t
%y is (T*1) dependent variable
%x is (T*p) independent variable whose param are time-varying
%z is (T*q) independent variable whose param are not time-varying
%w is (T*k) vector of instruments
%if no subsets then set z=0;
%t is the time of the break (NOT as fraction of sample size) -- ie two samples are 1:t-1,t:T
%const=0 there is no constant among w; =column where constant is if there is a constant
%the const is needed otherwise multicollinearity from splitting w in two subsets
%From Andrews, 1993

if z==0;
   
%Andrews' GMM formula
[T,p]=size(x); pi=t/T; [T,k]=size(w);
y1=y(1:t,:); x1=x(1:t,:); w1=w(1:t,:); 
y2=y(t+1:T,:); x2=x(t+1:T,:); w2=w(t+1:T,:); 
X=[x1,zeros(t,p);zeros(T-t,p),x2]; W=[w1,zeros(t,k);zeros(T-t,k),w2];
b=gmmbeta(y,X,W,heter); b1=b(1:p); b2=b(p+1:2*p);
e=gmmres(y,X,W,heter); 
M1=w1'*x1./t; 
M2=w2'*x2./(T-t); 
if heter==1;
S1=w1'*(diag(e(1:t).^2))*w1/t;
S2=w2'*(diag(e(t+1:T).^2))*w2/(T-t);
elseif heter==0; s2=(sum(e(1:T).^2))/(T); S1=s2*w1'*w1/t; S2=s2*w2'*w2/(T-t); 
end;
V1=inv(M1'*(inv(S1))*M1); V2=inv(M2'*(inv(S2))*M2);
wald=T*(b1-b2)'*(inv(V1/pi+V2/(1-pi)))*(b1-b2);

else;

%Andrews' GMM formula
[T,p]=size(x); [T,q]=size(z); pi=t/T; [T,k]=size(w);
y1=y(1:t,:); x1=x(1:t,:); z1=z(1:t,:); w1=w(1:t,:);
y2=y(t+1:T,:); x2=x(t+1:T,:); z2=z(t+1:T,:); w2=w(t+1:T,:);
X=[x1,zeros(t,p),z1;zeros(T-t,p),x2,z2]; W=[w1,zeros(t,k);zeros(T-t,k),w2];
b=gmmbeta(y,X,W,heter); b1=b(1:p); b2=b(p+1:2*p);
e=gmmres(y,X,W,heter);
M1=w1'*x1./t; 
M2=w2'*x2./(T-t);
if heter==1; S1=w1'*(diag(e(1:t).^2))*w1/t; S2=w2'*(diag(e(t+1:T).^2))*w2/(T-t);
elseif heter==0; s2=(sum(e(1:T).^2))/(T);S1=s2*w1'*w1/t;S2=s2*w2'*w2/(T-t);
end;
V1=inv(M1'*(inv(S1))*M1); V2=inv(M2'*(inv(S2))*M2);
wald=T*(b1-b2)'*(inv(V1/pi+V2/(1-pi)))*(b1-b2);

end;
