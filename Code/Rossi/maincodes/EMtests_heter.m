function result=EMtests_heter(y,x,z,types);
%This function implements Elliott-Muller's J-test
%INPUT: y=x*b(t)+z*d+e(t) -- if there are no z, put z=0
%types='J-stat' if retrieve Jstat value; 'reject' if want 1 or 0 depending on rej or not rej
%Assumptions: No serial correlation nor heterosked. ( = no HAC correction to Vx)
%OUTPUT: either Elliott-Muller's J-stat (types='J-stat') or 1(test rej.) (types='reject')
if cols(x)>10; result=NaN; %disp('Too many regressors');
else; 
    T=rows(x); if z==0; e=olsres(y,x); else; e=olsres(y,[x,z]); end;
    Vx=0; for s=1:T; Vx=Vx+[e(s,1)^2]*x(s,:)'*x(s,:)/T; end; 
    %sigmae2=cov(e); Vx=(x'*x/T)*sigmae2;
    u=[]; 
    for j=1:cols(x);
        u=[u,x(:,j).*e];
    end; 
    u=((inv(chol(Vx)))*u')';
    wkresv=[]; rbar=1-10/T; rbarvect=cumprod(ones(T,1).*rbar); 
    for j=1:cols(x);
        %u=[u,x(:,j).*e];
        wk=u(1); wkhere=wk; 
        for s=2:T; wkhere=rbar*wkhere+u(s,j)-u(s-1,j);
            wk=[wk;wkhere]; end;
        wkresv=[wkresv,(olsres(wk,rbarvect)).^2];
    end; 
    Jstat=sum(sum(wkresv))*rbar-sum(sum(u.^2)); 
    
    if types=='J-stat';
    result=Jstat;
    end;
    if types=='reject';
    Table1=[ -8.36; -14.32; -19.84; -25.28; -30.60; -35.74; -40.80; -46.18; -51.10; -56.14];
    if Jstat<Table1(cols(x),1); result=1; else; result=0; end;   
    end;
end;     
 