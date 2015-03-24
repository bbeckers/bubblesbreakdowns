
del=[0:.05:1];
ns = length(del);

rejfluc_pow1=zeros(ns,1);
rejdm_pow1=zeros(ns,1);
rejopt_pow1=zeros(ns,1);

rejfluc_pow2=zeros(ns,1);
rejdm_pow2=zeros(ns,1);
rejopt_pow2=zeros(ns,1);
for tt=1:ns
    tt
    
    [rejfluc_pow1(tt) rejdm_pow1(tt) rejopt_pow1(tt)]  = MCpower_onetime(150,150,.3*150,del(tt));
        [rejfluc_pow2(tt) rejdm_pow2(tt) rejopt_pow2(tt)]  = MCpower_onetime(150,150,.7*150,del(tt));
end


%Figures

h=plot(del,rejfluc_pow1(:),'--',del,rejdm_pow1(:),'r',del,rejopt_pow1(:),':');
legend('Fluctuation - \mu = .3','GW','One-Time')
axis([0.05 1 0 1]);
%title('TIME VARIATION IN RELATIVE PERFORMANCE');
ylabel('Rejection Frequency')
xlabel('Delta')


figure; 
h=plot(del,rejfluc_pow2(:),'--',del,rejdm_pow2(:),'r',del,rejopt_pow2(:),':');
legend('Fluctuation - \mu = .7','GW','One-Time')
axis([0.05 1 0 1]);
%title('TIME VARIATION IN RELATIVE PERFORMANCE');
ylabel('Rejection Frequency')
xlabel('Delta')


