sigpower=[.1:.1:1];

ns = length(sigpower);

rejfluc_powconst1=zeros(ns,1);
rejfluc_powconst1=zeros(ns,1);
rejdm_powconst1=zeros(ns,1);
rejdm_powconst2=zeros(ns,1);
rejopt_powconst1=zeros(ns,1);
rejopt_powconst2=zeros(ns,1);


load MCxdata.txt;
x = MCxdata;
b=0.05*ones(300,1);
for tt=1:ns
    tt
    [rejfluc_powconst1(tt) rejdm_powconst1(tt) rejopt_powconst1(tt)]  = MCpower_const(150,150,.3*150,x,b,sigpower(tt));
        [rejfluc_powconst2(tt) rejdm_powconst2(tt) rejopt_powconst2(tt)]  = MCpower_const(150,150,.7*150,x,b,sigpower(tt));
end
%figures

h=plot(sigpower,rejfluc_powconst1(:),'--',sigpower,rejdm_powconst1(:),'r',sigpower,rejopt_powconst1(:),':');
legend('Fluctuation - \mu = .3','GW','One-Time Reversal')
axis([0.05 1 0 1]);
%title('UNEQUAL BUT CONSTANT RELATIVE PERFORMANCE - MU =.3');
ylabel('Rejection Frequency')
xlabel('\sigma')


figure;
h=plot(sigpower,rejfluc_powconst2(:),'--',sigpower,rejdm_powconst2(:),'r',sigpower,rejopt_powconst2(:),':');
legend('Fluctuation - \mu = .7','GW','One-Time Reversal')
axis([0.05 1 0 1]);
%title('UNEQUAL BUT CONSTANT RELATIVE PERFORMANCE - \MU =.7');
ylabel('Rejection Frequency')
xlabel('\sigma')