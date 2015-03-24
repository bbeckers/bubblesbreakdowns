R=50;%[20 50 150];
nr=length(R);
n = 50;%[20 50 150];
nn=length(n);
mu=0.4;%[.1 .3 .5  .7 .9];

nmu=length(mu);

rejfluc_size=zeros(nr,nn,nmu);
rejdm_size=zeros(nr,nn);
rejopt_size=zeros(nr,nn);


for ii=1:nr
    for jj=1:nn
        for ss=1:nmu
            m = mu(ss)*n(jj);
       [rejfluc_size(ii,jj,ss) rejdm_size(ii,jj) rejopt_size(ii,jj)]  = MCsize_new(R(ii),n(jj),m);
          
          
        end

    end
end
rejfluc_size
rejdm_size
rejopt_size