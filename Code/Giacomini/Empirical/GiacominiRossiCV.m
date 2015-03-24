function result=GiacominiRossiCV(mu,alpha);
%alpha=0.05 or 0.10  
%mu=m/P


table=[3.393 3.170;
3.179 2.948;
3.012 2.766;
2.890 2.626;
2.779 2.500;
2.634 2.356;
2.560 2.252;
2.433 2.130;
2.248 1.950];


if alpha==0.05; result=table(round(mu*10),1);     
elseif alpha==0.10; result=table(round(mu*10),2);     
end