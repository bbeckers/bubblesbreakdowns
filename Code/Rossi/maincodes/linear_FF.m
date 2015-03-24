function result=linear_FF(y,x,R,h,types,bw); 
%Assumptions: OLS estimation, iid std errors
%Our correction 

%OUT-OF-SAMPLE FORECAST TESTS
f_in_seq=[];    f_oos_seq=[];   b_seq_v=[];   yf_seq=[]; 
f_in_split=[];  f_oos_split=[]; b_split_v=[]; yf_split=[]; 
f_in_roll=[];   f_oos_roll=[];  b_roll_v=[];  yf_roll=[];  true=[];

for j=1:rows(y)-R-(h-1); 
   b_seq    =olsbeta(  y(h:R+j-1,1),x(1:R+j-h,:)  );      b_seq_v=[b_seq_v;b_seq];
   yf_seq   =[yf_seq;   x(R+j,:)  *b_seq'      ];
   b_split  =olsbeta(  y(h:R,1), x(1:R-h+1,:) );          b_split_v=[b_split_v;b_split]; 
   yf_split =[yf_split; x(R+j,:)  *b_split'    ];
   b_roll   =olsbeta(  y(j+h-1:R+j-1,1),x(j:R+j-h,:)  );  b_roll_v=[b_roll_v;b_roll]; 
   yf_roll  =[yf_roll;  x(R+j,:)  *b_roll'     ];  
   true     =[true;     y(R+j+h-1,1)]; truetph=y(R+j+h-2,1);
   f_in_seq     = [f_in_seq;    media(( olsres(  y(h:R+j-1,1),x(1:R+j-h,:)  ) ).^2)];     
   f_oos_seq    = [f_oos_seq  ;(y(R+j+h-1,1)-x(R+j,:)  *b_seq')        .^2];    
   f_in_split   = [f_in_split;  media(( olsres(  y(h:R,1), x(1:R-h+1,:) )         ).^2)];       
   f_oos_split  = [f_oos_split;(y(R+j+h-1,1)-x(R+j,:)  *b_split')      .^2]; 
   f_in_roll    = [f_in_roll;   media(( olsres(  y(j+h-1:R+j-1,1),x(j:R+j-h,:)  ) ).^2)]; 
   f_oos_roll   = [f_oos_roll ;(y(R+j+h-1,1)-x(R+j,:)  *b_roll')       .^2]; 
end; 

    %if R==60; rollcv=0.6921; end;
    %if R==100; rollcv=0.9459; end;
    %if R==120; rollcv=1.0955; end;
    %if R==150; rollcv=1.3409; end;
    %if R==175; rollcv=1.5192; end;

% perform FF test
pai=rows(yf_roll)/R; 
[FFseq  ,pvFFseq  ] = FF_test(f_oos_seq,f_in_seq,1,bw);
[FFsplit,pvFFsplit] = FF_test(f_oos_split,f_in_split,1+pai,bw); 
if pai>1; 
[FFroll ,pvFFroll ] = FF_test(f_oos_roll,f_in_roll,2/(3*pai),bw); 
else; 
[FFroll ,pvFFroll ] = FF_test(f_oos_roll,f_in_roll,1-(1/3)*(pai^2),bw); 
end;

%plot(f_oos_roll    -f_in_roll); adsfsdf

%if pvFFseq<0.05;   pvFFseq = 0;   else; pvFFseq=1; end;
%if pvFFsplit<0.05; pvFFsplit = 0; else; pvFFsplit=1; end;
%if pvFFroll<0.05;  pvFFroll= 0;   else; pvFFroll=1; end;


%perform Mincer-Zarnowitz regression
%pvMZ = mincerzarn1(true-yf_roll, ones(rows(true),1) ); 

% if pai>1; 
% pvMZ = mincerzarn1hat_McC(true-yf_roll, ones(rows(true),1) ,2/(3*pai)   ); 
% else; 
% pvMZ = mincerzarn1hat_McC(true-yf_roll, ones(rows(true),1) ,1-(1/3)*(pai^2)   ); 
% end;

%MZ = mincerzarn1hat_McC(true-yf_roll, ones(rows(true),1) ,lambdahh); 
%if MZ>1.3341; pvMZ=0; else; pvMZ=1; end; 

if types=='sequ'; result=[FFseq  , pvFFseq  ]; end;
if types=='spli'; result=[FFsplit, pvFFsplit]; end;
if types=='roll'; result=[FFroll , pvFFroll ]; end;
if types=='allt'; result=[FFseq  , FFsplit,FFroll; pvFFseq, pvFFsplit, pvFFroll ]; end;
if types=='pval'; result=[pvFFseq, pvFFsplit, pvFFroll , pvMZ]; end;
if types=='tval'; result=[FFseq, FFsplit, FFroll , MZ]; end;

%FFroll, pvFFroll, plot(yf_roll); hold on; plot(true,'r'); adsfaf 
