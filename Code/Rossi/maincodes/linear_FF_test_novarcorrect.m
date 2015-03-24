function result=linear_FF_test_novarcorrect(y,x,R,h,types); 
%Assumptions: OLS estimation, iid std errors

%OUT-OF-SAMPLE FORECAST TESTS
f_in_seq=[];    f_oos_seq=[];   b_seq_v=[];   yf_seq=[]; 
f_in_split=[];  f_oos_split=[]; b_split_v=[]; yf_split=[]; 
f_in_roll=[];   f_oos_roll=[];  b_roll_v=[];  yf_roll=[];  true=[];

for j=1:rows(y)-R-(h-1); 
   b_seq    =olsbeta(  y(h:R+j-1,1),x(1:R+j-h,:)  );  b_seq_v=[b_seq_v;b_seq];
   yf_seq   =[yf_seq;   x(R+j-1,:)  *b_seq'      ];
   b_split  =olsbeta(  y(h:R,1), x(1:R-h+1,:) );          b_split_v=[b_split_v;b_split]; 
   yf_split =[yf_split; x(R+j-1,:)  *b_split'    ];
   b_roll   =olsbeta(  y(j+h-1:R+j-1,1),x(j:R+j-h,:)  );  b_roll_v=[b_roll_v;b_roll]; 
   yf_roll  =[yf_roll;  x(R+j-1,:)  *b_roll'     ];  
   true     =[true;     y(R+j+h-1,1)]; truetph=y(R+j+h-2,1);
   f_in_seq     = [f_in_seq;    media(( olsres(  y(h:R+j-1,1),x(1:R+j-h,:)  ) ).^2)];     
   f_oos_seq    = [f_oos_seq  ;(y(R+j+h-1,1)-x(R+j-1,:)  *b_seq')        .^2];    
   f_in_split   = [f_in_split;  media(( olsres(  y(h:R,1), x(1:R-h+1,:) )         ).^2)];       
   f_oos_split  = [f_oos_split;(y(R+j+h-1,1)-x(R+j-1,:)  *b_split')      .^2]; 
   f_in_roll    = [f_in_roll;   media(( olsres(  y(j+h-1:R+j-1,1),x(j:R+j-h,:)  ) ).^2)]; 
   f_oos_roll   = [f_oos_roll ;(y(R+j+h-1,1)-x(R+j-1,:)  *b_roll')       .^2]; 
end; 
    if R==60; rollcv=0.6921; end;
    if R==100; rollcv=0.9459; end;
    if R==120; rollcv=1.0955; end;
    if R==150; rollcv=1.3409; end;
    if R==175; rollcv=1.5192; end;

% perform FF test
[FFseq  ,pvFFseq  ] = FF_test_novarcorrect(f_oos_seq     -f_in_seq   );
[FFsplit,pvFFsplit] = FF_test_novarcorrect(f_oos_split   -f_in_split ); 
[FFroll ,pvFFroll ] = FF_test_novarcorrect(f_oos_roll    -f_in_roll  ); 
if FFseq>1.6619; pvFFseq = 0; else; pvFFseq=1; end;
if FFsplit>2.8210; pvFFsplit = 0; else; pvFFsplit=1; end;
if FFroll>rollcv; pvFFroll= 0; else; pvFFroll=1; end;


%perform Mincer-Zarnowitz regression
%pvMZ = mincerzarn1(true-yf_roll, ones(rows(true),1) ); 
MZ = mincerzarn1hat(true-yf_roll, ones(rows(true),1) ); 
if MZ>1.3341; pvMZ=0; else; pvMZ=1; end; 

if types=='sequ'; result=[FFseq  , pvFFseq  ]; end;
if types=='spli'; result=[FFsplit, pvFFsplit]; end;
if types=='roll'; result=[FFroll , pvFFroll ]; end;
if types=='allt'; result=[FFseq  , FFsplit,FFroll; pvFFseq, pvFFsplit, pvFFroll ]; end;
if types=='pval'; result=[pvFFseq, pvFFsplit, pvFFroll , pvMZ]; end;
if types=='tval'; result=[FFseq, FFsplit, FFroll , MZ]; end;
