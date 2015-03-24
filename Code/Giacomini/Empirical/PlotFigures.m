clear; 
warning off all; Fluc=1;

load Taylorrulefundamentalsdata_updated20080217.txt;
data=Taylorrulefundamentalsdata_updated20080217;
s_all=data(:,1:12);   m1_all=data(:,14:28);  i_all=data(:,30:39); 
y_all=data(:,41:53);  p_all=data(:,55:67);

country_name=['JAP';'CAN';'SWI';'GBP';'FRA';'GER';'ITA';'SWE';'AUS';'DEN';'NET';'POR'];
country_name1=['JAP ';'CAN ';'SWI ';'GBP ';'FRA ';'GER ';'ITA ';'SWE ';'AUS ';'DEN ';'NET ';'POR '];
country_name2=['Japan     ';'Canada    ';'Switzer.  ';'U. K.     ';'France    ';'Germany   ';'Italy     ';'Sweden    ';'Australia ';'Denmark   ';'Netherl.  ';'Portugal  '];
m_US=m1_all(:,1);  m_diff=-log(m1_all(:,4:15))+kron(log(m_US),ones(1,12)); 
y_US=y_all(:,1);   y_diff=log(y_all(:,2:13))-kron(log(y_US),ones(1,12)); y_diff_level=(y_all(:,2:13))-kron((y_US),ones(1,12));
i_US=i_all(:,1);   i_diff=[i_all(:,10),i_all(:,2:9)-kron(i_US,ones(1,8)),NaN*ones(rows(log(y_US)),3)];
p_US=p_all(:,1);   p_diff=-log(p_all(:,2:13))+kron(log(p_US),ones(1,12));
y_diff_alldata=y_diff; y_US_alldata=log(y_all(:,1)); y_nonUS_alldata=log(y_all(:,2:13));
y_US_alldata_level=(y_all(:,1)); y_nonUS_alldata_level=(y_all(:,2:13));

        %Calculate series of output gap
        ygap=[];  
        for j=1:cols(y_diff_alldata); ygapcountry=[]; 
            for s=1:rows(y_nonUS_alldata); X=[ones(s,1),cumsum(ones(s,1))];
              if j~=10; 
                yd=-y_nonUS_alldata(:,j)+y_US_alldata; 
                res=yd(s,:)-X(s,:)*(inv(X'*X))*X'*yd(1:s,:); 
              elseif j==10; 
                res=nan; 
              end;
             ygapcountry=[ygapcountry;res];  
            end; ygap=[ygap,ygapcountry];
        end;     

init=25; %Start from 1973:3
ygap=ygap(2+init:end-1,:);   
infl_US=log(p_US(2+init:end-1,:))-log(p_US(1+init-11:end-2-11,:)); infl_nonUS=log(p_all(2+init:end-1,2:end))-log(p_all(1+init-11:end-2-11,2:end));
infl_diff= kron(infl_US,ones(1,12))-infl_nonUS;
m_US=m_US(2+init:end-1);  m_diff=m_diff(2+init:end-1,:);  y_US=y_US(2+init:end-1);  y_diff=y_diff(2+init:end-1,:); y_diff_level=y_diff_level(2+init:end-1,:);
i_US=i_US(2+init:end-1);  i_diff=i_diff(2+init:end-1,:);  p_US=p_US(2+init:end-1);  p_diff=p_diff(2+init:end-1,:); 
Ds_tplus1=log(s_all(3+init:end,:))-log(s_all(2+init:end-1,:)); s_t=-log(s_all(2+init:end-1,:));
time=calendar(1973,3,418,'m'); tds=calendar_plot(1973,3,418,'m'); %R=108; %first estimation period is 1973:3-1982:2; 
startmonth=2; startyear=1983; %first estimation period is 1973:3-1982:2; 

TableMSE_all=[]; TableDMall=[]; TableCWall=[]; TableOPTall=[];
for model=[1,4]; 1:4;  %1=Taylor I, 2=monetary, 3=PPP, 4=UIP
    table_MSE=[]; table_DM=[]; table_CW=[]; NC=12; table_opt=[];
    for ncountry=1:NC; 
        e1=Ds_tplus1(:,ncountry); 
             
        if model==1; p1=[ygap(:,ncountry),infl_diff(:,ncountry)]; gap=1; end; 
        if model==2;   
            if ncountry==4; p1=[m_diff(:,ncountry)+s_t(:,ncountry)];
            elseif ncountry==9; p1=[m_diff(:,ncountry)+s_t(:,ncountry)];
            else;
            p1=[m_diff(:,ncountry)-s_t(:,ncountry)];;
            end; 
        end;
        if model==3; 
            if ncountry==4; p1=[p_diff(:,ncountry)+s_t(:,ncountry)];
            elseif ncountry==9; p1=[p_diff(:,ncountry)+s_t(:,ncountry)];
            else;
            p1=[p_diff(:,ncountry)-s_t(:,ncountry)];
            end; 
        end; 
        if model==4; p1=[i_diff(:,ncountry)]; end; 
        
        if gap==1; 
                %Calculate series of output gap forecasts based on the estimated trend
                trendp=[]; ygapcountry=[]; ygapcountryf=nan;
                    for s=1:rows(y_diff_alldata); X=[ones(s,1),cumsum(ones(s,1))];
                      if ncountry~=10; 
                        yd=-y_nonUS_alldata(:,ncountry)+y_US_alldata; 
                        res=yd(s,:)-X(s,:)*(inv(X'*X))*X'*yd(1:s,:); 
                        if s~=rows(y_diff_alldata);
                          resf=yd(s+1,:)-[1,s+1]*(inv(X'*X))*X'*yd(1:s,:); 
                        end; 
                      elseif ncountry==10; 
                        resf=nan; 
                      end;
                     ygapcountry=[ygapcountry;res]; 
                     if s~=rows(y_diff_alldata); ygapcountryf=[ygapcountryf;resf]; end;
                    end; 
                    ygapcountry=ygapcountry(2+init:end-1,:); 
                    ygapcountryf=ygapcountryf(2+init:end-1,:);            
        end;
        
        if model==1 & ncountry==3;  
            result=[nan, nan; nan, nan]; 
        elseif model==1 & ncountry==10;
            result=[nan, nan; nan, nan]; 
        elseif model==3 & ncountry==9;
            result=[nan, nan; nan, nan]; 
        elseif model==4 & ncountry==10;  
            result=[nan, nan; nan, nan]; 
        elseif model==4 & ncountry==11;  
            result=[nan, nan; nan, nan]; 
        elseif model==4 & ncountry==12;  
            result=[nan, nan; nan, nan]; 
        else;
            if model==1 & ncountry==9; p1=p1(:,2); result=[nan, nan; nan, nan]; end; %country 9 is missing output gap
            
            if ncountry~=9; m=50; 
               if Fluc==1;  
                   axis auto;
                result=testsoos3(e1,p1,startyear,startmonth,time,tds,m);  
                figure; plot(result(:,5),result(:,1)); titlestr=['title(''DMW in ',country_name2(ncountry,:),''');']; eval(titlestr); 
                hold on; plot(result(:,5),result(:,6),'r:','LineWidth',3); 
                ylabel('Relative Performance'); 
                if model==4; text(0,0,'UIRP model'); elseif model==1; text(0,0,'Taylor rule model'); end;
                xlabel('Time'); 
                legend('Fluctuation test','Fluctuation test Critical Value');
                hold off; 
                figure; plot(result(:,5),result(:,2)); titlestr=['title(''CW in ',country_name2(ncountry,:),''');']; eval(titlestr);
                hold on; plot(result(:,5),result(:,6),'r:','LineWidth',3);  
                %hold on; plot(result(:,5),result(:,7),'r-'); 
                if model==4; ylabel('UIRP model'); elseif model==1; ylabel('Taylor rule model'); end;
                xlabel('Time'); hold off; 
                legend('Fluctuation test Statistic','Fluctuation test Critical Value'); figure; 
               elseif Fluc==0;
                result=testsoosBREAK(e1,p1,startyear,startmonth,time,tds,m,ncountry,model,country_name);                     
               end;
            end;
            if model==4; if ncountry==4;  
                result22=testsoosBREAK(e1,p1,startyear,startmonth,time,tds,m,ncountry,model,country_name);                     
                hold off; figure; 
                result=testsoos3(e1,p1,startyear,startmonth,time,tds,m);  
                figure; plot(result(:,5),result(:,1),'LineWidth',3); titlestr=['title(''DMW in ',country_name(ncountry,:),' and model ',num2str(model),''');']; eval(titlestr); 
                hold on; 
                load tbreaksaved; yearbreak=floor(tbreaksaved); monthbreak=round((tbreaksaved-yearbreak)*10); 
                tbreak=find(result(:,5)==tbreaksaved); 
                plot(result(:,5),[mean(result(1:tbreak,1))*ones(tbreak,1); mean(result(tbreak+1:end,1))*ones(rows(result)-tbreak,1)],'r:','LineWidth',3);
                legend('Fluctuation test','Result based on Test against One-time Reversal');
                title('UIRP in the U.K.'); xlabel('Time'); 

                hold off; figure; 
                result=testsoos3(e1,p1,startyear,startmonth,time,tds,m);  
                figure; plot(result(:,5),result(:,1),'LineWidth',2);  
                legend('Rolling MSFE difference (standardized)'); 
                ylabel('Rolling MSFE difference (standardized)'); xlabel('Time');

                hold off; figure; 
                result=testsoos3(e1,p1,startyear,startmonth,time,tds,m);  
                figure; plot(result(:,5),result(:,1),'LineWidth',3);  
                hold on; plot(result(:,5),result(:,6),'r:','LineWidth',3); 
                load tbreaksaved; yearbreak=floor(tbreaksaved); monthbreak=round((tbreaksaved-yearbreak)*10); 
                tbreak=find(result(:,5)==tbreaksaved); 
                plot(result(:,5),[mean(result(1:tbreak,1))*ones(tbreak,1); mean(result(tbreak+1:end,1))*ones(rows(result)-tbreak,1)],'k--','LineWidth',3);
                legend('Fluctuation test Statistic','Fluctuation test Critical Value','Result based on Test against One-time Reversal');
                xlabel('Time');
                end; end;
            
        end; 
        table_opt=[table_opt,result(2,1)];
    end; 
TableOPTall =[TableOPTall ;   table_opt]; 
end;     
  