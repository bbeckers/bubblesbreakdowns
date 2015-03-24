clear
clc
% close all

SL = load('H:\git\bubblesbreakdowns\Results\SLarx.csv');
SLfit = load('H:\git\bubblesbreakdowns\Results\SLarxfit.csv');
SLCI = load('H:\git\bubblesbreakdowns\Results\SLarxCI.csv');
timeline = load('H:\git\bubblesbreakdowns\Results\timeline.csv');

SLCItrim = SLCI;
SLCItrim(SLCItrim<0) = 0;

[N,T] = size(SL);

breakdowns = (SLCI>0)*1;
avgbreakdowns = sum(breakdowns)/N;

x = ones(N,1)*timeline';
y = (1:N)'*ones(1,T);

figure(1)
surf(x,y,SL)
colormap([1  0  0; 0  0  1])
axis tight
title('Surprise losses','FontSize',14)
xlabel('Time','FontSize',12)
ylabel('Model','FontSize',12)
zlabel('Surprise Losses','FontSize',12)

figure(2)
surf(x,y,SLfit)
colormap([1  0  0; 0  0  1])
axis tight
title('Fitted surprise losses','FontSize',14)
xlabel('Time','FontSize',12)
ylabel('Model','FontSize',12)
zlabel('Fitted Surprise Losses','FontSize',12)

figure(3)
surf(x,y,SLCItrim)
colormap([1  0  0; 0  0  1])
axis tight
title('Breakdowns of 1-year IP forecasts','FontSize',14)
xlabel('Time','FontSize',12)
ylabel('Model','FontSize',12)
zlabel('95%-CI of Surprise Losses','FontSize',12)
xlabh = get(gca,'XLabel');
set(xlabh,'Position',get(xlabh,'Position') + [7.5 40 40])
ylabh = get(gca,'YLabel');
set(ylabh,'Position',get(ylabh,'Position') + [5 40 60])

figure(4)
plot(timeline,avgbreakdowns)
axis tight
title('Average number of model forecast breakdowns','FontSize',14)