%% scatter plot of cluster times

close all
clear

up=load('times_A7_hi.mat');
down=load('times_A6_hi.mat');
sr=25000;
timevec=[1/sr:1/sr:300];

up_1=up.cluster_class(up.cluster_class==1,2)./1e3;
up_2=up.cluster_class(up.cluster_class==2,2)./1e3;

down_1=down.cluster_class(down.cluster_class==1,2)./1e3;
down_2=down.cluster_class(down.cluster_class==2,2)./1e3;

xbin_edges=linspace(-1,1,51);

figure
[locA,locB]=ismembertol(up_1,down_1,0.001,'DataScale',1);
diff1=diff([down_1(locB(locB~=0)),up_1(locA)],1,2)*1e3;
diffvec=repmat(NaN,1,length(timevec));
[~,intersect_up,~]=intersect(round(timevec,10),round(up_1(locA),10));
diffvec(intersect_up)=diff1;
hold on
scatter(timevec,diffvec,'k')
scatter(up_1,ones(1,length(up_1)),'MarkerEdgeColor',[0.4660 0.6740 0.1880])
scatter(down_1,-ones(1,length(down_1)),'MarkerEdgeColor',[0.4940 0.1840 0.5560])
hold off
title('A7 Clus 1 A6 Clus 1')
xlabel('Time (s)')
ylabel('Delay from A7 to A6 (ms)')
ylim([-1,1])
set(gca,'FontSize',18)
%scatter(up_1(locA),down_1(locB(locB~=0)))
% xlabel('A7 Clus 1 time (s)')
% ylabel('A6 Clus 1 time (s)')
saveas(gcf,'A7Clus1A6Clus1hi_scatter.png')
figure
histogram(diff1,xbin_edges)
title('A7 Clus 1 A6 Clus 1')
xlabel('Delay (ms)')
ylabel('Spike Count')
set(gca,'FontSize',18)
saveas(gcf,'A7Clus1A6Clus1hi_hist.png')

figure
[locA,locB]=ismembertol(up_1,down_2,0.001,'DataScale',1);
diff2=diff([down_2(locB(locB~=0)),up_1(locA)],1,2)*1e3;
diffvec=repmat(NaN,1,length(timevec));
[~,intersect_up,~]=intersect(round(timevec,10),round(up_1(locA),10));
diffvec(intersect_up)=diff2;
hold on
scatter(timevec,diffvec,'k')
scatter(up_1,ones(1,length(up_1)),'MarkerEdgeColor',[0.4660 0.6740 0.1880])
scatter(down_2,-ones(1,length(down_2)),'MarkerEdgeColor',[0.4940 0.1840 0.5560])
hold off
title('A7 Clus 1 A6 Clus 2')
xlabel('Time (s)')
ylabel('Delay from A7 to A6 (ms)')
ylim([-1,1])
set(gca,'FontSize',18)
% scatter(up_1(locA),down_2(locB(locB~=0)))
% xlabel('A7 Clus 1 time (s)')
% ylabel('A6 Clus 2 time (s)')
saveas(gcf,'A7Clus1A6Clus2hi_scatter.png')
figure
histogram(diff2,xbin_edges)
title('A7 Clus 1 A6 Clus 2')
xlabel('Delay (ms)')
ylabel('Spike Count')
set(gca,'FontSize',18)
saveas(gcf,'A7Clus1A6Clus2hi_hist.png')

figure
[locA,locB]=ismembertol(up_2,down_1,0.001,'DataScale',1);
diff3=diff([down_1(locB(locB~=0)),up_2(locA)],1,2)*1e3;
diffvec=repmat(NaN,1,length(timevec));
[~,intersect_up,~]=intersect(round(timevec,10),round(up_2(locA),10));
diffvec(intersect_up)=diff3;
hold on
scatter(timevec,diffvec,'k')
scatter(up_2,ones(1,length(up_2)),'MarkerEdgeColor',[0.4660 0.6740 0.1880])
scatter(down_1,-ones(1,length(down_1)),'MarkerEdgeColor',[0.4940 0.1840 0.5560])
hold off
title('A7 Clus 2 A6 Clus 1')
xlabel('Time (s)')
ylabel('Delay from A7 to A6 (ms)')
ylim([-1,1])
set(gca,'FontSize',18)
% scatter(up_2(locA),down_1(locB(locB~=0)))
% xlabel('A7 Clus 2 time (s)')
% ylabel('A6 Clus 1 time (s)')
saveas(gcf,'A7Clus2A6Clus1hi_scatter.png')
figure
histogram(diff3,xbin_edges)
title('A7 Clus 2 A6 Clus 1')
xlabel('Delay (ms)')
ylabel('Spike Count')
set(gca,'FontSize',18)
saveas(gcf,'A7Clus2A6Clus1hi_hist.png')

figure
[locA,locB]=ismembertol(up_2,down_2,0.001,'DataScale',1);
diff4=diff([down_2(locB(locB~=0)),up_2(locA)],1,2)*1e3;
diffvec=repmat(NaN,1,length(timevec));
[~,intersect_up,~]=intersect(round(timevec,10),round(up_2(locA),10));
diffvec(intersect_up)=diff4;
hold on
scatter(timevec,diffvec,'k')
scatter(up_2,ones(1,length(up_2)),'MarkerEdgeColor',[0.4660 0.6740 0.1880])
scatter(down_2,-ones(1,length(down_2)),'MarkerEdgeColor',[0.4940 0.1840 0.5560])
hold off
title('A7 Clus 2 A6 Clus 2')
xlabel('Time (s)')
ylabel('Delay from A7 to A6 (ms)')
ylim([-1,1])
set(gca,'FontSize',18)
% scatter(up_2(locA),down_2(locB(locB~=0)))
% xlabel('A7 Clus 2 time (s)')
% ylabel('A6 Clus 1 time (s)')
saveas(gcf,'A7Clus2A6Clus2hi_scatter.png')
figure
histogram(diff4,xbin_edges)
title('A7 Clus 2 A6 Clus 2')
xlabel('Delay (ms)')
ylabel('Spike Count')
set(gca,'FontSize',18)
saveas(gcf,'A7Clus2A6Clus2hi_hist.png')

