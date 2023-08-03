%% Wave Clus Output Plotting
% load both the raw data and times files manually
%clc
%close all

format long
figure('Position',[0 0 1920 700])
hold on
data=highpass(data,300,25000);
x=[1:length(data)];
x=x./25000;
plot(x,data,'Color','k')
num_clus=(unique(cluster_class(:,1)));
if any(num_clus==0)
    num_clus(num_clus==0)=[];
end

%cluster range inputted by user to for single/multi cluster display
num_clus=[1];

%colors=distinguishable_colors(num_clus,{'w','b'});
colors=colororder;
xlim([0,x(end)])
y_low=min(data);
y_low=y_low(1)-10;
y_hi=max(data);
y_hi=y_hi(1)+10;
ylim([y_low,y_hi])
ylabel("uV")
xlabel("seconds")
%title(strcat(channels_to_plot(k),' FID',string(j)))
for clus=1:length(num_clus)
    pts=unique(cluster_class(cluster_class(:,1)==num_clus(clus),2)./1e3);
    [~,ia,~]=intersect(round(x,8),round(pts,8));
    plot(pts,data(ia),'.','Color',colors(num_clus(clus),:),'MarkerSize',10)
end
ax=gca;
ax.FontSize=20;
title('FID 2 A7')
hold off