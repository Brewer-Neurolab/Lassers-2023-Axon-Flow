%% Wave Clus Output Plotting
% load both the raw data and times files manually
%this is for wells
%clc
%close all

format long
figure('units','normalized','OuterPosition',[0 0 1 1])
hold on
data=highpass(data,300,25000);
x=[1:length(data)];
x=x./25000;
plot(x,data,'Color','k')

%cluster range inputted by user to for single/multi cluster display

%colors=distinguishable_colors(num_clus,{'w','b'});
xlim([0,x(end)])
y_low=min(data);
y_low=y_low(1)-10;
y_hi=max(data);
y_hi=y_hi(1)+10;
ylim([y_low,y_hi])
ylabel("uV")
xlabel("seconds")
yline(-threshold,'Color','r')
%title(strcat(channels_to_plot(k),' FID',string(j)))
pts=unique((index)./1e3);
[~,ia,~]=intersect(round(x,8),round(pts,8));
plot(pts,data(ia),'.','Color',[0 0 1],'MarkerSize',10)
ax=gca;
ax.FontSize=20;
%yline(threshold)
title('FID 1 C2')
hold off