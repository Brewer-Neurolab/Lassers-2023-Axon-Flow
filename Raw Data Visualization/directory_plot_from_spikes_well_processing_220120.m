%% Wave Clus Output Plotting
% load both the raw data and times files manually
%this is for wells
%clc
%close all

% current_dir=dir;
% spikes_files={current_dir.name};
% spikes_files=string(spikes_files);
% spikes_files=spikes_files(contains(spikes_files,'_spikes.mat'));
% for i=1:length(spikes_files)
%     split_name=strsplit(spikes_files(i),{'_'});
%     names(i)=split_name(1);
% end
clear
clc
close all
load('..\matching_table_wells_CW')
names=matching_table.electrode;

for i=1:length(names)
    load(strcat(names(i),'.mat'))
    load(strcat(names(i),'_spikes.mat'))
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
    yline(-threshold*10,'Color','r')
    %title(strcat(channels_to_plot(k),' FID',string(j)))
    pts=unique((index)./1e3);
    [~,ia,~]=intersect(round(x,8),round(pts,8));
    plot(pts,data(ia),'.','Color',[0 0 1],'MarkerSize',10)
    ax=gca;
    ax.FontSize=20;
    %yline(threshold)
    title(names(i))
    hold off
    %saveas(gcf,strcat('D:\Brewer lab data\HFS\No stim Raw Data Well Figs 5SD\FID10','\',names(i),'.png'))
    %saveas(gcf,strcat('D:\Brewer lab data\HFS\Theta Raw Data Well Figs 5SD\FID10','\',names(i),'.png'))
    saveas(gcf,strcat('D:\Brewer lab data\HFS\HFS Raw Data Well Figs 5SD\FID10','\',names(i),'.png'))    
end