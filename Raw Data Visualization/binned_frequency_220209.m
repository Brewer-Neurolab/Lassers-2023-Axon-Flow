%% Frequency in bins

clear
clc
% parent_dir="D:\Brewer lab data\HFS\No Stim\Wells_5SD_500maxSD";
parent_dir="D:\Brewer lab data\HFS\HFS Stim\Wells_5SD_500maxSD";
cd(parent_dir)
% spike_dyn=load('well_spike_dynamics_table_nostim.mat');
% spike_dyn=spike_dyn.well_spike_dynamics_table;

folders_in_dir=dir;
is_dir=[folders_in_dir.isdir];
folders_names=string({folders_in_dir.name});
folders=folders_names(is_dir);
folders=folders(3:end);

load('matching_table_wells_CW.mat')
chan_2_plot=matching_table.electrode;

fs=25000;
time_idx=[1/fs:1/fs:300];

%for i=1:length(folders)
for i=1:1
    cd(strcat(parent_dir,'\',folders(i)))
    mkdir('freq_bar')
    for j=1:length(chan_2_plot)
        try
            load(strcat(chan_2_plot(j),'.mat'));
            load(strcat(chan_2_plot(j),'_spikes.mat'));
        catch
            break
        end
        hold on
        plot(time_idx,data,'k')
        yline(-threshold,'r')
        yline(-threshold*100,'r')
        hold off
        
        ylim([min(data)-10,max(data)+10])
        
        % 100ms bin
        bins_edges_100ms=[time_idx(1):time_idx(0.1*fs):time_idx(end),time_idx(end)];
        index=index/1e3;
        vec_2_plot=[];
        for pts=1:length(bins_edges_100ms)-1
            freq=sum(index>bins_edges_100ms(pts)&index<=bins_edges_100ms(pts+1))/.1;
            vec_2_plot(pts)=freq;
        end
        
        figure
        centers=convert_edges_2_centers(bins_edges_100ms);
        bar(centers,vec_2_plot)
        title("100ms Bins")
        xlabel("Seconds")
        ylabel("Frequency")
%         ax=gca;
%         ax.YScale='log';
        saveas(gcf,strcat('freq_bar\',chan_2_plot(j),'_100ms.png'))
        
        % 200ms bin
        bins_edges_200ms=[time_idx(1):time_idx(0.2*fs):time_idx(end),time_idx(end)];
        %index=index/1e3;
        vec_2_plot=[];
        for pts=1:length(bins_edges_200ms)-1
            freq=sum(index>bins_edges_200ms(pts)&index<=bins_edges_200ms(pts+1))/.2;
            vec_2_plot(pts)=freq;
        end
        
        figure
        centers=convert_edges_2_centers(bins_edges_200ms);
        bar(centers,vec_2_plot)
        title("200ms Bins")
        xlabel("Seconds")
        ylabel("Frequency")
%         ax=gca;
%         ax.YScale='log';
        saveas(gcf,strcat('freq_bar\',chan_2_plot(j),'_200ms.png'))
        
        % 500ms bin
        bins_edges_500ms=[time_idx(1):time_idx(0.5*fs):time_idx(end),time_idx(end)];
        %index=index/1e3;
        vec_2_plot=[];
        for pts=1:length(bins_edges_500ms)-1
            freq=sum(index>bins_edges_500ms(pts)&index<=bins_edges_500ms(pts+1))/.5;
            vec_2_plot(pts)=freq;
        end
        
        figure
        centers=convert_edges_2_centers(bins_edges_500ms);
        bar(centers,vec_2_plot)
        title("500ms Bins")
        xlabel("Seconds")
        ylabel("Frequency")
%         ax=gca;
%         ax.YScale='log';
        saveas(gcf,strcat('freq_bar\',chan_2_plot(j),'_500ms.png'))
        
        % 1s bin
        bins_edges_1s=[time_idx(1):time_idx(1*fs):time_idx(end),time_idx(end)];
        %index=index/1e3;
        vec_2_plot=[];
        for pts=1:length(bins_edges_1s)-1
            freq=sum(index>bins_edges_1s(pts)&index<=bins_edges_1s(pts+1))/1;
            vec_2_plot(pts)=freq;
        end
        
        figure
        centers=convert_edges_2_centers(bins_edges_1s);
        bar(centers,vec_2_plot)
        title("1s Bins")
        xlabel("Seconds")
        ylabel("Frequency")
        saveas(gcf,strcat('freq_bar\',chan_2_plot(j),'_1s.png'))
        
        % 30s bin
        bins_edges_30s=[time_idx(1):time_idx(30*fs):time_idx(end),time_idx(end)];
        %index=index/1e3;
        vec_2_plot=[];
        for pts=1:length(bins_edges_30s)-1
            freq=sum(index>bins_edges_30s(pts)&index<=bins_edges_30s(pts+1))/30;
            vec_2_plot(pts)=freq;
        end
        
        figure
        centers=convert_edges_2_centers(bins_edges_30s);
        bar(centers,vec_2_plot)
        title("1s Bins")
        xlabel("Seconds")
        ylabel("Frequency")
        saveas(gcf,strcat('freq_bar\',chan_2_plot(j),'_30s.png'))
    end
    close all
end