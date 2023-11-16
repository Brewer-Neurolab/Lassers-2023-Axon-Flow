%% Frequency in bins

clear
clc
close all
% parent_dir="D:\Brewer lab data\HFS\No Stim\Wells_5SD_500maxSD";
parent_dir_nostim="D:\Brewer lab data\HFS\no stim\Wells_5SD_500maxSD";
parent_dir_theta="D:\Brewer lab data\HFS\Theta Stim\Wells_5SD_500maxSD";
parent_dir_HFS="D:\Brewer lab data\HFS\HFS Stim\Wells_5SD_500maxSD";
% spike_dyn=load('well_spike_dynamics_table_nostim.mat');
% spike_dyn=spike_dyn.well_spike_dynamics_table;

to_load_nostim=[1 2 3 6 7 9];
cd(parent_dir_nostim)
folders_in_dir=dir;
is_dir=[folders_in_dir.isdir];
folders_names=string({folders_in_dir.name});
folders=folders_names(is_dir);
folders_nostim=folders(3:end);
folders_nostim=folders_nostim(to_load_nostim);

cd(parent_dir_theta)
folders_in_dir=dir;
is_dir=[folders_in_dir.isdir];
folders_names=string({folders_in_dir.name});
folders=folders_names(is_dir);
folders_theta=folders(3:end);

cd(parent_dir_HFS)
folders_in_dir=dir;
is_dir=[folders_in_dir.isdir];
folders_names=string({folders_in_dir.name});
folders=folders_names(is_dir);
folders_HFS=folders(3:end);

load('matching_table_wells_CW.mat')
chan_2_plot=matching_table.electrode;

fs=25000;
time_idx=[1/fs:1/fs:300];

%for i=1:length(folders)
for i=1:1
    cd(strcat(parent_dir_nostim,'\',folders_nostim(i)))
    mkdir('freq_bar')
    cd(strcat(parent_dir_theta,'\',folders_theta(i)))
    mkdir('freq_bar')
    cd(strcat(parent_dir_HFS,'\',folders_HFS(i)))
    mkdir('freq_bar')
    for j=1:length(chan_2_plot)
        try
            nostim_raw=load(strcat(parent_dir_nostim,'\',folders_nostim(i),'\',chan_2_plot(j),'.mat'));
            nostim_spikes=load(strcat(parent_dir_nostim,'\',folders_nostim(i),'\',chan_2_plot(j),'_spikes.mat'));
            theta_raw=load(strcat(parent_dir_theta,'\',folders_theta(i),'\',chan_2_plot(j),'.mat'));
            theta_spikes=load(strcat(parent_dir_theta,'\',folders_theta(i),'\',chan_2_plot(j),'_spikes.mat'));
            HFS_raw=load(strcat(parent_dir_HFS,'\',folders_HFS(i),'\',chan_2_plot(j),'.mat'));
            HFS_spikes=load(strcat(parent_dir_HFS,'\',folders_HFS(i),'\',chan_2_plot(j),'_spikes.mat'));
        catch
            break
        end
        bin_width_vec=[0.1;0.5;1.0;5;10];
        boi_nostim=[];
        boi_5hfs=[];
        boi_40hfs=[];
        
        % 100 ms
        boi1=plot_freq_bins_bar_220210(time_idx,nostim_spikes.index,chan_2_plot(j),0.1,fs); %nostim
        ax1=gca;
        fig1=gcf;
        boi2=plot_freq_bins_bar_220210(time_idx,theta_spikes.index,chan_2_plot(j),0.1,fs); %theta
        ax2=gca;
        fig2=gcf;
        boi3=plot_freq_bins_bar_220210(time_idx,HFS_spikes.index,chan_2_plot(j),0.1,fs); %HFS
        ax3=gca;
        fig3=gcf;
        
        linkaxes([ax1,ax2,ax3],'y')
        
        boi_nostim{1,1}=boi1;
        boi_5hfs{1,1}=boi2;
        boi_40hfs{1,1}=boi3;
        
        saveas(fig1,strcat(parent_dir_nostim,'\',folders_nostim(i),'\freq_bar','\',chan_2_plot(j),'100ms.png'))
        saveas(fig2,strcat(parent_dir_theta,'\',folders_theta(i),'\freq_bar','\',chan_2_plot(j),'100ms.png'))
        saveas(fig3,strcat(parent_dir_HFS,'\',folders_HFS(i),'\freq_bar','\',chan_2_plot(j),'100ms.png'))
        
        % 500 ms
        boi1=plot_freq_bins_bar_220210(time_idx,nostim_spikes.index,chan_2_plot(j),0.5,fs); %nostim
        ax1=gca;
        fig1=gcf;
        boi2=plot_freq_bins_bar_220210(time_idx,theta_spikes.index,chan_2_plot(j),0.5,fs); %theta
        ax2=gca;
        fig2=gcf;
        boi3=plot_freq_bins_bar_220210(time_idx,HFS_spikes.index,chan_2_plot(j),0.5,fs); %HFS
        ax3=gca;
        fig3=gcf;
        
        linkaxes([ax1,ax2,ax3],'y')
        
        boi_nostim{2,1}=boi1;
        boi_5hfs{2,1}=boi2;
        boi_40hfs{2,1}=boi3;
        
        saveas(fig1,strcat(parent_dir_nostim,'\',folders_nostim(i),'\freq_bar','\',chan_2_plot(j),'500ms.png'))
        saveas(fig2,strcat(parent_dir_theta,'\',folders_theta(i),'\freq_bar','\',chan_2_plot(j),'500ms.png'))
        saveas(fig3,strcat(parent_dir_HFS,'\',folders_HFS(i),'\freq_bar','\',chan_2_plot(j),'500ms.png'))
        
        % 1s
        boi1=plot_freq_bins_bar_220210(time_idx,nostim_spikes.index,chan_2_plot(j),1,fs); %nostim
        ax1=gca;
        fig1=gcf;
        boi2=plot_freq_bins_bar_220210(time_idx,theta_spikes.index,chan_2_plot(j),1,fs); %theta
        ax2=gca;
        fig2=gcf;
        boi3=plot_freq_bins_bar_220210(time_idx,HFS_spikes.index,chan_2_plot(j),1,fs); %HFS
        ax3=gca;
        fig3=gcf;
        
        linkaxes([ax1,ax2,ax3],'y')
        
        boi_nostim{3,1}=boi1;
        boi_5hfs{3,1}=boi2;
        boi_40hfs{3,1}=boi3;
        
        saveas(fig1,strcat(parent_dir_nostim,'\',folders_nostim(i),'\freq_bar','\',chan_2_plot(j),'1s.png'))
        saveas(fig2,strcat(parent_dir_theta,'\',folders_theta(i),'\freq_bar','\',chan_2_plot(j),'1s.png'))
        saveas(fig3,strcat(parent_dir_HFS,'\',folders_HFS(i),'\freq_bar','\',chan_2_plot(j),'1s.png'))
        
        % 5s
        boi1=plot_freq_bins_bar_220210(time_idx,nostim_spikes.index,chan_2_plot(j),5,fs); %nostim
        ax1=gca;
        fig1=gcf;
        boi2=plot_freq_bins_bar_220210(time_idx,theta_spikes.index,chan_2_plot(j),5,fs); %theta
        ax2=gca;
        fig2=gcf;
        boi3=plot_freq_bins_bar_220210(time_idx,HFS_spikes.index,chan_2_plot(j),5,fs); %HFS
        ax3=gca;
        fig3=gcf;
        
        linkaxes([ax1,ax2,ax3],'y')
        
        boi_nostim{4,1}=boi1;
        boi_5hfs{4,1}=boi2;
        boi_40hfs{4,1}=boi3;
        
        saveas(fig1,strcat(parent_dir_nostim,'\',folders_nostim(i),'\freq_bar','\',chan_2_plot(j),'5s.png'))
        saveas(fig2,strcat(parent_dir_theta,'\',folders_theta(i),'\freq_bar','\',chan_2_plot(j),'5s.png'))
        saveas(fig3,strcat(parent_dir_HFS,'\',folders_HFS(i),'\freq_bar','\',chan_2_plot(j),'5s.png'))
        
        % 10s
        boi1=plot_freq_bins_bar_220210(time_idx,nostim_spikes.index,chan_2_plot(j),10,fs); %nostim
        ax1=gca;
        fig1=gcf;
        boi2=plot_freq_bins_bar_220210(time_idx,theta_spikes.index,chan_2_plot(j),10,fs); %theta
        ax2=gca;
        fig2=gcf;
        boi3=plot_freq_bins_bar_220210(time_idx,HFS_spikes.index,chan_2_plot(j),10,fs); %HFS
        ax3=gca;
        fig3=gcf;
        
        linkaxes([ax1,ax2,ax3],'y')
        
        boi_nostim{5,1}=boi1;
        boi_5hfs{5,1}=boi2;
        boi_40hfs{5,1}=boi3;
        
        saveas(fig1,strcat(parent_dir_nostim,'\',folders_nostim(i),'\freq_bar','\',chan_2_plot(j),'10s.png'))
        saveas(fig2,strcat(parent_dir_theta,'\',folders_theta(i),'\freq_bar','\',chan_2_plot(j),'10s.png'))
        saveas(fig3,strcat(parent_dir_HFS,'\',folders_HFS(i),'\freq_bar','\',chan_2_plot(j),'10s.png'))
        
        
        t=table(bin_width_vec,boi_nostim,boi_5hfs,boi_40hfs,'VariableNames',{'Bin widths','nostim','5hfs','40hfs'});
        save(strcat(parent_dir_nostim,'\',folders_nostim(i),'\freq_bar','\',chan_2_plot(j),'_table.mat'),'t')
    end
    close all
end