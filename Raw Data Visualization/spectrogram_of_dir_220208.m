%% Create spectograph of raw data from a directory input for WELLS
% source: http://www-users.med.cornell.edu/~jdvicto/pdfs/pubo08.pdf

%incomplete program

clear
clc
parent_dir="D:\Brewer lab data\HFS\No Stim\Wells_5SD_500maxSD";
cd(parent_dir)
spike_dyn=load('well_spike_dynamics_table_nostim.mat');
spike_dyn=spike_dyn.well_spike_dynamics_table;

folders_in_dir=dir;
is_dir=[folders_in_dir.isdir];
folders_names=string({folders_in_dir.name});
folders=folders_names(is_dir);
folders=folders(3:end);

load('matching_table_wells_CW.mat')
chan_2_plot=matching_table.electrode;

time_idx=[1/25000:1/25000:300];

for i=1:length(folders)
    cd(strcat(parent_dir,'\',folders(i)))
    for j=1:chan_2_plot
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
        % 1s bin
        target_hz=spike_dyn.SpikeRate(spike_dyn.fi==i && contains(spike_dyn.channel_name,chan_2_plot(j)));
        freq_vec_x
    end
end