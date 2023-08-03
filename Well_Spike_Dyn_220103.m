%% Well processing
% Sam Lassers
% Processes and computes spiking dynamics of MCS MEA wells

%% Directory Setup

clear
clc
recorded_date=date;

parent_dir='D:\Brewer lab data\HFS';
%nostim
% mcd_dir='D:\Brewer lab data\HFS\No Stim';
% well_dir='D:\Brewer lab data\HFS\No Stim\Wells_5SD_500maxSD';
%theta
mcd_dir='D:\Brewer lab data\HFS\Theta Stim';
well_dir='D:\Brewer lab data\HFS\Theta Stim\Wells_5SD_500maxSD';
%hfs
mcd_dir='D:\Brewer lab data\HFS\HFS Stim';
well_dir='D:\Brewer lab data\HFS\HFS Stim\Wells_5SD_500maxSD';

%if_cw=[0,0,0,0,1,1,1,1,1];
if_cw=[0,0,0,0,1,1];

%% get sampling frequency and file length

cd(mcd_dir)
filesHere_stim = dir('*.h5');
for i=1:length(filesHere_stim)
    stim_files(i)=McsHDF5.McsData(filesHere_stim(i).name);
    stim_fs(i)=double(1000000./stim_files(1, i).Recording{1, 1}.AnalogStream{1, 1}.Info.Tick(1));
    stim_duration(i)=stim_files(1, i).Recording{1, 1}.Duration/1000000; % In seconds
    stim_times_idx{i}=stim_files(1,i).Recording{1,1}.AnalogStream{1,1}.ChannelDataTimeStamps;
end
stim_scaling_factor=ones(1,length(stim_duration));
%% Wave Clus Batch

% tic
% disp("Calculating stimulated areas...")
% %[stim_areas,dated_output_folder_stim,NMI_Func_Version]=Wave_Clus_NMI_Func_210416(stim_dir,scaling_factor,stim_fs,max_clusters,3,1);
% [dated_output_folder_stim,NMI_Func_Version_stim,stim_times_dir]=...
%     Wave_Clus_batch_220103(mcd_dir,well_dir,stim_fs(1),1,1);
% toc
% 
% cd(dated_output_folder_stim)
% save("well_processing_info","NMI_Func_Version_stim","recorded_date")
% cd(mcd_dir)
% save("well_processing_info_latest","NMI_Func_Version_stim","recorded_date")

%% Convert H5 to mat
file_names=string({filesHere_stim.name});
for i=1:length(file_names)
    convert_h5_into_matfiles_210314(strcat(mcd_dir,'\',file_names(i)),strcat(well_dir,'\',file_names(i)),stim_fs(i))
end

%% Find spikes
cd(well_dir) 
current_dir=dir;
recordings=current_dir([current_dir.isdir]);
recordings=recordings(3:end); %removes relevant pathing
recording_names=string({recordings.name});
for i=1:length(recording_names)
    cd(strcat(well_dir,'\',recording_names(i)))
    current_dir=dir;
    mat_files=current_dir(contains(string({current_dir.name}),".mat") & ~contains(string({current_dir.name}),"_spikes.mat"));
    mat_files=string({mat_files.name});
    for j=1:length(mat_files)
        Get_spikes(char(mat_files(j)))
    end
end

%% Reload
cd(mcd_dir)
%load("well_processing_info_latest")
%% Plot spike shapes
num_rec=length(if_cw);
rec_dir=dir(well_dir);
rec_dir=string({rec_dir(cell2mat({rec_dir.isdir})).name});
rec_dir=rec_dir(3:end);
for i=1:num_rec
    well_plot_spike_shapes_220125(strcat(well_dir,'\',rec_dir(i)),well_dir,if_cw(i))
    close all
end
%% Peak train assembly and spike dynamics

[all_regions_unit_matched]=directory_peak_train_assembly_220103(well_dir,stim_times_idx,if_cw);
save(strcat(well_dir,'\',"allregion_unit_matched.mat"),"all_regions_unit_matched",'-v7.3')
[well_spike_dynamics_table]=well_compute_spike_burst_dynamics_axons_220105(all_regions_unit_matched);
cd(well_dir)
%save("well_spike_dynamics_table_nostim",'well_spike_dynamics_table')
%save("well_spike_dynamics_table_theta",'well_spike_dynamics_table')
save("well_spike_dynamics_table_hfs",'well_spike_dynamics_table')