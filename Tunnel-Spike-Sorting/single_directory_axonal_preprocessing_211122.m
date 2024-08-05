%% Comparison_Spontaneous_Stimulated_data
% Sam Lassers
% Uses inputs from h5 MCS files in order to perform a t test and linear
% regression. Input directories with h5 files converted from spontaneous
% and stimulated data sets from MCS.

% Change log starting 6/2/21 because this has gotten software sized
% added scaling functionality 6/2/21
% changed scaling to recording duration for frequency (spikes per second)
% to normalize across different recording lengths
% Implemented interunit peak finding 6/24/21
% Separated wave clus and NMI functions 6/24/21
%% Setup
clear
clc
recorded_date=date;
% stim_dir='D:\Brewer lab data\Yash\Filtered Stim Data MCD Files';
% spon_dir='D:\Brewer lab data\Yash\Filtered Spon MCD Files';

% stim_dir='D:\Brewer lab data\ttest files\Filtered Stim Data MCD Files';
% spon_dir='D:\Brewer lab data\ttest files\Filtered Spon MCD Files';

% stim_dir='D:\Brewer lab data\ttest HFS\Filtered Stim Data MCD Files';
% spon_dir='D:\Brewer lab data\ttest HFS\Filtered Spon MCD Files';

% stim_dir='D:\Brewer lab data\Theta Waves\Filtered Stim Data MCD Files';
% spon_dir='D:\Brewer lab data\Theta Waves\Filtered Spon MCD Files';

% spon_dir='D:\Brewer lab data\Ruiyi 5s validation\Filtered Spon MCD Files';
% stim_dir='D:\Brewer lab data\Ruiyi 5s validation\Filtered Stim Data MCD Files';

parent_dir='D:\Brewer lab data\HFS';
%stim_dir='D:\Brewer lab data\HFS\Filtered Stim Data MCD Files';
% stim_dir='D:\Brewer lab data\HFS\test\Filtered Stim Data MCD Files';
% spon_dir='D:\Brewer lab data\HFS\test\Filtered Spon MCD Files';
stim_dir='D:\Brewer lab data\HFS\No Stim';
% stim_dir='D:\Brewer lab data\HFS\HFS Stim';
% stim_dir='D:\Brewer lab data\HFS\Theta Stim';

%test dirs
% parent_dir='D:\Brewer lab data\HFS\Rep Dataset';
% spon_dir='D:\Brewer lab data\HFS\Rep Dataset\Filtered Spon MCD Files';
% stim_dir='D:\Brewer lab data\HFS\Rep Dataset\Filtered Stim Data MCD Files';

% parent_dir='D:\Brewer lab data\HFS\test 4';
% spon_dir='D:\Brewer lab data\HFS\test 4\Filtered Spon MCD Files';
% stim_dir='D:\Brewer lab data\HFS\test 4\Filtered Stim Data MCD Files';

% parent_dir='D:\Brewer lab data\test no stim';
% stim_dir='D:\Brewer lab data\test no stim\files';

low_thresh=5;
hi_thresh=500;

if_cw=[0,0,0,0,1,1,1,1,1];
% if_cw=[0,0,0,0,1,1];

%use for testing only:
%if_cw=[1];

%% Get sampling frequency and file length
cd(stim_dir)
filesHere_stim = dir('*.h5');
for i=1:length(filesHere_stim)
    stim_files(i)=McsHDF5.McsData(filesHere_stim(i).name);
    stim_fs(i)=double(1000000./stim_files(1, i).Recording{1, 1}.AnalogStream{1, 1}.Info.Tick(1));
    stim_duration(i)=stim_files(1, i).Recording{1, 1}.Duration/1000000; % In seconds
    stim_times_idx{i}=stim_files(1,i).Recording{1,1}.AnalogStream{1,1}.ChannelDataTimeStamps;
end

%scaling_factor=double(spon_duration(1)/stim_duration(1));
%scaling_factor=double(spon_duration./stim_duration);
%scaling_factor=5;

stim_scaling_factor=ones(1,length(stim_duration));
% for i=1:length(stim_files)
%     if stim_duration<=spon_duration(i)
%         stim_scaling_factor(i)=1;
%     else
%         stim_scaling_factor(i)=stim_duration(i)/spon_duration(i);
%     end
% end

%added 6/24/21
% spon_duration=double(spon_duration);
% stim_duration=double(stim_duration);

% spon_duration=[1,1,1,1,1];
% stim_duration=[1,1,1,1,1];

%Set max clusters to 4 for optimal NMI values
% Set to 200 just in case, large number
%max_clusters=200;
%% Notes:
% Switch back to old version, 210209 NMI is not working 2/09/11
% 210209 working now, don't know what happened 2/11/21

%% Calculate Areas
% Use tables instead of heirarchies, allows for vector indexing which gets
% rid of the for loops
%Wave_Clus_NMI_Func_211020
% Stimulated Area

tic
disp("Calculating stimulated areas...")
%[stim_areas,dated_output_folder_stim,NMI_Func_Version]=Wave_Clus_NMI_Func_210416(stim_dir,scaling_factor,stim_fs,max_clusters,3,1);
[stim_status_param,dated_output_folder_stim,NMI_Func_Version_stim,stim_times_dir,keep_stim]=...
    Wave_Clus_NMI_Func_211123(stim_dir,stim_scaling_factor,stim_fs,stim_duration,parent_dir,if_cw,2,1);
toc

% Interunit Peak Finding + Spike flagging no overlap
if length(stim_status_param)~=length(stim_status_param)
    error("Missing File in either Spon or Stim Folder")
end

total_bins=0;
overlap_bins=0;
total_area=0;
overlap_area=0;

for i=1:length(stim_status_param)
    %previous version 211027
    [areas,total_bins_calc,overlap_bins_calc,total_area_calc,overlap_area_clac]=...
        spike_flagging_single_211108(stim_status_param{i},stim_times_dir{i},keep_stim{i},stim_scaling_factor(i),parent_dir,if_cw(i));
    stim_areas{i}=areas{1};
    plot_spikes_from_NMI_211027(stim_areas{i}.channel_pair,stim_areas{i}.area,stim_areas{i}.up_times,stim_areas{i}.down_times,...
        strcat(dated_output_folder_stim,'\',stim_areas{i}.file_name))
    total_bins=total_bins+total_bins_calc;
    overlap_bins=overlap_bins+overlap_bins_calc;
    total_area=total_area+total_area_calc;
    overlap_area=overlap_area+overlap_area_clac;
end

%Where overlap bins is the number of bins where spikes need to be split
%into two axonal groups and the total number of bins are the number of bins
%that are considered to have their areas summed.

cd(dated_output_folder_stim)
save("stim_areas","stim_areas","NMI_Func_Version_stim","recorded_date")
cd(stim_dir)
save("stim_areas_latest","stim_areas","NMI_Func_Version_stim","recorded_date","dated_output_folder_stim")

%% All t-tests commented out 1/8/21: BE SURE TO UNCOMMENT AFTER RUN THROUGH OF ABOVE
%% T-test Setup
%load in data
% Check if stim data name contains spon data name
load('matching_table_template.mat')
cd(stim_dir)
load('stim_areas_latest.mat')
%%
channel_pairs=table2cell(matching_table(:,:));

%creates FID variables, toggle on and off as needed
%assign_fid_211118(dated_output_folder_stim,stim_dir,0)
cd(stim_dir)

% Adds FID
for i=1:length(stim_areas)
    fid_path=strcat(stim_dir,'/FID',string(i),'_Stim.mat');
    try
        load(fid_path)
        for j=1:length(filesHere_stim)
            if strcmp(stim_areas{j}.file_name,filename)
                stim_areas{j}.fid=fid;
            end
        end
    catch
        warning(strcat('The path ',fid_path,' does not exist.'))
        continue
    end
end

stim_file_names={filesHere_stim.name};
%% Reorganize areas files
% Set up data for t-test: Spon
% Reorganize loaded spon data

stim_areas=reorganize_NMI_outputs_210905(channel_pairs,stim_areas,if_cw,parent_dir);

has_up_stim=[1 1 1 1 1 1 1 1 1];

for i=1:length(stim_areas)
    
    if has_up_stim(i)
        stim_areas_table{i}=cell2table(stim_areas{i},'VariableNames',{'Subregion','Electrode Pairs','Spike Area','Conduction Time',...
        'Unit Pairs', 'Direction', 'Recording Name','fid','up_times','down_times'});
    else
        stim_areas_table{i}=cell2table(stim_areas{i},'VariableNames',{'Subregion','Electrode Pairs','Spike Area','Conduction Time',...
        'Unit Pairs', 'Direction', 'Recording Name','fid'});
    end
end

%% Order based on conduction time

for i=1:length(stim_areas_table)
    for j=1:length(stim_areas_table{i}.("Conduction Time"))
        if isempty(stim_areas_table{i}.("Conduction Time"){j})
            continue
        end        
        [sorted_times,sorted_indecies]=sort(stim_areas_table{i}.("Conduction Time"){j});
        stim_areas_table{i}.("Conduction Time"){j}=stim_areas_table{i}.("Conduction Time"){j}(sorted_indecies);
        %stim_areas_table{i}.("Unit Pairs"){j}=stim_areas_table{i}.("Unit Pairs"){j}(sorted_indecies);
        stim_areas_table{i}.("Direction"){j}=stim_areas_table{i}.("Direction"){j}(sorted_indecies);
        stim_areas_table{i}.("Spike Area"){j}=stim_areas_table{i}.("Spike Area"){j}(sorted_indecies);
        if strmatch('up_times',stim_areas_table{i}.Properties.VariableNames)~=0
            stim_areas_table{i}.("up_times"){j}=stim_areas_table{i}.("up_times"){j}(sorted_indecies);
            stim_areas_table{i}.("down_times"){j}=stim_areas_table{i}.("down_times"){j}(sorted_indecies);
        end
    end
end
%% Drop unit pairs that are less than 10% of max
% not reccomended for interunit spike finding

for i=1:length(stim_areas_table)
    for j=1:length(stim_areas_table{i}.("Spike Area"))
        if isempty(stim_areas_table{i}.("Spike Area"){j})
            continue
        end        
        for k=1:length(stim_areas_table{i}.("Spike Area"){j})
            [max_spike, max_spike_loc]=max(stim_areas_table{i}.("Spike Area"){j});
            ten_perc_max=0.1*stim_areas_table{i}.("Spike Area"){j}(max_spike_loc);
            stim_areas_table{i}.("Conduction Time"){j}=stim_areas_table{i}.("Conduction Time"){j}(stim_areas_table{i}.("Spike Area"){j}>ten_perc_max);
            stim_areas_table{i}.("Unit Pairs"){j}=stim_areas_table{i}.("Unit Pairs"){j}(stim_areas_table{i}.("Spike Area"){j}>ten_perc_max);
            %stim_areas_table{i}.("Unit Pairs"){j}=stim_areas_table{i}.("Unit Pairs"){j}(1:length(stim_areas_table{i}.("Conduction Time"){j}));
            stim_areas_table{i}.("Direction"){j}=stim_areas_table{i}.("Direction"){j}(stim_areas_table{i}.("Spike Area"){j}>ten_perc_max);
            stim_areas_table{i}.("Spike Area"){j}=stim_areas_table{i}.("Spike Area"){j}(stim_areas_table{i}.("Spike Area"){j}>ten_perc_max);
            if strmatch('up_times',stim_areas_table{i}.Properties.VariableNames)~=0
                stim_areas_table{i}.("up_times"){j}=stim_areas_table{i}.("up_times"){j}(stim_areas_table{i}.("Spike Area"){j}>ten_perc_max);
                stim_areas_table{i}.("down_times"){j}=stim_areas_table{i}.("down_times"){j}(stim_areas_table{i}.("Spike Area"){j}>ten_perc_max);
            end
        end            
    end
end

%% Split tables into feedforward and feedback
%Spontaneous Table

stim_areas_table_ff={};
stim_areas_table_fb={};
for i=1:length(stim_areas_table)
    if strmatch('up_times',stim_areas_table{i}.Properties.VariableNames)~=0
        stim_areas_table_ff{i}=table('Size',size(stim_areas_table{i}),'VariableTypes',{'string','string','cell','cell','cell','cell','cell','double','cell','cell'},'VariableNames',{'Subregion','Electrode Pairs','Spike Area','Conduction Time',...
        'Unit Pairs', 'Direction', 'Recording Name','fid','up_times','down_times'});
        stim_areas_table_fb{i}=table('Size',size(stim_areas_table{i}),'VariableTypes',{'string','string','cell','cell','cell','cell','cell','double','cell','cell'},'VariableNames',{'Subregion','Electrode Pairs','Spike Area','Conduction Time',...
        'Unit Pairs', 'Direction', 'Recording Name','fid','up_times','down_times'});
    else
        stim_areas_table_ff{i}=table('Size',size(stim_areas_table{i}),'VariableTypes',{'string','string','cell','cell','cell','cell','cell','double'},'VariableNames',{'Subregion','Electrode Pairs','Spike Area','Conduction Time',...
        'Unit Pairs', 'Direction', 'Recording Name','fid'});
        stim_areas_table_fb{i}=table('Size',size(stim_areas_table{i}),'VariableTypes',{'string','string','cell','cell','cell','cell','cell','double'},'VariableNames',{'Subregion','Electrode Pairs','Spike Area','Conduction Time',...
        'Unit Pairs', 'Direction', 'Recording Name','fid'});
    end
    for j=1:length(stim_areas_table{i}.("Direction"))
        if isempty(stim_areas_table{i}.("Direction"){j})
            stim_areas_table_ff{i}.("Subregion"){j}=stim_areas_table{i}.("Subregion"){j};
            stim_areas_table_fb{i}.("Subregion"){j}=stim_areas_table{i}.("Subregion"){j};
            stim_areas_table_ff{i}.("Electrode Pairs"){j}=stim_areas_table{i}.("Electrode Pairs"){j};
            stim_areas_table_fb{i}.("Electrode Pairs"){j}=stim_areas_table{i}.("Electrode Pairs"){j};
            stim_areas_table_ff{i}.("Spike Area"){j}=stim_areas_table{i}.("Spike Area"){j};
            stim_areas_table_fb{i}.("Spike Area"){j}=stim_areas_table{i}.("Spike Area"){j};
            stim_areas_table_ff{i}.("Conduction Time"){j}=stim_areas_table{i}.("Conduction Time"){j};
            stim_areas_table_fb{i}.("Conduction Time"){j}=stim_areas_table{i}.("Conduction Time"){j};
            stim_areas_table_ff{i}.("Unit Pairs"){j}=stim_areas_table{i}.("Unit Pairs"){j};
            stim_areas_table_fb{i}.("Unit Pairs"){j}=stim_areas_table{i}.("Unit Pairs"){j};
            stim_areas_table_ff{i}.("Direction"){j}=stim_areas_table{i}.("Direction"){j};
            stim_areas_table_fb{i}.("Direction"){j}=stim_areas_table{i}.("Direction"){j};
            stim_areas_table_ff{i}.("Recording Name"){j}=stim_areas_table{i}.("Recording Name"){j};
            stim_areas_table_fb{i}.("Recording Name"){j}=stim_areas_table{i}.("Recording Name"){j};
            stim_areas_table_ff{i}.("fid")(j)=stim_areas_table{i}.("fid")(j);
            stim_areas_table_fb{i}.("fid")(j)=stim_areas_table{i}.("fid")(j);
            
            %error found 8/16/21 left side set to spon areas table
            if strmatch('up_times',stim_areas_table{i}.Properties.VariableNames)~=0
                stim_areas_table_ff{i}.("up_times"){j}=stim_areas_table{i}.("up_times"){j};
                stim_areas_table_fb{i}.("up_times"){j}=stim_areas_table{i}.("up_times"){j};
                stim_areas_table_ff{i}.("down_times"){j}=stim_areas_table{i}.("down_times"){j};
                stim_areas_table_fb{i}.("down_times"){j}=stim_areas_table{i}.("down_times"){j};
            end
%             stim_areas_table_ff{i}(j,["Electrode Pairs",'Spike Area','Conduction Time','Unit Pairs','Direction','Recording Name'])=...
%                 stim_areas_table{i}(j,["Electrode Pairs",'Spike Area','Conduction Time','Unit Pairs','Direction','Recording Name']);
%             stim_areas_table_fb{i}(j,["Electrode Pairs",'Spike Area','Conduction Time','Unit Pairs','Direction','Recording Name'])=...
%                 stim_areas_table{i}(j,["Electrode Pairs",'Spike Area','Conduction Time','Unit Pairs','Direction','Recording Name']);
            
            continue
        end
        stim_areas_table_ff{i}.("Subregion"){j}=stim_areas_table{i}.("Subregion"){j};
        stim_areas_table_fb{i}.("Subregion"){j}=stim_areas_table{i}.("Subregion"){j};
        stim_areas_table_ff{i}.("Electrode Pairs"){j}=stim_areas_table{i}.("Electrode Pairs"){j};
        stim_areas_table_fb{i}.("Electrode Pairs"){j}=stim_areas_table{i}.("Electrode Pairs"){j};
        stim_areas_table_ff{i}.("fid")(j)=stim_areas_table{i}.("fid")(j);
        stim_areas_table_fb{i}.("fid")(j)=stim_areas_table{i}.("fid")(j);
        for k=1:length(stim_areas_table{i}.("Direction"){j})
            if stim_areas_table{i}.("Direction"){j}(k)==1
                stim_areas_table_ff{i}.("Spike Area"){j}=[stim_areas_table_ff{i}.("Spike Area"){j},stim_areas_table{i}.("Spike Area"){j}(k)];
                stim_areas_table_ff{i}.("Conduction Time"){j}=[stim_areas_table_ff{i}.("Conduction Time"){j}, stim_areas_table{i}.("Conduction Time"){j}(k)];
                stim_areas_table_ff{i}.("Unit Pairs"){j}=[stim_areas_table_ff{i}.("Unit Pairs"){j},stim_areas_table{i}.("Unit Pairs"){j}(k)];
                stim_areas_table_ff{i}.("Direction"){j}=[stim_areas_table_ff{i}.("Direction"){j}, stim_areas_table{i}.("Direction"){j}(k)];
                if strmatch('up_times',stim_areas_table{i}.Properties.VariableNames)~=0
                    stim_areas_table_ff{i}.("up_times"){j}=[stim_areas_table_ff{i}.("up_times"){j},stim_areas_table{i}.("up_times"){j}(k)];
                    stim_areas_table_ff{i}.("down_times"){j}=[stim_areas_table_ff{i}.("down_times"){j},stim_areas_table{i}.("down_times"){j}(k)];
                end
                
            elseif stim_areas_table{i}.("Direction"){j}(k)==0
                stim_areas_table_fb{i}.("Spike Area"){j}=[stim_areas_table_fb{i}.("Spike Area"){j},stim_areas_table{i}.("Spike Area"){j}(k)];
                stim_areas_table_fb{i}.("Conduction Time"){j}=[stim_areas_table_fb{i}.("Conduction Time"){j}, stim_areas_table{i}.("Conduction Time"){j}(k)];
                stim_areas_table_fb{i}.("Unit Pairs"){j}=[stim_areas_table_fb{i}.("Unit Pairs"){j}, stim_areas_table{i}.("Unit Pairs"){j}(k)];
                stim_areas_table_fb{i}.("Direction"){j}=[stim_areas_table_fb{i}.("Direction"){j}, stim_areas_table{i}.("Direction"){j}(k)];
                if strmatch('up_times',stim_areas_table{i}.Properties.VariableNames)~=0
                    stim_areas_table_fb{i}.("up_times"){j}=[stim_areas_table_fb{i}.("up_times"){j},stim_areas_table{i}.("up_times"){j}(k)];
                    stim_areas_table_fb{i}.("down_times"){j}=[stim_areas_table_fb{i}.("down_times"){j},stim_areas_table{i}.("down_times"){j}(k)];
                end
                
            end

        end
        stim_areas_table_ff{i}.("Recording Name"){j}=stim_areas_table{i}.("Recording Name"){j};
        stim_areas_table_fb{i}.("Recording Name"){j}=stim_areas_table{i}.("Recording Name"){j};
    end
end
%% Combine Units No Wobble
for i=1:length(stim_areas_table)
    stim_areas_table_ff{i}=combine_similar_times_nowobble_210719(stim_areas_table_ff{i});
    stim_areas_table_fb{i}=combine_similar_times_nowobble_210719(stim_areas_table_fb{i});
end

%% Combine Units with Wobble
% If there are equal numbers of units being compared across spon and stim,
% then do not combine 2/16/21 uncleared assumption
for i=1:length(stim_areas_table)
    stim_areas_table_ff{i}=combine_similar_times_210719(stim_areas_table_ff{i});
    stim_areas_table_fb{i}=combine_similar_times_210719(stim_areas_table_fb{i});
end

%% Network Stats

%stim
network_stats_211110(stim_areas_table_ff,stim_areas_table_fb)

%% Peak train assembly
[allregion_unit_matched_stim]=peak_train_assembly_210914(stim_areas_table_ff,stim_areas_table_fb,stim_times_idx);
[spike_burst_dyn_table_stim]=compute_spike_burst_dynamics_axons_220822(allregion_unit_matched_stim);
% save(strcat(parent_dir,'\','allregion_unit_matched_stim'),'allregion_unit_matched_stim','-v7.3')
% save(strcat(parent_dir,'\','spike_burst_dyn_table_stim'),'spike_burst_dyn_table_stim')
% save(strcat(dated_output_folder_stim,'\','allregion_unit_matched_stim'),'allregion_unit_matched_stim','-v7.3')
% save(strcat(dated_output_folder_stim,'\','spike_burst_dyn_table_stim'),'spike_burst_dyn_table_stim')