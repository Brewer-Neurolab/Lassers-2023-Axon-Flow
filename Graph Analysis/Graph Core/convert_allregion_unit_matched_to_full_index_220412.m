%% Two to One Electrode tunnels
% To be used with the AFR Computation scripts
% This program transplants the whole spike index over two electrodes if
% there is only one axon. This is to cut down on errors evolved from poor
% coupling of axon to electrode. Current allregion_unit_matched files have
% outputs from NMI, load original spike index before NMI

%% Load Prerequisites

clear
clc

% 11SD
% parent_dir="D:\Brewer lab data\HFS\No Stim\23-Nov-2021_A";
% pseudo_dir="D:\Brewer lab data\HFS\No Stim\xcorr_pseudo_times";
% load('D:\Brewer lab data\HFS\No Stim\23-Nov-2021_A\allregion_unit_matched_stim.mat')
% output_dir="D:\Brewer lab data\HFS\Temporal Analysis\No Stim\";

% parent_dir="D:\Brewer lab data\HFS\Theta Stim\03-Mar-2022_B";
% parent_dir="D:\Brewer lab data\HFS\Theta Stim\10-May-2022_A";
% pseudo_dir="D:\Brewer lab data\HFS\Theta Stim\xcorr_pseudo_times";
% load(parent_dir+'\allregion_unit_matched_stim.mat');
% output_dir="D:\Brewer lab data\HFS\Temporal Analysis\5 HFS\"; %40 HFS previously was a mistake

% parent_dir="D:\Brewer lab data\HFS\HFS Stim\24-Nov-2021_A";
% pseudo_dir="D:\Brewer lab data\HFS\HFS Stim\xcorr_pseudo_times";
% load(parent_dir+'\allregion_unit_matched_stim.mat')
% output_dir="D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\"; %5 HFS previously was a mistake

% 5SD
% parent_dir="D:\Brewer lab data\HFS\No Stim\02-Oct-2022_A";
% pseudo_dir="D:\Brewer lab data\HFS\No Stim\xcorr_pseudo_times";
% load('D:\Brewer lab data\HFS\No Stim\02-Oct-2022_A\allregion_unit_matched_stim.mat')
% output_dir="D:\Brewer lab data\HFS\Temporal Analysis\No Stim\";

% parent_dir="D:\Brewer lab data\HFS\Theta Stim\27-Sep-2022_A";
% pseudo_dir="D:\Brewer lab data\HFS\Theta Stim\xcorr_pseudo_times";
% load(parent_dir+'\allregion_unit_matched_stim.mat');
% output_dir="D:\Brewer lab data\HFS\Temporal Analysis\5 HFS\"; 

parent_dir="D:\Brewer lab data\HFS\HFS Stim\03-Oct-2022_A";
pseudo_dir="D:\Brewer lab data\HFS\HFS Stim\xcorr_pseudo_times";
load(parent_dir+'\allregion_unit_matched_stim.mat')
output_dir="D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\"; 

fs=25000;

time_idx=[1/fs:1/fs:300];

% get subfolders
times_dir=dir(parent_dir);
times_folders=string({times_dir([times_dir.isdir]).name});
times_folders=times_folders(3:end);

% copy subfolders to new dir for modification
top_dir=erase(pseudo_dir,"xcorr_pseudo_times");
if exist(strcat(top_dir,"full_index_pseudo_times"),'dir')
    rmdir(strcat(top_dir,"full_index_pseudo_times"),'s')
    mkdir(strcat(top_dir,"full_index_pseudo_times"))
    for i=1:length(times_folders)
        mkdir(strcat(top_dir,"full_index_pseudo_times\",times_folders(i)))
        
        filenames=dir(pseudo_dir);
        filenames=filenames(~[filenames.isdir]);
        for j=1:length(filenames)
            copyfile(pseudo_dir+'\'+filenames(j).name,strcat(top_dir,"full_index_pseudo_times"))
        end
        
        filenames=dir(pseudo_dir+'\'+times_folders(i));
        for j=3:length(filenames)
            copyfile(pseudo_dir+'\'+times_folders(i)+'\'+filenames(j).name,strcat(top_dir,"full_index_pseudo_times\",times_folders(i)))
        end
    end
else
    mkdir(strcat(top_dir,"full_index_pseudo_times"))
    for i=1:length(times_folders)
        mkdir(strcat(top_dir,"full_index_pseudo_times\",times_folders(i)))
        
        filenames=dir(pseudo_dir);
        filenames=filenames(~[filenames.isdir]);
        for j=1:length(filenames)
            copyfile(pseudo_dir+'\'+filenames(j).name,strcat(top_dir,"full_index_pseudo_times"))
        end
        
        filenames=dir(pseudo_dir+'\'+times_folders(i));
        for j=3:length(filenames)
            copyfile(pseudo_dir+'\'+times_folders(i)+'\'+filenames(j).name,strcat(top_dir,"full_index_pseudo_times\",times_folders(i)),'f')
        end
    end
end

full_idx_dir=strcat(top_dir,"full_index_pseudo_times");

%NodeTable = create_node_table_200820;
regList = ["EC","DG","CA3","CA1"];
source_regList_labels = ["EC->tunnels","DG->tunnels","CA3->tunnels","CA1->tunnels"];
target_regList_labels = ["tunnels->EC","tunnels->DG","tunnels->CA3","tunnels->CA1"];

if_cw=[0,0,0,0,1,1,1,1,1];
%%
for i=1:length(allregion_unit_matched_stim)
    for j=1:height(allregion_unit_matched_stim{i})
        %get electrode names
        electrodes=strsplit(allregion_unit_matched_stim{i}.("Electrode Pairs"){j},{'-'});
        % get num axons
        num_axons=length([allregion_unit_matched_stim{i}.up_ff{j},allregion_unit_matched_stim{i}.up_fb{j}]);
        
        %load associated files
        try
            elec1=load(parent_dir+'\'+times_folders(i)+'\'+'times_'+electrodes(1)+'.mat');
            elec1_idx=elec1.cluster_class(:,2);
        catch
            elec1_idx=[];
        end
        
        try
            elec2=load(parent_dir+'\'+times_folders(i)+'\'+'times_'+electrodes(2)+'.mat');
            elec2_idx=elec2.cluster_class(:,2);
        catch
            elec2_idx=[];
        end
        
        % replace spikes with full spike index if num axons is 1
        if num_axons==1
            if length(elec1_idx)>=length(elec2_idx)
                %construct full index
                full_index=ismembertol(time_idx,elec1_idx/1000,1e-10,'DataScale',1);
                allregion_unit_matched_stim{i}.OriginalChan{j}=electrodes{1};
            else
                %construct full index
                full_index=ismembertol(time_idx,elec2_idx/1000,1e-10,'DataScale',1);
                allregion_unit_matched_stim{i}.OriginalChan{j}=electrodes{2};
            end
            if ~isempty(allregion_unit_matched_stim{i}.up_ff{j}) && isempty(allregion_unit_matched_stim{i}.up_fb{j})
                allregion_unit_matched_stim{i}.up_ff{j}={full_index};
                allregion_unit_matched_stim{i}.down_ff{j}={full_index};
                write_new_tunnel_pseudo_times(full_idx_dir+'\'+times_folders(i), full_index, time_idx, 1, electrodes)
            elseif isempty(allregion_unit_matched_stim{i}.up_ff{j}) && ~isempty(allregion_unit_matched_stim{i}.up_fb{j})
                allregion_unit_matched_stim{i}.up_fb{j}={full_index};
                allregion_unit_matched_stim{i}.down_fb{j}={full_index};
                write_new_tunnel_pseudo_times(full_idx_dir+'\'+times_folders(i), full_index, time_idx, 0, electrodes)
            else
                error("More than one axon.-SBL")
            end
        % if there is no axon identified yet still spikes in tunnels,
        % correlate active channel with wells on either side of tunnels.
        % A majorty of the tunnels are unidirectional, will assume all are
        % for now to decide which column the spike indecies are placed in.
        elseif num_axons==0
            active_elec=[];
            if length(elec1_idx)>length(elec2_idx)
                active_elec=elec1_idx;
                active_elec_name=electrodes{1};
                allregion_unit_matched_stim{i}.OriginalChan{j}=active_elec_name;
            elseif length(elec1_idx)<length(elec2_idx)
                active_elec=elec2_idx;
                active_elec_name=electrodes{2};
                allregion_unit_matched_stim{i}.OriginalChan{j}=active_elec_name;
            else
                continue
            end
            
            regions=strsplit(allregion_unit_matched_stim{i}.Subregion{j},{'-'});
            reg1=load(pseudo_dir+'\mea_'+regions(1)+'.mat');
            reg2=load(pseudo_dir+'\mea_'+regions(2)+'.mat');
            
            if if_cw(i)
                reg1=reg1.cw;
                reg2=reg2.cw;
            else
                reg1=reg1.ccw;
                reg2=reg2.ccw;
            end
            regs=[reg1;reg2];
            vote_dir=[];
            vote_r=[];
            for reg=1:2
                for chani=1:length(regs(reg).channel_names)
                    load(pseudo_dir+'\'+times_folders(i)+'\'+'times_'+string(regs(reg).channel_names(chani))+'_well.mat')
                    well_elec=cluster_class(:,2);
                    if reg==1
                        [direction,r,if_save]=seed_tunnels_with_direction_220321(well_elec,active_elec,fs);
                    else
                        [direction,r,if_save]=seed_tunnels_with_direction_220321(active_elec,well_elec,fs);
                    end
                    title(["Well: "+string(regs(reg).channel_names(chani)),"Tunnel: "+string(active_elec_name)])
                    
                    if ~exist(full_idx_dir+"\"+times_folders(i)+"\"+"xcorr_pics",'dir')
                        mkdir(full_idx_dir+"\"+times_folders(i)+"\"+"xcorr_pics")
                    end
                    
                    if if_save
                    saveas(gcf,full_idx_dir+'\'+times_folders(i)+'\'+'xcorr_pics'+...
                        '\'+string(regs(reg).channel_names(chani))+'_'+active_elec_name+'.png')
                    end
                    vote_dir=[vote_dir, direction]; 
                    vote_r=[vote_r,r];
                    close(gcf)
                end
            end
            if isempty(vote_dir)
                continue
            end
%             num_ff=length(find(vote_dir==1));
%             num_fb=length(find(vote_dir==0));
            idx_ff=(vote_dir==1);
            idx_fb=(vote_dir==0);
            ff_r=sum(vote_r(idx_ff));
            fb_r=sum(vote_r(idx_fb));
            
            if ff_r>fb_r
                full_index=ismembertol(time_idx,active_elec/1000,1e-10,'DataScale',1);
                allregion_unit_matched_stim{i}.up_ff{j}={full_index};
                allregion_unit_matched_stim{i}.down_ff{j}={full_index};
                write_new_tunnel_pseudo_times(full_idx_dir+'\'+times_folders(i), full_index, time_idx, 1, electrodes)
            elseif ff_r<fb_r
                full_index=ismembertol(time_idx,active_elec/1000,1e-10,'DataScale',1);
                allregion_unit_matched_stim{i}.up_fb{j}={full_index};
                allregion_unit_matched_stim{i}.down_fb{j}={full_index};
                write_new_tunnel_pseudo_times(full_idx_dir+'\'+times_folders(i), full_index, time_idx, 0, electrodes)
            end
        end
    end
    disp(string(i)+" file(s) completed.")
end
% 
% for i=1:length(allregion_unit_matched_stim)
%     for j=1:length(regList)
%         load(info_dir+'\mea_'+regList(j)+'.mat')
%         if if_cw(i)
%             
%         else
%             
%         end
%     end
% end

save(output_dir+"full_idx_allregion_unit_matched_stim.mat", "allregion_unit_matched_stim",'-v7.3')