%% Create mea tables
% Creates tables compatible with
% compute_allregion_AFR_table_at_P1s_203011.m ultimately for creation of
% connectivity movie.
% In order to use, directories must be populated with matlab variables
% conaining wave_clus cluster_class variables. different files need to be
% created for feedforward and feedback axons. related file for genration is
% create_pseudo_times_files_220216

%% initialize
clear
clc

input_folder="D:\Brewer lab data\HFS\No Stim\xcorr_pseudo_times";
% input_folder="D:\Brewer lab data\HFS\Theta Stim\xcorr_pseudo_times";
% input_folder="D:\Brewer lab data\HFS\HFS Stim\xcorr_pseudo_times";
% output_folder="";
recordings_dir=dir(input_folder);
recordings_isdir=[recordings_dir.isdir];
recording_names=string({recordings_dir(recordings_isdir).name});
recording_names=recording_names(3:end);

fs=25000;
rec_length=300; %seconds

if_cw=[0,0,0,0,1,1,1,1,1];
% if_cw=[0,0,0,0,1,1];
%% Create tunnel mea
mea_subregion='tunnels';
for i=1:length(recording_names)
    current_folder=strcat(input_folder,'\',recording_names(i));
    sorted_mea_tunnels_ff=[];
    sorted_mea_tunnels_fb=[];
    par=[];
    par.fs=fs;
    par.subregion=mea_subregion;
    par.name=recording_names(i);
    if if_cw(i)
        par.orientation='cw';
        load(strcat(input_folder,"\cw_vars_tunnel.mat"))
    elseif ~if_cw(i)
        par.orientation='ccw';
        load(strcat(input_folder,"\ccw_vars_tunnel.mat"))
    else
        error("Please load directions correctly as a logical vector. -SBL")
    end
    
    %ff sorted_spike data
    times_files=dir(current_folder);
    times_names=string({times_files.name});
    times_names=times_names(3:end);
    sorted_spike_data_ff=[];
    sorted_spike_data_fb=[];
    for j=1:length(channel_names)
        idx_to_load=find(contains(times_names,"_ff") & contains(times_names,channel_names{j}));
        if isempty(idx_to_load)
            sorted_spike_data_ff{j,1}=[];
            sorted_spike_data_fb{j,1}=[];
            continue
        end
        file_to_load=times_names{idx_to_load};
        load(strcat(current_folder,'\',file_to_load))
        num_clus=max(unique(cluster_class(:,1),'stable'));
        spikes=[];
        count=[];
        for clus=1:num_clus
            vec=[1/fs:1/fs:rec_length]; 
            spike_times=cluster_class(cluster_class(:,1)==clus,2);
            logical_spike_vec=double(ismembertol(vec,spike_times./1e3,1e-10,'DataScale',1))';
            spikes{clus}=sparse(logical_spike_vec);
            count(clus)=sum((logical_spike_vec));
        end
        sorted_spike_data_ff{j,1}=spikes;
        sorted_spike_data_ff{j,2}=count;
    end
    %fb sorted_spike data
    for j=1:length(channel_names)
        idx_to_load=find(contains(times_names,"_fb") & contains(times_names,channel_names{j}));
        if isempty(idx_to_load)
            sorted_spike_data_fb{j,1}=[];
            sorted_spike_data_fb{j,1}=[];
            continue
        end
        file_to_load=times_names{idx_to_load};
        load(strcat(current_folder,'\',file_to_load))
        num_clus=max(unique(cluster_class(:,1),'stable'));
        spikes=[];
        count=[];
        for clus=1:num_clus
            vec=[1/fs:1/fs:rec_length]; 
            spike_times=cluster_class(cluster_class(:,1)==clus,2);
            logical_spike_vec=double(ismembertol(vec,spike_times./1e3,1e-10,'DataScale',1))';
            spikes{clus}=sparse(logical_spike_vec);
            count(clus)=sum((logical_spike_vec));
        end
        sorted_spike_data_fb{j,1}=spikes;
        sorted_spike_data_fb{j,2}=count;
    end
    
    sorted_mea_tunnels_ff.channel_names=channel_names;
    sorted_mea_tunnels_ff.par=par;
    sorted_mea_tunnels_ff.coordinates=coordinates;
    sorted_mea_tunnels_ff.subregion=subregions;
    sorted_mea_tunnels_ff.sorted_spike_data=sorted_spike_data_ff;
    
    sorted_mea_tunnels_fb.channel_names=channel_names;
    sorted_mea_tunnels_fb.par=par;
    sorted_mea_tunnels_fb.coordinates=coordinates;
    sorted_mea_tunnels_fb.subregion=subregions;
    sorted_mea_tunnels_fb.sorted_spike_data=sorted_spike_data_fb;
    
    mea=sorted_mea_tunnels_ff;
    save(strcat(current_folder,'\mea_tunnels_ff'),"mea")
    mea=sorted_mea_tunnels_fb;
    save(strcat(current_folder,'\mea_tunnels_fb'),"mea")
end

%% Create Subregion MEA

for i=1:length(recording_names)
    subregions=["EC","DG","CA3","CA1"];
    for regi=1:4
        current_folder=strcat(input_folder,'\',recording_names(i));
        sorted_mea_tunnels=[];
        par=[];
        par.fs=fs;
        %par.subregion=mea_subregion;
        par.name=recording_names(i);
        if if_cw(i)
            par.orientation='cw';
            load(strcat(input_folder,"\mea_",subregions(regi)),"cw")
            channel_names=cw.channel_names;
            coordinates=cw.coordinates;
            subregion=cw.subregion;
        elseif ~if_cw(i)
            par.orientation='ccw';
            load(strcat(input_folder,"\mea_",subregions(regi)),"ccw")
            channel_names=ccw.channel_names;
            coordinates=ccw.coordinates;
            subregion=ccw.subregion;
        else
            error("Please load directions correctly as a logical vector. -SBL")
        end

        %ff sorted_spike data
        times_files=dir(current_folder);
        times_names=string({times_files.name});
        times_names=times_names(3:end);
        times_names=times_names(contains(times_names,"times_") & contains(times_names,"_well"));
        for j=1:length(times_names)
            new_name=strsplit(times_names(j),{'_'});
            times_names(j)=new_name(2);
        end
        spike_data=[];
        for j=1:length(channel_names)
            idx_to_load=find(strcmp(times_names,channel_names{j}));
            if isempty(idx_to_load)
                spike_data{j,1}=[];
                spike_data{j,1}=[];
                continue
            end
            file_to_load="times_"+times_names{idx_to_load}+"_well.mat";
            load(strcat(current_folder,'\',file_to_load))
            %num_clus=max(unique(cluster_class(:,1),'stable'));
            spikes=[];
            count=[];
            for clus=1:1 %only one cluster for wells
                vec=[1/fs:1/fs:rec_length]; 
                spike_times=cluster_class(cluster_class(:,1)==clus,2);
                logical_spike_vec=double(ismembertol(vec,spike_times./1e3,1e-10,'DataScale',1))';
                spikes=sparse(logical_spike_vec);
                %count(clus)=sum((logical_spike_vec));
            end
            spike_data{j}=spikes;
            %spike_data{j,2}=count;
        end

        sorted_mea_tunnels.channel_names=channel_names;
        sorted_mea_tunnels.par=par;
        sorted_mea_tunnels.coordinates=coordinates;
        sorted_mea_tunnels.subregion=subregions;
        sorted_mea_tunnels.spike_data=spike_data;
        mea=sorted_mea_tunnels;
        save(strcat(current_folder,'\mea_',subregions(regi)),"mea")
    end
end