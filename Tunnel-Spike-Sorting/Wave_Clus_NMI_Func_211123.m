function [status_params,dated_output_folder,NMI_version, folder_address, keep_times]=...
    Wave_Clus_NMI_Func_211123(directory,scalingFactor,fs,duration,parent_dir,if_cw,skip_to,reuse_last_NMI)
% Reuse last NMI only works if skip_to is not 0, reduces redundancy, saves memory 

%Change log starting 6/2/21
%added multiple scaling factor support
% Removed peak finding and output status and file names
% Changed primary output to include status, current tunnel pair, and ISI
% Added folder address output to function 7/7/21
%% Load and Convert h5 files into mat files
%areas=[];
status_params=[];
cd(directory) %Change this to file you want
filesHere = dir('*.h5');
keep_times=[];

%outputs NMI Version
st=dbstack;
NMI_version=st.name;

% Hard Code: Subject to change as needed
% These are the tunnel electrodes on an MEA 60
% tunnelElectrodes = {'A6';'A7';'B6';'B7';'C6';'C7';'D6';'D7';'F1';'G1';'F2';'G2';'F3';'G3';'F4';'G4'...
%     ;'F9';'G9';'F10';'G10';'F11';'G11';'F12';'G12';...
%     'J6';'J7';'K6';'K7';'L6';'L7';'M6';'M7'};

% These are the tunnel electrodes on an MEA 120
tunnelElectrodes = {'A6';'A7';'B6';'B7';'C6';'C7';'D6';'D7';'E6';'E7';...
    'F1';'G1';'F2';'G2';'F3';'G3';'F4';'G4';'F5';'G5';...
    'F8';'G8';'F9';'G9';'F10';'G10';'F11';'G11';'F12';'G12';...
    'H6';'H7';'J6';'J7';'K6';'K7';'L6';'L7';'M6';'M7'};

%test electrodes
% tunnelElectrodes = {'A6';'A7';'B6';'B7';...
%     'F4';'G4';'F5';'G5';...
%     'F11';'G11';'F12';'G12';...
%     'L6';'L7';'M6';'M7'};
%tunnelElectrodes = {'A6';'A7'};
% [filterParam, flag] = filter_getParam();


if skip_to~=0
    %Creates dated folders for saving analyses 
    current_date=date;
    date_files=dir;
    dirFolders=[date_files.isdir];
    subFolders=date_files(dirFolders);
    subFolders=string({subFolders.name});
    subFolders=subFolders(contains(subFolders,current_date));
    next_folder=strcat(current_date,'_',char(65+length(subFolders)));
    dated_output_folder=strcat(directory,'\',next_folder);
    
    for i=1:length(filesHere)
        
        % Creates directory and converts h5 files into mat files for current
        % specified directory 
        src_file_addr=strcat(filesHere(i).folder,'\',filesHere(i).name);
        folder_address{i} = [filesHere(i).folder,'\',next_folder,'\' ,filesHere(i).name(1:end-3), '_mat_files'];
        mkdir(folder_address{i})
        
        % If user has specified, this copies over all the times_ files from
        % the last time times_ files appeared in a subfolder of the same
        % filesHere(i).folder directory. This prevents having to rerun the
        % spike and clustering code
        if reuse_last_NMI==1
            subFolders=date_files(dirFolders);
            folder_dates=string({subFolders.date});
            time_diff=datenum(datetime)-datenum(folder_dates);
            sorted_diff=sort(time_diff);
            for diff=1:length(sorted_diff)
                subFolders_names=string({subFolders.name});
                last_folder=subFolders_names(time_diff==sorted_diff(diff));
                last_folder_fullpath=strcat(filesHere(i).folder,'\',last_folder,'\' ,filesHere(i).name(1:end-3),'_mat_files');
                dir_last_times=dir(strcat(last_folder_fullpath,'\','*.mat'));
                dir_last_times=string({dir_last_times.name});
                dir_last_times=dir_last_times(contains(dir_last_times,'times_'));
                %Changed 2/9/21 to check for if the previous folder
                %contains as many times files as there are electroeds in
                %the tunnel pairs we are examaning
                if ~isempty(dir_last_times)
%                     copyfile(strcat(last_folder_fullpath,'\','*.mat'),folder_address{i})
%                     copyfile(strcat(last_folder_fullpath,'\','*.txt'),folder_address{i})
                    copyfile(strcat(last_folder_fullpath,'\','times_*'),folder_address{i})
                    %Warns if there is not the max number of spike channels
                    %present in the folder with the max number of times
                    %files. This is a possibility if there was low spiking
                    %activity in that channel. However, the user could have
                    %stopped the program in the middle of running not
                    %generating most of the files
                    if length(dir_last_times)<length(tunnelElectrodes)
                        warning(strcat('Only ', string(length(dir_last_times)),...
                            ' times files found in last direcectory where times files exist.'))
                    end
                    break
                end                  
            end
            
        end
        
    end
        
else
    disp('Converting h5 to mat...')
    
    %Creates dated folders for saving analyses 
    current_date=date;
    date_files=dir;
    dirFolders=[date_files.isdir];
    subFolders=date_files(dirFolders);
    subFolders=string({subFolders.name});
    subFolders=subFolders(contains(subFolders,current_date));
    next_folder=strcat(current_date,'_',char(65+length(subFolders)));
    dated_output_folder=strcat(directory,'\',next_folder);
    
    for i=1:length(filesHere)

        % Creates directory and converts h5 files into mat files for current
        % specified directory 
        src_file_addr=strcat(filesHere(i).folder,'\',filesHere(i).name);
        folder_address{i} = [filesHere(i).folder,'\',next_folder,'\' ,filesHere(i).name(1:end-3), '_mat_files'];
        mkdir(folder_address{i})
        convert_h5_into_matfiles_210314(src_file_addr,folder_address{i},fs)
        
        %Creates txt for batch spikes
        fid = fopen('Files.txt','w+');
        for j=1:length(tunnelElectrodes)
            fprintf(fid, [tunnelElectrodes{j},'.mat']);
            if j~=length(tunnelElectrodes)
                fprintf(fid, '\n');
            end
        end
        fclose(fid);

        %Creates Spikes txt for batch clustering
        fid=fopen('Files_spikes.txt','w+');
        for j=1:length(tunnelElectrodes)
            fprintf(fid, [tunnelElectrodes{j},'_spikes','.mat']);
            if j~=length(tunnelElectrodes)
                fprintf(fid, '\n');
            end
        end

    %     %% Filters with SpyCode Functions and Down Samples
    %     
    %     outputMessage = filter_comput(src_file_addr,dst_file_addr{i},filterParam);
    %     if ~outputMessage
    %         errordlg('Error computing filtering!','!!Error!!')
    %         return
    %     end

    end
    cd ../
end

% filter tunnel data
% for i =1:length(folder_address)
%     for j=1:length(tunnelElectrodes)
%         load(strcat(folder_address{i},'\',tunnelElectrodes{j},'.mat'),'data','sr')
%         data=highpass(data,400,sr);
%         save(strcat(folder_address{i},'\',tunnelElectrodes{j},'.mat'),'data')
%     end
% end

%% Get Spikes and do clustering
% Note to self: need sampling rates from h5 file
% Set current to 40khz
% Use the GUI for this section!!!!!!!!!!!!!!!!!
pause on
if skip_to==1||skip_to==0
    hi_lo_wave_clus_batch_211020(folder_address,tunnelElectrodes,duration)
end

%% Create Conduction Time Histogram and Fitting Spontaneous peak to stimulation histogram and get stim values(Updated)

% Use the new "create conduction time histogram'' script, make sure to use 
% the new function, the 'bincountmatrix' will give information of bincount 
% numbers (as row vectors' for each histogram. And the record time is a new 
% variable to calculate the averaged response: eg. the initial recording 
% length is 300 s, and we want to average to 5 s and compare with another 5 s 
% recording, then the factor is 300/5 = 60, you need to put 60 at the second 
% input variable. It will generate both new histogram and bincount info.

%Use the script 'Sponpeakarea' to fit each row vector to find the primary 
% peak and the area underneath it. The purpose is to determine a single event 
% and its values conducted in the tunnel. Notice: 1.you will need to define 
% the direction of the primary peak, if the primary peak is at positive 
% conduction time, the input "z =1", on the other hand, negative 
% conduction time, use "z=0". 2. Match the single row vector with the 
% right name of the tunnel. Eg: the first vector may match with "F1-G1'', 
% but the layout is different, please rename the stored pics with right 
% channel info.
if skip_to>=1||skip_to==0
    for i=1:length(folder_address)
        cd(folder_address{i})
        %cluster_frequency_rejection_211021(folder_address{i},duration)
        combine_all_clusters_211026(folder_address{i})
        fprintf('Calculating NMI...\n')
        compute_NMI_210726(folder_address{i},parent_dir,if_cw(i))
        
        load('matching_table_wNMI.mat')
        
        % Loads matching table (assumes current folder is correct directory to load
        % from) and deletes all png files in the current directory containing the
        % tunnel pairs
        fileshere=dir('*.png');
        size_matching_table=size(matching_table);
        for matching_table_idx=1:size_matching_table
            for files_idx=1:length(fileshere)
                if contains(fileshere(files_idx).name,matching_table.tunnel_pairs(matching_table_idx))
                    delete(fileshere(files_idx).name)
                end
            end
        end
        
        sizeMatchingTable=size(matching_table);
        
        % Labels each areas structure with the name of the folder the data
        % was from
        last_slash=strfind(folder_address{i},'\');
        folder_name=folder_address{i}(last_slash(end)+1:end);
        status_params{i}.file_name=folder_name;
        
        for j=1:sizeMatchingTable(1)
            current_Tunnel_Pair = matching_table.tunnel_pairs{j};
            if isnan(matching_table.NMI{j})
                continue
            end
            NMI_mat=matching_table.NMI{j};
%             unit_match_mat=match_units_to_find_axons_203012(NMI_mat,.2);
            unit_match_mat=match_units_to_find_axons_211019(NMI_mat,.2);
%             [maxNMIVal, maxNMIIdx]=max(NMI_mat);
%             Unit_pairs=[];
%             for k=1:length(maxNMIIdx)
%                 Unit_pairs=[Unit_pairs;[k,maxNMIIdx(k)]];
%             end
             matching_table.unit_pairs{j}=unit_match_mat;
        end
        save('matching_table_wNMI.mat', 'matching_table')
        %Below is what Ruiyi recommended 
        %Make sure the bins are some multiple of 25us because of 40kHz sampling
        %rate
        % the second input is a scaling parameter, number of spike pairs per
        % second
        % The following function is a wrapper function, the meat is later on
        % and if I i have to do anything in the future, it might be pertinent
        % to use that instead, use tables, not nested cells :)
        %status=create_conduction_time_histogram_200721(folder_address{i},60);%IMPORTANT NOTE: the 60 is hardcoded in because each recording is 300s and there are 5 recordings, 300/5=60
        fprintf('Calculating Conduction Time Histogram...\n')
        [status, current_Tunnel_Pair, ISI,kept_times]=create_conduction_time_histogram_210712(folder_address{i},scalingFactor(i));
        status_params{i}.status=status;
        status_params{i}.current_Tunnel_Pair=current_Tunnel_Pair;
        status_params{i}.ISI=ISI;
        keep_times{i}=kept_times;
    end
end
end