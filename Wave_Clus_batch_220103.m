function [dated_output_folder,NMI_version, folder_address, keep_times]=...
    Wave_Clus_batch_220103(directory,well_dir,fs,skip_to,reuse_last_NMI)
% Reuse last NMI only works if skip_to is not 0, reduces redundancy, saves memory 

%Change log starting 6/2/21
%added multiple scaling factor support
% Removed peak finding and output status and file names
% Changed primary output to include status, current tunnel pair, and ISI
% Added folder address output to function 7/7/21
%% Load and Convert h5 files into mat files
%areas=[];
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
    dated_output_folder=strcat(well_dir,'\',next_folder);
    
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
    dated_output_folder=strcat(well_dir,'\',next_folder);
    
    for i=1:length(filesHere)

        % Creates directory and converts h5 files into mat files for current
        % specified directory 
        src_file_addr=strcat(filesHere(i).folder,'\',filesHere(i).name);
        folder_address{i} = [filesHere(i).folder,'\','Wells','\',next_folder,'\' ,filesHere(i).name(1:end-3), '_mat_files'];
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

%% Get Spikes
% Note to self: need sampling rates from h5 file
% Set current to 40khz
% Use the GUI for this section!!!!!!!!!!!!!!!!!
if skip_to==1||skip_to==0
    cd(well_dir)
    current_dir=dir;
    recordings=current_dir(current_dir.isdir);
    recordings=recordings(3:end); %removes relevant pathing
    recording_names=string({recordings.name});
    for i=1:length(recording_names)
        cd(strcat(well_dir,recording_names(i)))
        current_dir=dir;
        mat_files=current_dir(contains(string({current_dir.name}),".mat"));
        mat_files=string({mat_files.name});
        for j=1:length(mat_files)
            Get_spikes(char(mat_files(j)))
        end
    end
end

end