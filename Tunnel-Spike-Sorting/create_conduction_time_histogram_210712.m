function [bincountmatrix, channel_Out, xbin_ISI,keep_times]=  create_conduction_time_histogram_210712( folder_name, length_rec_sec)

%Edit made to original: added output channel name
channel_Out=[];

%Creates conduction time histogram for matched units present in
%matching_table_wNMI present in folder indicated by folder name
%
%Once unit-matching has been performed based on Normalized Matching Index, 
%this function creates conduction time histogram for the unit pairs
%indicated by unit pairs column in matching_table. The number of histograms
%for a tunnel is equal to number of rows in unit_pairs variable in the same
%row.
%Conduction time histogram gives 2 information. It gives the velocity of
%axonal communication captured by the tunnels and gives the direction of
%axonal communication.
%
%FORMAT: status=create_conduction_time_histogram_200701(folder_name)
%
%INPUT: folder_name contains address of the address of folder containing
%the output of wave_clus of the format times_xx, where xx-A6,A7 etc.
%e.g. - "C:\Users\yashv\Desktop\Yash\Ruiyi_data\op_folder\"
%
%OUTPUT: A file named 'matching_table_wNMI' is saved in the output folder

%check prerequisites

cd(char(folder_name))
fs=40e3;                                                                    %change samplig freq if needed
try
    load('matching_table_wNMI.mat')
catch
    error('Matching_table_wNMI.mat not found. Please run compute_NMI_200701 on the wave_clus output first');
end

if ~any(cellfun(@isempty, matching_table.unit_pairs))
    error('Please match unit_pairs manually')
end

no_channels = length(matching_table.tunnel_pairs);
if_flip = 0;if_hist = 1;

bincountmatrix = [];
keep_times=[];
% Loop through channels
% Ask about this
for chani=2:no_channels
    if isempty(matching_table.unit_pairs{chani})
        warning(char("No unit info found for "+matching_table.tunnel_pairs(chani)))
        continue
    end
    % Load tunnel channels
    split_channel_name = strsplit(matching_table.tunnel_pairs(chani),"-");
    downstream_channel_name = split_channel_name(2);
    upstream_channel_name = split_channel_name(1);
    
    % Warning if can't find files
    try
        down_var = load(char("times_"+downstream_channel_name), 'cluster_class');
    catch
        warning(char("Could not find "+downstream_channel_name))
        continue
    end

    try
        up_var = load(char("times_"+upstream_channel_name), 'cluster_class');
    catch
        warning(char("Could not find "+upstream_channel_name))
        continue
    end
   
    no_units = size(matching_table.unit_pairs{chani},1);   
    %Loop through unit-pairs
    for ui=1:no_units
        axon_unit_pair = matching_table.unit_pairs{chani}(ui,:);
        down_clusI = axon_unit_pair(1);
        up_clusI =  axon_unit_pair(2);
        %plot conduction time histogram
        down_channel = down_var.cluster_class(down_var.cluster_class(:,1) == down_clusI, 2);
        up_channel = up_var.cluster_class(up_var.cluster_class(:,1) == up_clusI, 2);
        %xbin_ISI should be the same every time, ok to leave it not
        %updating
        [bincountvec, ~, xbin_ISI,~,kept_times] = get_spike_delay_hist_210712(down_channel, up_channel, if_flip, if_hist, length_rec_sec);
        %saving
        bincountmatrix = [bincountmatrix;bincountvec];
        keep_times=[keep_times;kept_times];
        title(strcat(char(matching_table.tunnel_pairs(chani)+"-u"+ui),' up ',string(up_clusI),' down ',string(down_clusI)))
        saveas(gcf,char(matching_table.tunnel_pairs(chani)+"-u"+ui),'png')
        channel_Out=[channel_Out; {char(matching_table.tunnel_pairs(chani)+"-u"+ui)}];
    end
    disp(chani+" processed")
end
    
end

