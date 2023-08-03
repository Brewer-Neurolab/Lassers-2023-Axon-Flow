function compute_NMI_210726(folder_name, parent_dir, if_cw)
%Computes Normalized Matching Index matrix(NMI) for all output of wave_clus 
%tunnel pairs in the folder indicated by folder_name
%
%Normalized Matching Index(NMI) is used to indentify which using in
%tunnel-1 corresponds to unit in tunnel-2 to form an axon, where units are 
%created using spike sorting with wave_clus. It quantifies
%the similarity in spike-timings between the 2 tunnel channels, similarity
%defined as no. of spikes within a time window (based on physiological
%estimate of velocity).
%
%This function requires 'matching_table_template to be in the workspace. It
%creates a matrix of NMI with size = [no. of units in tunnel-1, no. of
%units in tunnel-2]. It creates a new column in the matching_table called
%NMI where it stores NMI matrix for every col.
%
%Range of NMI = [0:1], If NaN, check if the appropriate channel wiles are
%present in the folder
%
%NMI = No. of paired spikes / max(no. of spikes in upstream tunnel, no. of
%                               spikes in the downstream tunnel)
%
%FORMAT: compute_NMI_200701(folder_name)
%
%INPUT: folder_name contains address of the address of folder containing
%the output of wave_clus of the format times_xx, where xx-A6,A7 etc.
%e.g. - "C:\Users\yashv\Desktop\Yash\Ruiyi_data\op_folder\"
%
%OUTPUT: A file named 'matching_table_wNMI' is saved in the output folder

% Loading prerequisites
fs=40e3;                                                                    %change samplig freq if needed
try 
    if if_cw
        load(strcat(parent_dir,'\','matching_table_cw'))
    elseif ~if_cw
        load(strcat(parent_dir,'\','matching_table_ccw'))
    else
        error("Please define if the array is cw or ccw with 1 and 0.")
    end
catch
    error('Could not find matching_table_template. Please include it in the path')
end

cd(char(folder_name))

no_channels = length(matching_table.tunnel_pairs);
NMI = cell(no_channels,1);
unit_pairs = cell(no_channels,1);

% Loop through channels
for chani=1:no_channels
    % Load tunnel channels
    split_channel_name = strsplit(matching_table.tunnel_pairs(chani),"-");
    
        downstream_channel_name = split_channel_name(2);
        upstream_channel_name = split_channel_name(1);
    % Warning if can't find files
    try
        down_var = load(char("times_"+downstream_channel_name), 'cluster_class');
    catch
        warning(char("Could not find "+downstream_channel_name))
        NMI{chani} = NaN;
        continue
    end

    try
        up_var = load(char("times_"+upstream_channel_name), 'cluster_class');
    catch
        warning(char("Could not find "+upstream_channel_name))
        NMI{chani} = NaN;
        continue
    end
    % Compute NMI_mat
    no_clus_down = max(down_var.cluster_class(:,1));
    no_clus_up = max(up_var.cluster_class(:,1));
    if_flip = 0;if_hist = 1; %Changed 11/10/20 from 0->1
    NMI_mat = zeros(no_clus_down, no_clus_up);
    
    for down_clusI = 1:no_clus_down
        for up_clusI = 1:no_clus_up
            down_channel = down_var.cluster_class(down_var.cluster_class(:,1) == down_clusI, 2);
            up_channel = up_var.cluster_class(up_var.cluster_class(:,1) == up_clusI, 2);
            [~, NMI_val] = get_spike_delay_hist_200721(down_channel, up_channel, if_flip, if_hist, 1);
            NMI_mat(down_clusI, up_clusI) = NMI_val;
        end
    end
    % Store in table
    NMI{chani} = NMI_mat;
    disp(chani+" processed")
end

% Append NMI and unit_pairs
matching_table = [matching_table, table(NMI, unit_pairs)];

% save in the folder
save('./matching_table_wNMI','matching_table')
end

