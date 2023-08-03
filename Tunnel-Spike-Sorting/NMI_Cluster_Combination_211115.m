function NMI_Cluster_Combination_211115(folder_address,if_cw,parent_dir,length_rec_sec)

cd(folder_address)
times_dir=dir('times_*.mat');
times_name=string({times_dir.name});

fprintf('Calculating NMI...\n')
compute_NMI_210726(folder_address,parent_dir,if_cw)

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
last_slash=strfind(folder_address,'\');
folder_name=folder_address(last_slash(end)+1:end);
status_params.file_name=folder_name;

for j=1:sizeMatchingTable(1)
    current_Tunnel_Pair = matching_table.tunnel_pairs{j};
    if isnan(matching_table.NMI{j})
        continue
    end
    NMI_mat=matching_table.NMI{j};
%             unit_match_mat=match_units_to_find_axons_203012(NMI_mat,.2);
    unit_match_mat=match_units_to_find_axons_211019(NMI_mat,.05);
%             [maxNMIVal, maxNMIIdx]=max(NMI_mat);
%             Unit_pairs=[];
%             for k=1:length(maxNMIIdx)
%                 Unit_pairs=[Unit_pairs;[k,maxNMIIdx(k)]];
%             end
     matching_table.unit_pairs{j}=unit_match_mat;
end
save('matching_table_wNMI.mat', 'matching_table')

tunnel_merge_211124(folder_name,length_rec_sec)

end