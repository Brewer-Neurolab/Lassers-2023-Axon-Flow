%% Plots a raster plot of spikes for correlations

%% Load prerequisite data
clear
clc
load("D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\dataInfo.mat")
load("D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\full_idx_allregion_unit_matched_stim.mat")
data_folder_addr="D:\Brewer lab data\HFS\HFS stim\full_index_pseudo_times";
source=load("D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\AFR Output Full Index\concatenated_source_table.mat");
target=load("D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\AFR Output Full Index\\concatenated_target_table.mat");

data_folder_dir=dir(data_folder_addr);
data_folder_isdir=[data_folder_dir.isdir];
data_folder_names=string({data_folder_dir(data_folder_isdir).name});
data_folder_names=data_folder_names(3:end)';
all_region_order=[];
for i=1:length(data_folder_names)
    all_region_order(i)=find(contains(data_folder_names,dataInfo.meaName(i)));
end
allregion_unit_matched=allregion_unit_matched(all_region_order);
NodeTable = create_node_table_200820;
regList = ["EC","DG","CA3","CA1"];
source_regList_labels = ["EC->tunnels","DG->tunnels","CA3->tunnels","CA1->tunnels"];
target_regList_labels = ["tunnels->EC","tunnels->DG","tunnels->CA3","tunnels->CA1"];

%% Subregion to tunnels plot

slope_thresh=0.1;
rsq_thresh=0.2;
