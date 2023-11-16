%% Plot from Allregions

%% Initialization

clear
clc

subregion="";
interRegion="";

%no stim
raw_dir="D:\Brewer lab data\HFS\No Stim\23-Nov-2021_B";
tunnel_allregion=load("D:\Brewer lab data\HFS\No Stim\23-Nov-2021_B\allregion_unit_matched_stim.mat");
well_allregion=load("D:\Brewer lab data\HFS\No Stim\Wells_5SD_500maxSD\allregion_unit_matched.mat");
load("D:\Brewer lab data\HFS\Temporal Analysis\No Stim\dataInfo.mat")

data_folder_dir=dir(raw_dir);
data_folder_isdir=[data_folder_dir.isdir];
data_folder_names=string({data_folder_dir(data_folder_isdir).name});
data_folder_names=data_folder_names(3:end)';
% data_folder_names=erase(data_folder_names,"_mat_files")';
all_region_order=[];
for i=1:length(data_folder_names)
    all_region_order(i)=find(contains(data_folder_names,dataInfo.meaName(i)));
end
tunnel_allregion=tunnel_allregion.allregion_unit_matched_stim(all_region_order);
well_allregion=well_allregion.all_regions_unit_matched(all_region_order);

fs=25000;

time_idx=[1/fs:1/fs:300];
%% Plotting

for fi=1:9
    dataInfo_idx=dataInfo.s_no==fi;
    f=figure('units','normalized','OuterPosition',[0 0 1 1]);
    tunnels=tunnel_allregion{fi}(tunnel_allregion{fi}.Subregion==interRegion);%...
        %& (~isempty(tunnel_allregion{fi}.up_ff) | ~isempty(tunnel_allregion{fi}.up_fb)));
    wells=well_allregion{fi}(well_allregion{fi}.Subregion==subregion);
    
    t=tiledlayout('flow','TileSpacing','tight','Padding','tight');
    maxTiles=1000;
    
    for regi=1:num_tunnels
        ax=nexttile([1,maxTiles]);
        axes=[axes,ax];
        hold on
        load(strcat(raw_dir,'\',dataInfo.meaName(dataInfo_idx),'\',tunnels(tun),".mat"))
        hold off
    end
    
    for regi=1:num_wells
        
    end
end