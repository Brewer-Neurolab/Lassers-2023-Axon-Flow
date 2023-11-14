%% Plots raw data for tunnels and wells using the full spiking index created for the AFR/edge plotting scripts

%% Prerequisites
clear 
clc

tunnel_reg="CA3-CA1";
well_reg="";

wave_color=load("D:\Brewer lab data\HFS\wave_clus_colororder.mat");
wave_color=wave_color.wave_clus_color_order;

%no stim
% save_to="D:\Brewer lab data\HFS\Temporal Analysis\No Stim\raw data\";
% raw_dir="D:\Brewer lab data\HFS\No Stim\23-Nov-2021_B";
% spikes_dir="D:\Brewer lab data\HFS\No Stim\full_index_pseudo_times";
% tunnel_allregion=load("D:\Brewer lab data\HFS\Temporal Analysis\No Stim\full_idx_allregion_unit_matched_stim.mat");
% well_allregion=load("D:\Brewer lab data\HFS\No Stim\Wells_5SD_500maxSD\allregion_unit_matched.mat");
% load("D:\Brewer lab data\HFS\Temporal Analysis\No Stim\dataInfo.mat")

%5 hfs
save_to="D:\Brewer lab data\HFS\Temporal Analysis\5 HFS\raw data\";
raw_dir="D:\Brewer lab data\HFS\Theta Stim\03-Mar-2022_A";
spikes_dir="D:\Brewer lab data\HFS\Theta Stim\full_index_pseudo_times";
tunnel_allregion=load("D:\Brewer lab data\HFS\Temporal Analysis\5 HFS\full_idx_allregion_unit_matched_stim.mat");
well_allregion=load("D:\Brewer lab data\HFS\Theta Stim\Wells_5SD_500maxSD\allregion_unit_matched.mat");
load("D:\Brewer lab data\HFS\Temporal Analysis\5 HFS\dataInfo.mat")

%40 HFS
% save_to="D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\raw data\";
% raw_dir="D:\Brewer lab data\HFS\HFS Stim\24-Nov-2021_A";
% spikes_dir="D:\Brewer lab data\HFS\HFS Stim\full_index_pseudo_times";
% tunnel_allregion=load("D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\full_idx_allregion_unit_matched_stim.mat");
% well_allregion=load("D:\Brewer lab data\HFS\HFS Stim\Wells_5SD_500maxSD\allregion_unit_matched.mat");
% load("D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\dataInfo.mat")

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

for fi=6%:length(well_allregion)
    
    dataInfo_idx=dataInfo.s_no==fi;
    
    f=figure('units','normalized','OuterPosition',[0 0 1 1]);
    num_tunnels=sum(tunnel_allregion{dataInfo_idx}.Subregion==tunnel_reg & ...
        ~isempty(cellfun(@isempty,tunnel_allregion{dataInfo_idx}.up_ff)) & ...
        ~isempty(cellfun(@isempty,tunnel_allregion{dataInfo_idx}.up_fb)));
    is_pop = @(spike_array) sum(spike_array)>100;
    num_well=sum(well_allregion{dataInfo_idx}.Subregion==well_reg & ...
        cellfun(is_pop,well_allregion{dataInfo_idx}.("Spike Train")));
    
%     t=tiledlayout(num_tunnels+num_well,1,'TileSpacing','tight','Padding','tight');
    t=tiledlayout('flow','TileSpacing','tight','Padding','tight');
    maxTiles=1000;

    single_axon_tunnels=string(tunnel_allregion{dataInfo_idx}.OriginalChan...
        (tunnel_allregion{dataInfo_idx}.Subregion==tunnel_reg & ...
        ~cellfun(@isempty,tunnel_allregion{dataInfo_idx}.OriginalChan)));
    multi_axon_tunnels=tunnel_allregion{dataInfo_idx}.("Electrode Pairs")(cellfun(@isempty,tunnel_allregion{dataInfo_idx}.OriginalChan) ...
        & tunnel_allregion{dataInfo_idx}.Subregion==tunnel_reg);
    splitter=@(x) strsplit(x,{'-'});
    multi_axon_tunnels=(cellfun(splitter,multi_axon_tunnels,'UniformOutput',false));
    multi_axon_tunnels=string(reshape([multi_axon_tunnels{:}],[],2))';
    if ~isempty(multi_axon_tunnels)
        multi_axon_tunnels=multi_axon_tunnels(:,1);
        tunnels=[single_axon_tunnels;multi_axon_tunnels];
    else
        tunnels=[single_axon_tunnels];
    end
    
    axes=[];
    
    
    [~,idx]=sort_nat(tunnel_allregion{dataInfo_idx}.("Electrode Pairs")(tunnel_allregion{dataInfo_idx}.Subregion==tunnel_reg));
    tunnels=tunnels(idx);
    
    for tun=1:num_tunnels
        ax=nexttile([1,maxTiles]);
        axes=[axes,ax];
        hold on
%         load(strcat(raw_dir,'\',dataInfo.meaName(dataInfo_idx),'_mat_files','\',tunnels(tun),".mat"))
        load(strcat(raw_dir,'\',dataInfo.meaName(dataInfo_idx),'\',tunnels(tun),".mat"))
        plot(time_idx,data,'k')
        
        files=string(ls(strcat(spikes_dir,'\',dataInfo.meaName(dataInfo_idx),'_mat_files','\')));
        to_load=files(contains(files,tunnels(tun)+"_"));
        color_idx=1;
        for i=1:length(to_load)
%             load(strcat(spikes_dir,'\',dataInfo.meaName(dataInfo_idx),'_mat_files','\',to_load(i)))
            load(strcat(spikes_dir,'\',dataInfo.meaName(dataInfo_idx),'\',to_load(i)))
            for j=1:length(max(cluster_class(:,1)))
                spike_idx=cluster_class(cluster_class(:,1)==j,2)/1000;
                y_spike=data(ismembertol(time_idx,spike_idx,'DataScale',1));
                plot(spike_idx,y_spike,'.','Color',wave_color(color_idx,:))
                color_idx=color_idx+1;
            end
        end
        
        title("Tunnel: " + tunnels(tun)+ " Reg: "+tunnel_reg)
        hold off
    end
    
    wells=well_allregion{dataInfo_idx}.Electrode(well_allregion{dataInfo_idx}.Subregion==well_reg & ...
        cellfun(is_pop,well_allregion{dataInfo_idx}.("Spike Train")));
    
    for well=1:num_well
        ax=nexttile([1,maxTiles]);
        axes=[axes,ax];
        hold on
        
        title("Well: " + wells(well) + " Reg: "+well_reg)
        
%         load(strcat(raw_dir,'\',dataInfo.meaName(dataInfo_idx),'_mat_files','\',wells(well),".mat"))
        load(strcat(raw_dir,'\',dataInfo.meaName(dataInfo_idx),'\',wells(well),".mat"))
        plot(time_idx,data,'k')
        
%         load(strcat(spikes_dir,'\',dataInfo.meaName(dataInfo_idx),'_mat_files','\',"times_",wells(well),"_well.mat"))
        load(strcat(spikes_dir,'\',dataInfo.meaName(dataInfo_idx),'\',"times_",wells(well),"_well.mat"))

        color_idx=1;
        for j=1:length(max(cluster_class(:,1)))
            spike_idx=cluster_class(cluster_class(:,1)==j,2)/1000;
            y_spike=data(ismembertol(time_idx,spike_idx,'DataScale',1));
            plot(spike_idx,y_spike,'.','Color',wave_color(color_idx,:))
            color_idx=color_idx+1;
        end
        
        hold off
    end
    
%     x_slide_scale(f,axes)
    
    linkaxes(axes,'x')
    
    to_del=0;
    while ~isempty(to_del)
        to_del=input("Plot to delete: ");
        if ~isempty(to_del)
            delete(nexttile(to_del))
        end
    end
    
    x=" ";
    more=" ";
    
    while ~(more=="no" || more=="No" || more=="")
        while ~(x=="yes" || x=="Yes" || x=="")
            try
                lo=input("Enter Lower Bounds: ");
                hi=input("Enter Upper Bounds: ");
                xlim([lo,hi])
                ylim("auto")
                to_del=0;
                while ~isempty(to_del)
                    to_del=input("Plot to delete: ");
                    if ~isempty(to_del)
                        delete(nexttile(to_del))
                    end
                end
                xlim([lo,hi])
                ylim("auto")
                x=input("Done? ","s");
                x=string(x);
            catch
                lo=0;
                hi=300;
                xlim([lo,hi])
                ylim("auto")
                x="";
            end
        end
        
        x=" ";
        saveas(gcf,save_to+tunnel_reg+" "+well_reg+" "+"lo " +string(lo)+ " hi " + string(hi)+".png")
        more=input("Save more?","s");
    end
    
end

%% Related functions

% function x_slide_scale(fig,axes)
% 
% % slider=figure('units','normalized','OuterPosition',[0.75 0.75 0.25 0.25]);
% p1=uipanel(fig,'Position',[0.75,0.75,0.1,0.1]);
% p2=uipanel(fig,'Position',[0.85,0.75,0.1,0.1]);
% linkaxes(axes,'x')
% low_x=uicontrol(p1,'Style','slider','String',"Lower Bounds");
% hi_x=uicontrol(p2,'Style','slider','String',"Upper Bounds");
% 
% end