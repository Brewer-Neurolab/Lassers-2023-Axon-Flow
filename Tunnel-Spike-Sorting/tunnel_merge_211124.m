function tunnel_merge_211124(folder_name,length_rec_sec)

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


%keep_times=[];
% Loop through channels
% Ask about this
for chani=2:no_channels
    bincountmatrix = [];
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
        down_var = load(char("times_"+downstream_channel_name), 'cluster_class','par','spikes');
    catch
        warning(char("Could not find "+downstream_channel_name))
        continue
    end

    try
        up_var = load(char("times_"+upstream_channel_name), 'cluster_class','par','spikes');
    catch
        warning(char("Could not find "+upstream_channel_name))
        continue
    end
   
    no_units = size(matching_table.unit_pairs{chani},1); 
    down_clusI=[];
    up_clusI=[];
    nmi_vals=[];
    %Loop through unit-pairs
    for ui=1:no_units
        axon_unit_pair = matching_table.unit_pairs{chani}(ui,:);
        nmi_vals(ui)=matching_table.NMI{chani}(ui);
        down_clusI(ui) = axon_unit_pair(1);
        up_clusI(ui) =  axon_unit_pair(2);
        %plot conduction time histogram
        down_channel = down_var.cluster_class(down_var.cluster_class(:,1) == down_clusI(ui), 2);
        up_channel = up_var.cluster_class(up_var.cluster_class(:,1) == up_clusI(ui), 2);
        %xbin_ISI should be the same every time, ok to leave it not
        %updating
        [bincountvec, ~, xbin_ISI,~,kept_times] = get_spike_delay_hist_210712(down_channel, up_channel, if_flip, if_hist, length_rec_sec);
        %saving
        bincountmatrix = [bincountmatrix;bincountvec];
        %keep_times=[keep_times;kept_times];
        title(strcat(char(matching_table.tunnel_pairs(chani)+"-u"+ui),' up ',string(up_clusI(ui)),' down ',string(down_clusI(ui))))
        saveas(gcf,char(matching_table.tunnel_pairs(chani)+"-u"+ui),'png')
        channel_Out=[channel_Out; {char(matching_table.tunnel_pairs(chani)+"-u"+ui)}];
    end
    disp(chani+" processed")
    
    peakprom_mult=0.12;
    std_num=1.9;
    
    max_pks=[];
    max_locs=[];
    to_delete=[];
    %find primary peak
    for ui=1:length(bincountmatrix(:,1))
        peakprom=peakprom_mult*max(bincountmatrix(ui,:));
        min_height=std_num*std(bincountmatrix(ui,:));
        peak_distance=abs(diff(xbin_ISI));
        peak_distance=peak_distance(1)*2;
        [pks,locs]=findpeaks(bincountmatrix(ui,:),'MinPeakProminence',peakprom,'MinPeakHeight',min_height,...
            'MinPeakDistance',peak_distance);
        if isempty(pks)
            to_delete=[to_delete,ui]
        end
        max_pks=max(pks);
        max_pks(ui)=max_pks(1);
        max_locs(ui)=locs(pks==max_pks(ui));
    end
    
    %delete no peaks
    nmi_vals(to_delete)=[];
    down_clusI(to_delete)=[];
    up_clusI(to_delete)=[];
    bincountmatrix(to_delete,:)=[];
    matching_table.unit_pairs{chani}(to_delete,:)=[];
    matching_table.NMI{chani}(to_delete)=[];
    
    %find similar peaks
    [sorted_locs,sorted_locs_I]=sort(max_locs);
    bins_apart=find(diff([0 sorted_locs])>1);
    p = diff([1 bins_apart numel(sorted_locs)+1]);
    to_merge_locs=mat2cell(sorted_locs_I,1,p);
    nmi_vals=mat2cell(nmi_vals(sorted_locs_I),1,p);
    down_clusI=mat2cell(down_clusI(sorted_locs_I),1,p);
    up_clusI=mat2cell(up_clusI(sorted_locs_I),1,p);
    
    %ensure clusters are not duplicated across different conduction times,
    %elimnate based on higher/lower NMI
    unique_down=[];
    unique_up=[];
    unique_down_length=[];
    unique_up_length=[];
    %check each conduction time cell for unique cluster classes, eliminate
    %those that have low NMI till only one remains
    for i=1:length(down_clusI)
        unique_down=[unique_down,unique(down_clusI{i},'stable')];
        unique_down_length(i)=length(unique(down_clusI{i},'stable'));
        unique_up=[unique_up,unique(up_clusI{i},'stable')];
        unique_up_length(i)=length(unique(up_clusI{i},'stable'));
    end
    down_counts=histcounts(unique_down);
    up_counts=histcounts(unique_up);
    
    multi_time_clusters_down=find(down_counts>=2);
    multi_time_clusters_up=find(up_counts>=2);
    %check downstream clusters for duplicates
    for i=1:length(multi_time_clusters_down)
        max_nmi_down=[];
        for j=1:length(down_clusI)
            max_nmi=max(nmi_vals{j}(down_clusI{j}==multi_time_clusters_down));
            if ~isempty(max_nmi)
                max_nmi=Nan;
            end
            max_nmi_down(j)=max_nmi;
        end
        [~,max_cdt_NMI]=max(max_nmi_down);
        down_clusI_idx=[1:length(down_clusI)];
        not_max_cdt_NMI=find(down_clusI_idx~=max_cdt_NMI);
        low_NMI_cdt=[down_clusI{not_max_cdt_NMI}];
        to_drop_clus=low_NMI_cdt==multi_time_clusters_down(i);
        to_drop_idx=to_merge_locs{not_max_cdt_NMI};
    end
    
    
    %check if cluster heights are the same
    down_std=[];
    up_std=[];
    down_mean=[];
    up_mean=[];
    for ui=1:length(bincountmatrix(:,1))
        down_spike_height=down_var.spikes(down_var.cluster_class(:,1) == down_clusI(ui), down_var.par.w_pre);
        up_spike_height=up_var.spikes(up_var.cluster_class(:,1) == up_clusI(ui), up_var.par.w_pre);
        down_std=[down_std;[std(down_spike_height),-std(down_spike_height)]];
        up_std=[up_std;[std(up_spike_height),-std(up_spike_height)]];
        down_mean=[down_mean;mean(down_spike_height)];
        up_mean=[up_mean;mean(up_spike_height)];
    end
    %organize down/up
    down_std=mat2cell(down_std(sorted_locs_I,:),1,p);
    up_std=mat2cell(up_std(sorted_locs_I,:),1,p);
    down_mean=mat2cell(down_mean(sorted_locs_I),1,p);
    up_mean=mat2cell(up_mean(sorted_locs_I),1,p);
    
    %allow merging within cells if mean falls within std
    %check all are the same length
    if_same_length=(length(down_std)==length(up_std))&&(length(down_std)==length(up_std))&&(length(up_mean)==length(down_mean));
    if ~if_same_length
        error("Cannot compare std and means.")
    end
    
    i=1;
    while i<=length(down_std)
        current_down_std=down_std{i};
        %current_up_std=up_std{i};
        current_down_mean=down_mean{i};
        current_down_clus=down_clus{i};
        merge_matrix=zeros(length(current_down_mean));
        %j is row values, representing different cluster means
        for j=1:length(current_down_mean)
            %k is column values, iterates through std vals
            for k=1:length(current_down_std)
                %add 1 if inside std and not the same square
                if current_down_mean(j)>current_down_std(j,1) && current_down_mean(j)<current_down_std(j,2) && j~=k
                    merge_matrix(j,k)=1;
                end
            end
        end
        
        for j=1:length(current_down_mean)
            for k=1:length(current_down_std)
                if merge_matrix(j,k) && j~=k
                    current_down_clus(k)=i;
                    current_down_clus(j)=i;
                end
            end
        end
        
    end
    
    
    
end
    
end